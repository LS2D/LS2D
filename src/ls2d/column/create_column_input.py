#
# This file is part of LS2D.
#
# Copyright (c) 2017-2026 Wageningen University & Research
# Author: Bart van Stratum (WUR)
#
# LS2D is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LS2D is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with LS2D.  If not, see <http://www.gnu.org/licenses/>.
#

# Third party modules
import numpy as np
import xarray as xr

# LS2D modules
import ls2d.core.spatial_tools as spatial
from ls2d.column.spec import validate
from ls2d.core.messages import *

r_earth = 6.37e6
omega_earth = 7.2921e-5


def _interp_extrap(x, xp, fp):
    """
    Linear interpolation, with linear extrapolation outside `xp`.
    """
    i = np.argsort(xp)
    xp, fp = xp[i], fp[i]
    y = np.interp(x, xp, fp)
    lo, hi = x < xp[0], x > xp[-1]
    y[lo] = fp[0] + (x[lo] - xp[0]) * (fp[1] - fp[0]) / (xp[1] - xp[0])
    y[hi] = fp[-1] + (x[hi] - xp[-1]) * (fp[-1] - fp[-2]) / (xp[-1] - xp[-2])
    return y


def create_column_input(ds, z, n_av=0):
    """
    Calculate the large-scale forcings and column input for LES/SCM
    from a generic 3D dataset (see `ls2d.column.spec`).

    Arguments:
        ds : xarray.Dataset
            Generic 3D dataset, e.g. from `ls2d.read_era5_arco()`.
        z : np.ndarray
            Full level heights LES grid (m).
        n_av : int
            Number of grid points (+/-) over which the fields and forcings are averaged.

    Returns:
        xarray.Dataset with the LES input.
    """

    header('Calculating large-scale forcings')

    validate(ds)

    clat = ds.attrs['central_lat']
    clon = ds.attrs['central_lon']

    # Nearest column.
    i = int(np.abs(ds.longitude.values - clon).argmin())
    j = int(np.abs(ds.latitude.values - clat).argmin())

    lat, lon = float(ds.latitude[j]), float(ds.longitude[i])
    distance = spatial.haversine(lon, lat, clon, clat)
    message(f'Using nearest lat/lon = {lat:.2f}/{lon:.2f} (requested = {clat:.2f}/{clon:.2f}), distance ~= {distance / 1000:.1f} km')

    # Averaging window, plus one grid point on each side for the gradients.
    if min(i, j) - n_av < 1 or i + n_av > ds.sizes['longitude'] - 2 or j + n_av > ds.sizes['latitude'] - 2:
        error(f'Domain too small for n_av={n_av}; download a larger area.')

    sub = ds.isel(latitude=slice(j - n_av - 1, j + n_av + 2), longitude=slice(i - n_av - 1, i + n_av + 2))
    inner = dict(latitude=slice(1, -1), longitude=slice(1, -1))

    dlon = float(ds.longitude[1] - ds.longitude[0])
    dlat = float(ds.latitude[1] - ds.latitude[0])
    area = f'{(1 + 2 * n_av) * dlon:.2f}°×{(1 + 2 * n_av) * dlat:.2f}°'
    message(f'Averaging over a {area} spatial area.')

    def mean(da):
        return da.isel(inner).mean(('latitude', 'longitude'))

    # Horizontal derivatives (2nd order), from per-degree to per-meter.
    m_per_deg = r_earth * np.pi / 180

    def ddx(da):
        return da.differentiate('longitude') / (m_per_deg * np.cos(np.deg2rad(sub.latitude)))

    def ddy(da):
        return da.differentiate('latitude') / m_per_deg

    def advec(da):
        return mean(-sub.u * ddx(da) - sub.v * ddy(da))

    # Geostrophic wind on pressure levels, interpolated to the mean model level pressure.
    fc = 2 * omega_earth * np.sin(np.deg2rad(clat))
    col = mean(sub)

    def to_p(da):
        return xr.apply_ufunc(
            _interp_extrap, col.p, col.pressure_level, da,
            input_core_dims=[['level'], ['pressure_level'], ['pressure_level']],
            output_core_dims=[['level']], vectorize=True,
        )

    ug = to_p(mean(-ddy(sub.phi_p) / fc))
    vg = to_p(mean(ddx(sub.phi_p) / fc))

    # Half level temperature for radiation, extrapolated at the surface and top.
    T, Tz, zh = col['T'].values, col.z.values, col.zh.values
    Th = np.empty_like(zh)
    Th[:, 1:-1] = 0.5 * (T[:, 1:] + T[:, :-1])
    Th[:, 0] = T[:, 0] - (Th[:, 1] - T[:, 0]) / (zh[:, 1] - Tz[:, 0]) * Tz[:, 0]
    Th[:, -1] = T[:, -1] + (T[:, -1] - Th[:, -2]) / (Tz[:, -1] - zh[:, -2]) * (zh[:, -1] - Tz[:, -1])

    #
    # Output, interpolated to LES grid.
    #
    z = xr.DataArray(z, dims='z', coords=dict(z=z), attrs=dict(long_name='full level height LES', units='m'))

    def to_les(da):
        return xr.apply_ufunc(
            np.interp, z, col.z, da,
            input_core_dims=[['z'], ['level'], ['level']], output_core_dims=[['z']], vectorize=True,
        )

    lay = dict(level='lay')
    lev = dict(level_half='lev')

    variables = {
        'thl': (to_les(col.thl), 'liquid water potential temperature', 'K'),
        'qt': (to_les(col.qt), 'total specific humidity', 'kg kg-1'),
        'u': (to_les(col.u), 'zonal wind component', 'm s-1'),
        'v': (to_les(col.v), 'meridional wind component', 'm s-1'),
        'wls': (to_les(col.w), 'vertical wind component', 'm s-1'),
        'p': (to_les(col.p), 'air pressure', 'Pa'),
        'dtthl_advec': (to_les(advec(sub.thl)), 'advective tendency liquid water potential temperature', 'K s-1'),
        'dtqt_advec': (to_les(advec(sub.qt)), 'advective tendency total specific humidity', 'kg kg-1 s-1'),
        'dtu_advec': (to_les(advec(sub.u)), 'advective tendency zonal wind', 'm s-2'),
        'dtv_advec': (to_les(advec(sub.v)), 'advective tendency meridional wind', 'm s-2'),
        'ug': (to_les(ug), 'geostrophic wind component zonal wind', 'm s-1'),
        'vg': (to_les(vg), 'geostrophic wind component meridional wind', 'm s-1'),
        'z_lay': (col.z.rename(lay), 'Full level heights radiation', 'm'),
        'z_lev': (col.zh.rename(lev), 'Half level heights radiation', 'm'),
        'p_lay': (col.p.rename(lay), 'full level pressure radiation', 'Pa'),
        'p_lev': (col.ph.rename(lev), 'half level pressure radiation', 'Pa'),
        't_lay': (col['T'].rename(lay), 'full level temperature radiation', 'K'),
        't_lev': (xr.DataArray(Th, dims=('time', 'lev')), 'half level temperature radiation', 'K'),
        'h2o_lay': (col.h2o.rename(lay), 'moisture volume mixing ratio', ''),
        'ps': (col.ps, 'surface pressure', 'Pa'),
        'ts': (col.ts, 'surface (skin) temperature', 'K'),
        'wth': (col.wth, 'surface sensible heat flux', 'K m s-1'),
        'wq': (col.wq, 'surface latent heat flux', 'kg kg-1 m s-1'),
    }

    # Optional variables.
    if 'o3' in col:
        variables['o3'] = (to_les(col.o3) * 1e6, 'ozone volume mixing ratio', 'ppmv')
        variables['o3_lay'] = (col.o3.rename(lay) * 1e6, 'ozone volume mixing ratio radiation', 'ppmv')

    optional = {
        'sst': ('sst', 'sea surface temperature', 'K'),
        'lai_low_veg': ('lai_low', 'LAI low vegetation', '-'),
        'lai_high_veg': ('lai_high', 'LAI high vegetation', '-'),
        'c_low_veg': ('c_low_veg', 'fraction low vegetation', '-'),
        'c_high_veg': ('c_high_veg', 'fraction high vegetation', '-'),
        'z0m': ('z0m', 'roughness length momentum', 'm'),
        'z0h': ('z0h', 'roughness length scalars', 'm'),
        't_soil': ('t_soil', 'soil temperature', 'K'),
        'theta_soil': ('theta_soil', 'soil moisture content', 'm3 m-3'),
    }
    for name, (src, long_name, units) in optional.items():
        if src in col:
            variables[name] = (col[src], long_name, units)

    # Nearest-neighbour land surface types and root fractions.
    nn = ds.isel(latitude=j, longitude=i)
    if 'type_soil' in ds:
        is_land = int(nn.type_soil[0]) != 0
        if is_land:
            message('Selected grid point is over land.')
        else:
            warning('Selected grid point is water/sea! Setting vegetation/soil indexes to 1e9.')

        for name, long_name in [
            ('type_soil', 'soil type (Fortran indexing!)'),
            ('type_low_veg', 'low vegetation type (Fortran indexing!)'),
            ('type_high_veg', 'high vegetation type (Fortran indexing!)'),
        ]:
            value = int(nn[name][0]) if is_land else int(1e9)
            variables[name] = (xr.DataArray(value), long_name, '-')

    for name, long_name in [
        ('root_frac_low_veg', 'root fraction low vegetation'),
        ('root_frac_high_veg', 'root fraction high vegetation'),
    ]:
        if name in ds:
            variables[name] = (nn[name], long_name, '-')

    out = xr.Dataset()
    for name, (da, long_name, units) in variables.items():
        out[name] = da.rename(soil_layer='zs') if 'soil_layer' in da.dims else da
        out[name].attrs = dict(long_name=long_name, units=units)

    out = out.drop_vars([c for c in out.coords if c not in out.dims])
    if 'zs' in out.dims:
        out = out.assign_coords(zs=('zs', ds.z_soil.values, dict(long_name='full level depth soil', units='m')))
    out = out.assign_coords(lay=np.arange(out.sizes['lay']), lev=np.arange(out.sizes['lev']))
    out['z'].attrs = z.attrs

    out['time_sec'] = ('time', (out.time.values - out.time.values[0]) / np.timedelta64(1, 's'))
    out['time_sec'].attrs = dict(long_name='seconds since start of experiment', units='s')

    out.attrs = {
        'fc': fc,
        'central_lon': clon,
        'central_lat': clat,
        'area': f'{area} spatial average',
        'source': f'{ds.attrs.get("source", "unknown")} + (LS)²D',
        'description': 'Generated by (LS)²D: https://github.com/LS2D & https://pypi.org/project/ls2d)',
        'reference': (
            'van Stratum et al. (2023). The benefits and challenges of downscaling a global reanalysis with '
            'doubly-periodic large-eddy simulations. JAMES, https://doi.org/10.1029/2023MS003750'
        ),
    }

    return out
