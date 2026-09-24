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

# Python modules
import os

# Third party modules
import numpy as np
import xarray as xr

# LS2D modules
import ls2d.ecmwf.era_tools as era_tools
import ls2d.ecmwf.htessel as htessel
from ls2d.ecmwf.IFS_tools import IFS_tools
from ls2d.column.spec import validate
from ls2d.core.messages import *

ifs = IFS_tools('L137')

# Molar mass ratio dry air / ozone.
_md_mo3 = 28.9644 / 47.9982


def read_era5_arco(settings):
    """
    Read ERA5 files downloaded by `download_era5_arco()`, and
    return them as a generic 3D dataset (see `ls2d.column.spec`).

    Arguments:
        settings : dict
            Dictionary with keys `central_lat`, `central_lon`, `era5_path`,
            `case_name`, `start_date`, and `end_date`.
    """

    start = era_tools.lower_to_hour(settings['start_date'])
    end = era_tools.lower_to_hour(settings['end_date'])

    header(f'Reading ERA5 (ARCO) from {start} to {end}')

    files = [
        era_tools.era5_file_path(d.year, d.month, d.day, settings['era5_path'], settings['case_name'], 'era5_arco', False)
        for d in era_tools.get_required_analysis(start, end)
    ]
    for f in files:
        if not os.path.exists(f):
            error(f'File "{f}" does not exist. Run `ls2d.download_era5_arco()` first.')

    era = xr.concat([xr.open_dataset(f) for f in files], dim='time').sel(time=slice(start, end)).load()

    # Levels surface->top, latitude south->north, pressure levels surface->top.
    era = era.isel(level=slice(None, None, -1), latitude=slice(None, None, -1), pressure_level=slice(None, None, -1))

    T = era.t.values
    qv = era.q.values
    ql = era.clwc.values + era.ciwc.values + era.crwc.values + era.cswc.values
    qt = qv + ql
    Tv = ifs.calc_virtual_temp(T, qv, era.clwc.values, era.ciwc.values, era.crwc.values, era.cswc.values)

    # Half level pressure and height from the L137 hybrid coefficients.
    ps = era.sp.values
    ph = ifs.a[None, :, None, None] + ifs.b[None, :, None, None] * ps[:, None]
    ph[:, -1] = 0.34  # As in `IFS_tools.calc_half_level_pressure()`

    dzh = -ifs.Rd * Tv * np.log(ph[:, 1:] / ph[:, :-1]) / ifs.grav
    zh = np.concatenate((np.zeros_like(ps)[:, None], np.cumsum(dzh, axis=1)), axis=1)

    p = 0.5 * (ph[:, 1:] + ph[:, :-1])
    z = 0.5 * (zh[:, 1:] + zh[:, :-1])

    exn = ifs.calc_exner(p)
    thl = T / exn - ifs.Lv / (ifs.cpd * exn) * ql
    rho = p / (ifs.Rd * Tv)

    # Surface fluxes (ERA5 = positive downward) to upward kinematic fluxes.
    rhos = ps / (ifs.Rd * ifs.calc_virtual_temp(era.skt.values, qv[:, 0]))
    wth = -era.ishf.values / (rhos * ifs.cpd * ifs.calc_exner(ps))
    wq = -era.ie.values / rhos

    # Land surface. Vegetation/soil types and root fractions are time invariant.
    types = {name: np.round(era[name].values).astype(np.int32) for name in ('slt', 'tvl', 'tvh')}
    is_sea = types['slt'][0] == 0
    root_frac = {
        name: np.where(is_sea, -1.0, htessel.root_fraction(np.maximum(types[name][0], 1))) for name in ('tvl', 'tvh')
    }

    d3 = ('time', 'level', 'latitude', 'longitude')
    d3h = ('time', 'level_half', 'latitude', 'longitude')
    d3s = ('time', 'soil_layer', 'latitude', 'longitude')
    d2 = ('time', 'latitude', 'longitude')

    ds = xr.Dataset(
        {
            'p': (d3, p),
            'z': (d3, z),
            'T': (d3, T),
            'thl': (d3, thl),
            'qt': (d3, qt),
            'u': (d3, era.u.values),
            'v': (d3, era.v.values),
            'w': (d3, -era.w.values / (rho * ifs.grav)),
            'h2o': (d3, qv / ((ifs.Rd / ifs.Rv) * (1 - qt))),
            'o3': (d3, era.o3.values * _md_mo3),
            'ph': (d3h, ph),
            'zh': (d3h, zh),
            'phi_p': (('time', 'pressure_level', 'latitude', 'longitude'), era.z.values),
            'ps': (d2, ps),
            'ts': (d2, era.skt.values),
            'sst': (d2, era.sst.values),
            'wth': (d2, wth),
            'wq': (d2, wq),
            'z0m': (d2, era.fsr.values),
            'z0h': (d2, np.exp(era.flsr.values)),
            'lai_low': (d2, era.lai_lv.values),
            'lai_high': (d2, era.lai_hv.values),
            'c_low_veg': (d2, era.cvl.values),
            'c_high_veg': (d2, era.cvh.values),
            'type_soil': (d2, types['slt']),
            'type_low_veg': (d2, types['tvl']),
            'type_high_veg': (d2, types['tvh']),
            't_soil': (d3s, np.stack([era[f'stl{i}'].values for i in range(1, 5)], axis=1)),
            'theta_soil': (d3s, np.stack([era[f'swvl{i}'].values for i in range(1, 5)], axis=1)),
            'root_frac_low_veg': (d3s[1:], root_frac['tvl']),
            'root_frac_high_veg': (d3s[1:], root_frac['tvh']),
        },
        coords={
            'time': era.time.values,
            'pressure_level': era.pressure_level.values * 100.0,
            'latitude': era.latitude.values,
            'longitude': era.longitude.values,
            'z_soil': ('soil_layer', htessel.z_soil),
        },
        attrs={
            'central_lat': settings['central_lat'],
            'central_lon': settings['central_lon'],
            'source': 'ERA5 (Google ARCO)',
        },
    )

    validate(ds)

    return ds
