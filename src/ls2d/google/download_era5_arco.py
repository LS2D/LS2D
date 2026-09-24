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
import datetime
import time
import os

# Third party modules
import numpy as np
import xarray as xr
import requests
import gcsfs

# LS2D modules
import ls2d.ecmwf.era_tools as era_tools
from ls2d.google.arco_tools import get_layout, read_rows
from ls2d.core.messages import *

_bucket = 'gcp-public-data-arco-era5/ar'
_store_ml = f'{_bucket}/model-level-1h-0p25deg.zarr-v1'
_store_sl = f'{_bucket}/full_37-1h-0p25deg-chunk-1.zarr-v3'

# ARCO variable name -> output (ECMWF short) name.
# Model level analysis:
_vars_ml = {
    'temperature': 't',
    'u_component_of_wind': 'u',
    'v_component_of_wind': 'v',
    'vertical_velocity': 'w',
    'specific_humidity': 'q',
    'specific_cloud_liquid_water_content': 'clwc',
    'specific_cloud_ice_water_content': 'ciwc',
    'specific_rain_water_content': 'crwc',
    'specific_snow_water_content': 'cswc',
    'ozone_mass_mixing_ratio': 'o3',
}

# Pressure level analysis:
_vars_pl = {
    'geopotential': 'z',
}

# Surface analysis:
_vars_sfc = {
    'surface_pressure': 'sp',
    'skin_temperature': 'skt',
    'sea_surface_temperature': 'sst',
    'instantaneous_surface_sensible_heat_flux': 'ishf',
    'instantaneous_moisture_flux': 'ie',
    'forecast_surface_roughness': 'fsr',
    'forecast_logarithm_of_surface_roughness_for_heat': 'flsr',
    'soil_type': 'slt',
    'type_of_low_vegetation': 'tvl',
    'type_of_high_vegetation': 'tvh',
    'leaf_area_index_low_vegetation': 'lai_lv',
    'leaf_area_index_high_vegetation': 'lai_hv',
    'low_vegetation_cover': 'cvl',
    'high_vegetation_cover': 'cvh',
    'soil_temperature_level_1': 'stl1',
    'soil_temperature_level_2': 'stl2',
    'soil_temperature_level_3': 'stl3',
    'soil_temperature_level_4': 'stl4',
    'volumetric_soil_water_layer_1': 'swvl1',
    'volumetric_soil_water_layer_2': 'swvl2',
    'volumetric_soil_water_layer_3': 'swvl3',
    'volumetric_soil_water_layer_4': 'swvl4',
}


def _open_metadata(store):
    """
    Open ARCO store for coordinates and attributes only; no data is read through xarray/zarr.
    """
    return xr.open_zarr(f'gs://{store}', chunks=None, storage_options=dict(token='anon'))


def _last_valid_date(store):
    """
    Last day with valid data in ARCO store, including ERA5T if available.
    Uses plain HTTPS, so that we can exit cleanly before opening any `gcsfs` sessions.
    """
    attrs = requests.get(f'https://storage.googleapis.com/{store}/.zattrs', timeout=30).json()
    stop = attrs.get('valid_time_stop_era5t', attrs['valid_time_stop'])
    return datetime.datetime.strptime(stop, '%Y-%m-%d')


def _clean_attrs(attrs):
    """
    Remove GRIB attributes, which are not relevant after the ARCO conversion.
    """
    return {k: v for k, v in attrs.items() if not k.startswith('GRIB')}


def download_era5_arco(settings, batch_size=512):
    """
    Download all required ERA5 fields for an experiment between
    `start_date` and `end_date` from the Google ARCO-ERA5 archive.

    All model level, pressure level, and surface fields are saved
    as one NetCDF file per day (00 UTC to (including) 23 UTC), on the native
    0.25 deg ERA5 grid, as: `era5_path/case_name/yyyy/mm/dd/era5_arco.nc`.

    Arguments:
        settings : dict
            Dictionary with keys:
                central_lat, central_lon : requested latitude and longitude
                area_size : download an area of lat+/-size, lon+/-size (degrees)
                era5_path : absolute or relative path to save the NetCDF data
                case_name : case name used in file path of NetCDF files
                start_date, end_date : datetime objects with start/end of experiment
        batch_size : int
            Number of concurrent HTTP range requests.
    """

    header(f'Downloading ERA5 from Google ARCO for period: {settings["start_date"]} to {settings["end_date"]}')

    # Check if output directory exists.
    if not os.path.isdir(settings['era5_path']):
        error(f'Output directory "{settings["era5_path"]}" does not exist!')

    # Round date/time to full hours, and get list of days to download.
    start = era_tools.lower_to_hour(settings['start_date'])
    end = era_tools.lower_to_hour(settings['end_date'])
    an_dates = era_tools.get_required_analysis(start, end)

    download_dates = []
    for date in an_dates:
        era_dir, era_file = era_tools.era5_file_path(
            date.year, date.month, date.day, settings['era5_path'], settings['case_name'], 'era5_arco'
        )
        if os.path.isfile(era_file):
            message(f'Found {era_file} local')
        else:
            download_dates.append((date, era_dir, era_file))

    if len(download_dates) == 0:
        return True

    # Check if data is available for all requested days.
    last_valid = min(_last_valid_date(_store_ml), _last_valid_date(_store_sl))
    for date, _, _ in download_dates:
        if date > last_valid:
            error(f'ERA5 data at {date:%Y-%m-%d} is not (yet) available in ARCO. Last available day: {last_valid:%Y-%m-%d}')

    fs = gcsfs.GCSFileSystem(token='anon')
    ds_ml = _open_metadata(_store_ml)
    ds_sl = _open_metadata(_store_sl)

    # Both stores share the same lat/lon grid and time axis.
    lats = ds_sl.latitude.values
    lons = ds_sl.longitude.values
    if not (np.array_equal(lats, ds_ml.latitude.values) and np.array_equal(lons, ds_ml.longitude.values)):
        error('Model and single/pressure level ARCO stores have different grids!')

    # Select box on native grid, including one extra grid point for gradients on pressure levels.
    # Latitude = contiguous rows. Longitude: full rows are read anyway, so any
    # (also 0-deg crossing) selection is free. Output longitudes are -180..180, west->east,
    # except for boxes crossing the 180 deg meridian, which keep 0..360.
    half = settings['area_size'] + 0.25 + 1e-6
    ilat = np.where(np.abs(lats - settings['central_lat']) <= half)[0]
    ilat0, nlat = ilat.min(), ilat.size

    dlon = ((lons - settings['central_lon'] % 360 + 180) % 360) - 180
    ilon = np.where(np.abs(dlon) <= half)[0]
    ilon = ilon[np.argsort(dlon[ilon])]
    lons_out = ((lons[ilon] + 180) % 360) - 180
    if np.any(np.diff(lons_out) < 0):
        lons_out = lons[ilon]

    # Variables to read, with their chunk layout.
    fields = [(_store_ml, name, get_layout(ds_ml, name)) for name in _vars_ml]
    fields += [(_store_sl, name, get_layout(ds_sl, name)) for name in _vars_pl]
    fields += [(_store_sl, name, get_layout(ds_sl, name)) for name in _vars_sfc]

    arco_ds = {**{name: ds_ml for name in _vars_ml}, **{name: ds_sl for name in {**_vars_pl, **_vars_sfc}}}
    short_names = {**_vars_ml, **_vars_pl, **_vars_sfc}

    for date, era_dir, era_file in download_dates:
        header(f'Downloading {date:%Y-%m-%d}')

        if not os.path.exists(era_dir):
            message(f'Creating output directory {era_dir}')
            os.makedirs(era_dir)

        times = np.array([np.datetime64(date + datetime.timedelta(hours=h), 'ns') for h in range(24)])
        itimes = np.searchsorted(ds_sl.time.values, times)
        if not (np.array_equal(ds_sl.time.values[itimes], times) and np.array_equal(ds_ml.time.values[itimes], times)):
            error('Requested times not found in ARCO time axis!')

        data = {name: np.empty((24, lay['nlev'], nlat, ilon.size), np.float32) for _, name, lay in fields}

        t_start = time.perf_counter()
        for n, it in enumerate(itimes):
            rows = read_rows(fs, fields, it, ilat0, nlat, lats.size, lons.size, batch_size)
            for name in data:
                data[name][n] = rows[name][:, :, ilon]

        # Create combined dataset and save to NetCDF.
        variables = {}
        for _, name, _ in fields:
            da = arco_ds[name][name]
            attrs = _clean_attrs(da.attrs)
            if name in _vars_ml:
                variables[short_names[name]] = (('time', 'level', 'latitude', 'longitude'), data[name], attrs)
            elif name in _vars_pl:
                variables[short_names[name]] = (('time', 'pressure_level', 'latitude', 'longitude'), data[name], attrs)
            else:
                variables[short_names[name]] = (('time', 'latitude', 'longitude'), data[name][:, 0], attrs)

        ds = xr.Dataset(
            variables,
            coords=dict(
                time=times,
                level=('level', ds_ml.hybrid.values.astype(np.int32), dict(long_name='model level')),
                pressure_level=('pressure_level', ds_sl.level.values.astype(np.int32), ds_sl.level.attrs),
                latitude=('latitude', lats[ilat0 : ilat0 + nlat], ds_sl.latitude.attrs),
                longitude=('longitude', lons_out, ds_sl.longitude.attrs),
            ),
            attrs=dict(
                source='ERA5 from Google ARCO-ERA5 (https://github.com/google-research/arco-era5)',
                stores=f'gs://{_store_ml}, gs://{_store_sl}',
                history=f'Downloaded by (LS)2D on {datetime.datetime.now():%Y-%m-%d %H:%M}',
            ),
        )

        # Write to temporary file first, so that an interrupted download does not leave a valid looking file.
        tmp_file = f'{era_file}.tmp'
        ds.to_netcdf(tmp_file)
        os.replace(tmp_file, era_file)
        message(f'Saved {era_file} in {time.perf_counter() - t_start:.0f} sec')

    return True
