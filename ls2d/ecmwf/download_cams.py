#
# This file is part of LS2D.
#
# Copyright (c) 2017-2024 Wageningen University & Research
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
import subprocess as sp
import datetime
import sys
import os
import pickle
import requests
import yaml

# Third party modules
import xarray as xr
import numpy as np

# LS2D modules
import ls2d.ecmwf.era_tools as era_tools
from ls2d.src.messages import *

# Yikes, but necessary (?) if you want to use
# MARS downloads without the Python CDS api installed?
try:
    import cdsapi
except ImportError:
    cdsapi = None


def regrid(nc_file, central_lon, central_lat, resolution):
    """
    Re-grid regular lat/lon NetCDF file (like ERA5 or CAMS)
    onto new grid / resolution. Note: this function
    overwrites the original `nc_file`.

    Arguments:
        nc_file : str
            Path to NetCDF file.
        central_lon : float
            Central longitude of new grid
        central_lat : float
            Central latitude of new grid
        resolution : float
            New output resolution
    """

    ds = xr.open_dataset(nc_file)

    # Find start/end lon/lat to stay within bounds of input dataset.
    lon_frac = np.abs((ds.longitude.values - central_lon) / resolution)
    lat_frac = np.abs((ds.latitude.values  - central_lat) / resolution)

    lon_0 = central_lon - np.floor(lon_frac[0 ]) * resolution
    lon_1 = central_lon + np.floor(lon_frac[-1]) * resolution

    lat_0 = central_lat - np.floor(lat_frac[0 ]) * resolution
    lat_1 = central_lat + np.floor(lat_frac[-1]) * resolution

    lon_out = np.arange(lon_0, lon_1+1e-12, resolution)
    lat_out = np.arange(lat_0, lat_1+1e-12, resolution)

    # Interpolate. Extrapolation is necessary if the first or last
    # new coordinate is exactly equal to the start/end coordinate
    # of the input NetCDF file. Otherwise, those values are set to
    # NaN (no idea why....).
    dsi = ds.interp(longitude=lon_out, latitude=lat_out, kwargs={'fill_value': 'extrapolate'})

    # Save back.
    ds.close()
    dsi.to_netcdf(nc_file)


def _download_cams_file(settings, variables, grid):
    """
    Download single CAMS file.
    """
    header('Downloading: {} - {}'.format(settings['date'], settings['ftype']))

    variables = variables[settings['ftype']]

    # Read ADS url/key from `.cdsapirc` file.
    with open(settings['cdsapirc'], 'r') as f:
        credentials = yaml.safe_load(f)

    if 'url_ads' not in credentials.keys() or 'key_ads' not in credentials:
        error('You need to specify `url_ads` and `key_ads` in your `.cdsapirc` file!')

    # Keep track of CDS downloads which are finished:
    finished = False

    # Output file name
    nc_dir, nc_file = era_tools.era5_file_path(
            settings['date'].year, settings['date'].month, settings['date'].day,
            settings['cams_path'], settings['case_name'], settings['ftype'])

    # Write CDS API prints to log file (NetCDF file path/name appended with .out/.err)
    if settings['write_log']:
        out_file   = '{}.out'.format(nc_file[:-3])
        err_file   = '{}.err'.format(nc_file[:-3])
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = open(out_file, 'w')
        sys.stderr = open(err_file, 'w')

    # Bounds of domain
    lat_n = settings['central_lat']+settings['area_size']
    lat_s = settings['central_lat']-settings['area_size']
    lon_w = settings['central_lon']-settings['area_size']
    lon_e = settings['central_lon']+settings['area_size']

    # Check if pickle with previous request is available.
    # If so, try to download NetCDF file, if not, submit new request
    pickle_file = '{}.pickle'.format(nc_file[:-3])

    if os.path.isfile(pickle_file):
        message('Found previous CDS request!')

        with open(pickle_file, 'rb') as f:
            cds_request = pickle.load(f)

            try:
                cds_request.update()
            except requests.exceptions.HTTPError:
                error('CDS request is no longer available online!', exit=False)
                error('To continue, delete the previous request: {}'.format(pickle_file))

            state = cds_request.reply['state']

            if state == 'completed':
                message('Request finished, downloading NetCDF file')

                cds_request.download(nc_file)
                os.remove(pickle_file)

                if grid is not None:
                    message(f'Re-gridding NetCDF to {grid:.2f}°×{grid:.2f}° degree grid.')
                    regrid(nc_file, settings['central_lon'], settings['central_lat'], grid)

                finished = True

            elif state in ('queued', 'running'):
                message('Request not finished, current status = \"{}\"'.format(state))

            else:
                error('Request failed, status = \"{}\"'.format(state), exit=False)
                message('Error message = {}'.format(cds_request.reply['error'].get('message')))
                message('Error reason = {}'.format(cds_request.reply['error'].get('reason')))


    else:
        message('No previous CDS request, submitting new one')

        # Create instance of CDS API
        server = cdsapi.Client(
                url=credentials['url_ads'], key=credentials['key_ads'], verify=True,
                wait_until_complete=False, delete=False)

        model_level = [str(x) for x in range(1,61)]
        analysis_times = ['{0:02d}:00'.format(i) for i in range(0, 22, 3)]
        steps = ['{}'.format(i) for i in range(0, 22, 3)]
        area = [lat_n, lon_w, lat_s, lon_e]
        date = settings['date'].strftime('%Y-%m-%d')

        request = {
            'format': 'netcdf',
            'variable': variables,
            'date': date,
            'area': area,
            #'grid': '0.25/0.25'   # NOTE to self: not supported in the ADS!
        }

        if settings['ftype'] == 'eac4_ml' or settings['ftype'] == 'eac4_sfc':
            request.update({'time': analysis_times})
        elif settings['ftype'] == 'egg4_ml':
            request.update({'step': steps})

        if settings['ftype'] == 'eac4_ml' or settings['ftype'] == 'egg4_ml':
            request.update({'model_level': model_level})

        if settings['ftype'] == 'eac4_ml' or settings['ftype'] == 'eac4_sfc':
            cds_request = server.retrieve('cams-global-reanalysis-eac4', request)
        elif settings['ftype'] == 'egg4_ml':
            cds_request = server.retrieve('cams-global-ghg-reanalysis-egg4', request)

        # Save pickle for later processing/download
        with open(pickle_file, 'wb') as f:
            pickle.dump(cds_request, f)

    if settings['write_log']:
        sys.stdout = old_stdout
        sys.stderr = old_stderr

    return finished


def download_cams(settings, variables, grid=None):
    """
    Download all required CAMS fields.
    Only `data_source = 'CDS' is supported for now!

    Arguments:
        settings : dictionary
            Dictionary with keys:
                date : datetime object with date to download (always retrieves 00:00 to 23:00 UTC)
                lat, lon : requested latitude and longitude
                size : download an area of lat+/-size, lon+/-size (degrees)
                cams_path : absolute or relative path to save the NetCDF data
                case : case name used in file name of NetCDF files
                cdsapirc : absolute path to `.cdsapirc` file.
        variables : dictionary
            Nested dictionary with variables for each CAMS type (see examples).
        grid : float, optional (default = None)
            Option to regrid from 0.75 degree regular lat/lon grid (CAMS default) to another resolution.
    """

    # Checks!
    if settings['data_source'] != 'CDS':
        error('CAMS downloads only support CDS (for now...)!')

    if 'cdsapirc' not in settings.keys():
        error('You need to specify the location of your `.cdsapirc` file in `settings`!')

    header('Downloading CAMS for period: {} to {}'.format(settings['start_date'], settings['end_date']))

    # Check if output directory exists, and ends with '/'
    if not os.path.isdir(settings['cams_path']):
        error('Output directory \"{}\" does not exist!'.format(settings['cams_path']))
    if settings['cams_path'][-1] != '/':
        settings['cams_path'] += '/'

    if cdsapi is None:
        error('CDS API is not installed. See: https://cds.climate.copernicus.eu/api-how-to')

    # Round date/time to full hours
    start = era_tools.lower_to_hour(settings['start_date'])
    end   = era_tools.lower_to_hour(settings['end_date']  )

    # Get list of required forecast and analysis times
    an_dates = era_tools.get_required_analysis(start, end, freq=3)

    # Loop over all required files, check if there is a local version, if not add to download queue
    download_settings = settings.copy()
    download_queue = []

    for date in an_dates:
        for ftype in variables.keys():

            if ftype == 'egg4_ml':
                error('Downloading CAMS EGG4 data currently does not work due to an open ADS bug.')

            era_dir, era_file = era_tools.era5_file_path(
                    date.year, date.month, date.day, settings['cams_path'], settings['case_name'], ftype)

            if not os.path.exists(era_dir):
                message('Creating output directory {}'.format(era_dir))
                os.makedirs(era_dir)

            if os.path.isfile(era_file):
                message('Found {} - {} local'.format(date, ftype))
            else:
                settings_tmp = download_settings.copy()
                settings_tmp.update({'date': date, 'ftype': ftype})
                download_queue.append(settings_tmp)

    finished = True
    for req in download_queue:
        if not _download_cams_file(req, variables, grid):
            finished = False

    if not finished:
        print(' --------------------------------------------------------------')
        print(' | One or more requests are not finished.                     |')
        print(' | For ADS request, you can monitor the progress at:          |')
        print(' | https://ads.atmosphere.copernicus.eu/cdsapp#!/yourrequests |')
        print(' | This script will stop now, you can restart it              |')
        print(' | at any time to retry, or download the results.             |')
        print(' --------------------------------------------------------------')

        sys.exit(1)
