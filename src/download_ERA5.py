#
# This file is part of LS2D.
#
# Copyright (c) 2017-2018 Bart van Stratum
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

# Standard Python packages
import sys
import os
import datetime
import numpy as np

import subprocess as sp
import multiprocessing as mp

# Custom tools (in src subdirectory)
import time_tools as tt
from messages import *


def retrieve_from_MARS(request, settings, nc_dir, nc_file, qos):
    """
    Retrieve file from MARS
    """

    def execute(task):
        sp.call(task, shell=True, executable='/bin/bash')

    clean_name = nc_file[:-3]
    mars_req   = '{}.mars' .format(clean_name)
    grib_file  = '{}.grib' .format(clean_name)
    slurm_job  = '{}.slurm'.format(clean_name)

    # Wall clock limit
    wc_lim = '03:00:00' if qos == 'express' else '06:00:00'

    # Create MARS request
    f = open(mars_req, 'w')
    f.write('retrieve,\n')
    for key, value in request.items():
        f.write('{}={},\n'.format(key,value))
    f.write('target=\"{}\"\n'.format(grib_file))
    f.close()

    # Create SLURM job file
    date  = settings['date']
    ftype = settings['ftype'].split('_')
    jobname = '{0:04d}{1:02d}{2:02d}{3:}{4:}'.format(date.year, date.month, date.day, ftype[1], ftype[0])

    f = open(slurm_job, 'w')
    f.write('#!/bin/ksh\n')
    f.write('#SBATCH --qos={}\n'.format(qos))
    f.write('#SBATCH --job-name={}\n'.format(jobname))
    f.write('#SBATCH --output={}.%N.%j.out\n'.format(slurm_job))
    f.write('#SBATCH --error={}.%N.%j.err\n'.format(slurm_job))
    f.write('#SBATCH --workdir={}\n'.format(nc_dir))
    f.write('#SBATCH --time={}\n\n'.format(wc_lim))

    f.write('mars {}\n'.format(mars_req))
    f.write('grib_to_netcdf -o {} {}'.format(nc_file, grib_file))
    f.close()

    # Submit job
    execute('sbatch {}'.format(slurm_job))


def ERA5_file_path(year, month, day, path, case, ftype, return_dir=True):
    """
    Return saving path of files in format `path/yyyy/mm/dd/type.nc`
    """

    ERA_dir = "{0}/{1}/{2:04d}/{3:02d}/{4:02d}".format(path, case, year, month, day)
    ERA_file = "{0}/{1}.nc".format(ERA_dir, ftype)

    if return_dir:
        return ERA_dir, ERA_file
    else:
        return ERA_file


def download_ERA5_file(settings):
    """
    Download ERA5 analysis or forecasts on surface, model or pressure levels
    Requested parameters are hardcoded and chosen for the specific use of LS2D

    Arguments:
        settings : dictionary
            Dictionary with keys:
                date : datetime object with date to download
                lat, lon : requested latitude and longitude
                size : download an area of lat+/-size, lon+/-size (degrees)
                path : absolute or relative path to save the NetCDF data
                case : case name used in file name of NetCDF files
                ftype : level/forecast/analysis switch (in: [model_an, model_fc, pressure_an, surface_an])
    """

    message('Downloading: {} - {}'.format(settings['date'], settings['ftype']))

    # Output file name
    nc_dir, nc_file = ERA5_file_path(settings['date'].year, settings['date'].month, settings['date'].day,\
                                     settings['base_path'], settings['case_name'], settings['ftype'])

    # Write CDS API prints to log file (NetCDF file path/name appended with .log)
    if settings['write_log']:
        log_file   = '{}.log'.format(nc_file)
        old_stdout = sys.stdout
        sys.stdout = open(log_file, 'w')

    # Monitor the required download time
    start = datetime.datetime.now()


    # Shared set of CDS Python API settings for all download types:
    request = {
        'class'   : 'ea',
        'expver'  : '1',
        #'expver'  : '5',
        'stream'  : 'oper',
        'date'    : '{0:04d}-{1:02d}-{2:02d}'.format(settings['date'].year, settings['date'].month, settings['date'].day),
        'area'    : '{}/{}/{}/{}'.format(settings['central_lat']+settings['area_size'], settings['central_lon']-settings['area_size'],\
                                         settings['central_lat']-settings['area_size'], settings['central_lon']+settings['area_size']),
        'grid'    : '0.25/0.25',
        'format'  : 'netcdf',
    }

    # Model levels and time steps to retrieve
    model_levels = '1/to/137/by/1'
    press_levels = '1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/\
500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000'

    an_times = '0/to/23/by/1'
    fc_times = '06:00:00/18:00:00'
    fc_steps = '1/to/12/by/1'

    # Update request based on level/analysis/forecast:
    if settings['ftype'] == 'model_an':
        qos = 'normal'
        request.update({
            'levtype'  : 'ml',
            'type'     : 'an',
            'levelist' : model_levels,
            'time'     : an_times,
            'param'    : '75/76/129/130/131/132/133/135/152/246/247/248/203'
        })

    elif settings['ftype'] == 'model_fc':
        qos = 'normal'
        request.update({
            'levtype'  : 'ml',
            'type'     : 'fc',
            'levelist' : model_levels,
            'step'     : fc_steps,
            'time'     : fc_times,
            'param'    : '129/152/235001/235002/235003/235004/235005/235006/235007/235008'
        })

    elif settings['ftype'] == 'pressure_an':
        qos = 'normal'
        request.update({
            'levtype'  : 'pl',
            'type'     : 'an',
            'levelist' : press_levels,
            'time'     : an_times,
            'param'    : '129.128'
        })

    elif settings['ftype'] == 'surface_an':
        qos = 'normal'
        request.update({
            'levtype'  : 'sfc',
            'type'     : 'an',
            'time'     : an_times,
            'param'    : '15.128/16.128/17.128/18.128/27.128/28.128/29.128/30.128/34.128/35.128/36.128/37.128/38.128/39.128/40.128/41.128/42.128/43.128/66.128/67.128/74.128/78.128/79.128/89.228/90.228/129.128/134.128/136.128/137.128/139.128/151.128/160.128/161.128/162.128/163.128/164.128/165.128/166.128/167.128/168.128/170.128/172.128/183.128/186.128/187.128/188.128/198.128/229.128/230.128/231.128/232.128/235.128/236.128/243.128/244.128/245.128'
        })

    elif settings['ftype'] == 'surface_fc':
        qos = 'normal'
        request.update({
            'levtype'  : 'sfc',
            'type'     : 'fc',
            'step'     : fc_steps,
            'time'     : fc_times,
            'param'    : '21.228/22.228/35.235/36.235/37.235/38.235/39.235/40.235/49.235/50.235/51.235/52.235/53.235/58.235/59.235/68.235/69.235/129.228/130.228/169.128/175.128/176.128/177.128/178.128/179.128/208.128/209.128/210.128/211.128'
        })

    # Retrieve NetCDF file from CDS or MARS:
    if settings['data_source'] == 'CDS':
        import cdsapi
        server = cdsapi.Client()
        server.retrieve('reanalysis-era5-complete', request, nc_file)
    if settings['data_source'] == 'MARS':
        retrieve_from_MARS(request, settings, nc_dir, nc_file, qos)

    # Restore printing to screen
    if settings['write_log']:
        sys.stdout = old_stdout

    message('Finished: {} - {} in {}'.format(settings['date'], settings['ftype'], datetime.datetime.now()-start))


def download_ERA5(settings):
    """
    Download all required ERA5 fields for an experiment
    between `starttime` and `endtime`

    Analysis and forecasts are downloaded as 24 hour blocks:
        Analysis: 00 UTC to (including) 23 UTC
        Forecast: 06 UTC to (including) 05 UTC next day

    Arguments:
        start : datetime object
            Start date+time of experiment
        end : datetime object
            End date+time of experiment
        lat, lon : float
            Requested center latitude and longitude
        size : float
            Download an area of lat+/-size, lon+/-size degrees
        path : string
            Directory to save files
        case : string
            Case name used in file name of NetCDF files
    """

    header('Downloading ERA5 for period: {} to {}'.format(settings['start_date'], settings['end_date']))

    # Check if output directory exists, and ends with '/'
    if not os.path.isdir(settings['base_path']):
        error('Output directory \"{}\" does not exist!'.format(settings['base_path']))
    if settings['base_path'][-1] != '/':
        settings['base_path'] += '/'

    # Round date/time to full hours
    start = tt.lower_to_hour(settings['start_date'])
    end   = tt.lower_to_hour(settings['end_date']  )

    # Get list of required forecast and analysis times
    an_dates = tt.get_required_analysis(start, end)
    fc_dates = tt.get_required_forecast(start, end)

    # Base dictionary to pass to download function. In Python >3.3, multiprocessings Pool() can accept
    # multiple arguments. For now, keep it generic for older versions by passing all arguments inside a dict.
    download_settings = settings.copy()
    download_queue = []

    # Loop over all required files, check if there is a local version, if not add to download queue
    # Analysis files:
    for date in an_dates:
        for ftype in ['model_an', 'pressure_an', 'surface_an']:

            ERA_dir, ERA_file = ERA5_file_path(date.year, date.month, date.day, settings['base_path'], settings['case_name'], ftype)

            if not os.path.exists(ERA_dir):
                message('Creating output directory {}'.format(ERA_dir))
                os.makedirs(ERA_dir)

            if os.path.isfile(ERA_file):
                message('Found {} - {} local'.format(date, ftype))
            else:
                settings_tmp = download_settings.copy()
                settings_tmp.update({'date': date, 'ftype':ftype})
                download_queue.append(settings_tmp)

    # Forecast files
    for date in fc_dates:
        for ftype in ['model_fc', 'surface_fc']:

            ERA_dir, ERA_file = ERA5_file_path(date.year, date.month, date.day, settings['base_path'], settings['case_name'], ftype)

            if not os.path.exists(ERA_dir):
                message('Creating output directory {}'.format(ERA_dir))
                os.makedirs(ERA_dir)

            if os.path.isfile(ERA_file):
                message('Found {} - {} local'.format(date, ftype))
            else:
                settings_tmp = download_settings.copy()
                settings_tmp.update({'date': date, 'ftype':ftype})
                download_queue.append(settings_tmp)

    if settings['ntasks'] > 1:
        # Create download Pool:
        pool = mp.Pool(processes=settings['ntasks'])
        pool.map(download_ERA5_file, download_queue)
    else:
        # Simple serial retrieve:
        for req in download_queue:
            download_ERA5_file(req)



if __name__ == "__main__":
    """ Test / example, only executed if script is called directly """
    #case = 'cabauw'
    case = 'atto'
    #case = 'memo'
    #case = 'barbados'
    #case = 'k34'

    if case == 'atto':
        # ATTO
        settings = {
            'central_lat' : -2.1457,
            'central_lon' : -59.0048,
            'area_size'   : 2,
            'case_name'   : 'atto',
            'base_path'   : '/scratch/ms/nl/nkbs/LS2D/',            # ECMWF
            'start_date'  : datetime.datetime(year=2019, month=9,  day=1, hour=0),
            'end_date'    : datetime.datetime(year=2019, month=10, day=1, hour=0),
            'write_log'   : False,
            'data_source' : 'MARS',
            'ntasks'      : 1
            }

    if case == 'k34':
        # K34 tower
        settings = {
            'central_lat' : -2.609097222,
            'central_lon' : -60.20929722,
            'area_size'   : 2,
            'case_name'   : 'k34',
            'base_path'   : '/scratch/ms/nl/nkbs/LS2D/',            # ECMWF
            'start_date'  : datetime.datetime(year=2014, month=8, day=20, hour=0),
            'end_date'    : datetime.datetime(year=2014, month=9, day=21, hour=0),
            'write_log'   : False,
            'data_source' : 'MARS',
            'ntasks'      : 1
            }

    elif case == 'cabauw':
        # Cabauw
        settings = {
            'central_lat' : 51.971,
            'central_lon' : 4.927,
            'area_size'   : 1.5,
            'case_name'   : 'cabauw',
            'base_path'   : '/scratch/ms/nl/nkbs/LS2D/',            # ECMWF
            'start_date'  : datetime.datetime(year=2017, month=8, day=19, hour=0),
            'end_date'    : datetime.datetime(year=2017, month=8, day=20, hour=0),
            'write_log'   : False,
            'data_source' : 'MARS',
            'ntasks'      : 1
            }

    elif case == 'memo':
        # NEMO field experiment
        settings = {
            'central_lat' : 44.9579,
            'central_lon' : 25.7794,
            'area_size'   : 1,
            'case_name'   : 'romania',
            'base_path'   : '/scratch/ms/nl/nkbs/LS2D/',            # ECMWF
            'start_date'  : datetime.datetime(year=2019, month=10, day=17, hour=0),
            'end_date'    : datetime.datetime(year=2019, month=10, day=19, hour=0),
            'write_log'   : False,
            'data_source' : 'MARS',
            'ntasks'      : 1
            }

    elif case == 'barbados':
        # Barbados
        settings = {
            'central_lat' : 13.3,
            'central_lon' : -57.6,
            'area_size'   : 4,
            'case_name'   : 'barbados',
            'base_path'   : '/scratch/ms/nl/nkbs/LS2D/',            # ECMWF
            'start_date'  : datetime.datetime(year=2020, month=1, day=19, hour=0),
            'end_date'    : datetime.datetime(year=2020, month=2, day=21, hour=0),
            'write_log'   : False,
            'data_source' : 'MARS',
            'ntasks'      : 1
            }


    # Download the ERA5 data (or check whether it is available local)
    download_ERA5(settings)
