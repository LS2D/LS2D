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
from multiprocessing import Process
import multiprocessing

try:
    from ecmwfapi import ECMWFDataServer
except:
    sys.exit('ERROR: Can\'t find the ECMWF Python api....\nSee https://software.ecmwf.int/wiki/display/WEBAPI/ECMWF+Web+API+Home')

# Custom tools (in src subdirectory)
import time_tools as tt
from messages import *

def get_download_path(year, month, day, path, case, type):
    """
    Return saving path of files in format `path/yyyymmdd_case_type.nc`
    """

    return "{0:}{1:04d}{2:02d}{3:02d}_{4:}_{5:}.nc".format(path, year, month, day, case, type)


def get_ERA5(settings):
    """
    Download ERA5 analysis or forecasts on surface, model or pressure levels
    Requested parameters are hardcoded and chosen for the specific use of LS2D

    Arguments:
        settings : dictionary
            Dictionary with keys:
                date : datetime object with date to download
                lat, lon : requested latitude and longitude
                size : download an area of lat+/-size, lon+/-size
                path : absolute or relative path to save the NetCDF data
                case : case name used in file name of NetCDF files
                ftype : level/forecast/analysis switch (in: [model_an, model_fc, pressure_an, surface_an])
    """

    message('Downloading: {} - {}'.format(settings['date'], settings['ftype']))

    # Mute the ECMWF API prints.....
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    # Time required download time
    start = datetime.datetime.now()

    # Create new instance of ECMWF Python api
    server = ECMWFDataServer()

    # Shared set of ECMWF Python api settings:
    request = {
        "class"   : "ea",
        "dataset" : "era5",
        "date"    : "{0:04d}-{1:02d}-{2:02d}".format(settings['date'].year, settings['date'].month, settings['date'].day),
        "expver"  : "1",
        "stream"  : "oper",
        "grid"    : "0.3/0.3",
        "area"    : "{}/{}/{}/{}".format(settings['lat']+settings['size'], settings['lon']-settings['size'], settings['lat']-settings['size'], settings['lon']+settings['size']),
        "format"  : "netcdf",
        "target"  : get_download_path(settings['date'].year, settings['date'].month, settings['date'].day, settings['path'], settings['case'], settings['ftype'])
    }

    # Update request based on level/analysis/forecast:
    if settings['ftype'] == 'model_an':
        request.update({
            "levtype"  : "ml",
            "levelist" : "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137",
            "param"    : "75/76/129/130/131/132/133/135/152/246/247/248",
            "time"     : "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
            "type"     : "an"
        })
    elif settings['ftype'] == "model_fc":
        request.update({
            "levtype"  : "ml",
            "levelist" : "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137",
            "param"    : "129/152/235001/235002/235003/235004",
            "step"     : "0/1/2/3/4/5/6/7/8/9/10/11",
            "time"     : "06:00:00/18:00:00",
            "type"     : "fc"
        })
    elif settings['ftype'] == 'pressure_an':
        request.update({
            "levtype"  : "pl",
            "levelist" : "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
            "param"    : "129.128",
            "time"     : "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
            "type"     : "an"
        })
    elif settings['ftype'] == 'surface_an':
        request.update({
            "levtype"  : "sfc",
            "param"    : "78.128/79.128/89.228/90.228/134.128/136.128/137.128/151.128/159.128/164.128/165.128/166.128/167.128/168.128/186.128/187.128/188.128/229.128/230.128/231.128/232.128/235.128/244.128/245.128/246.228/247.228",
            "time"     : "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
            "type"     : "an"
        })

    # Retrieve file
    server.retrieve(request)

    # Restore print
    sys.stdout = old_stdout

    message('Finished: {} - {} in {}'.format(settings['date'], settings['ftype'], datetime.datetime.now()-start))



def download_ERA5_period(start, end, lat, lon, size, path, case):
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

    header('Dowloading ERA5 for period: {} to {}'.format(start, end))

    # Check if output directory exists, and ends with '/'
    if not os.path.isdir(path):
        error('ERA5 output directory \"{}\" does not exist!'.format(path))
        sys.exit()
    if path[-1] != '/':
        path += '/'

    # Round date/time to full hours
    start = tt.lower_to_hour(start)
    end   = tt.lower_to_hour(end)

    # Get list of required forecast and analysis times
    an_dates = tt.get_required_analysis(start, end)
    fc_dates = tt.get_required_forecast(start, end)

    # Base dictionary to pass to download function. In Python >3.3, multiprocessings Pool() can accept
    # multiple arguments. For now, keep it generic for older versions by passing all arguments inside a dict
    download_settings = {'lat':lat, 'lon':lon, 'size':size, 'path':path, 'case':case}
    download_queue = []

    # Loop over all required files, check if there is a local version, if not add to download queue
    # Analysis files:
    for date in an_dates:
        for ftype in ['model_an', 'pressure_an', 'surface_an']:

            file_name = get_download_path(date.year, date.month, date.day, path, case, ftype)

            if os.path.isfile(file_name):
                message('Found {} - {} local'.format(date, ftype))
            else:
                settings = download_settings.copy()
                settings.update({'date': date, 'ftype':ftype})
                download_queue.append(settings)

    # Forecast files
    for date in fc_dates:
        for ftype in ['model_fc']:

            file_name = get_download_path(date.year, date.month, date.day, path, case, ftype)

            if os.path.isfile(file_name):
                message('Found {} - {} local'.format(date, ftype))
            else:
                settings = download_settings.copy()
                settings.update({'date': date, 'ftype':ftype})
                download_queue.append(settings)

    # Create download Pool with 3 threads (ECMWF allows up to 3 parallel requests):
    pool = multiprocessing.Pool(processes=5)
    pool.map(get_ERA5, download_queue)


if __name__ == "__main__":
    """ Test / example, only executed if script is called directly """

    lat   = 51.971
    lon   = 4.927
    size  = 2
    case  = 'cabauw'
    #path  = '/home/scratch1/meteo_data/ERA5/LS2D/'
    path  = '/nobackup/users/stratum/ERA5/LS2D/'
    #path  = '/Users/bart/meteo/data/ERA5/LS2D/'

    start = datetime.datetime(year=2016, month=5, day=4, hour=0)
    end   = datetime.datetime(year=2016, month=5, day=4, hour=23, minute=30)

    download_ERA5_period(start, end, lat, lon, size, path, case)
