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

import sys
import os
from multiprocessing import Process

try:
    from ecmwfapii import ECMWFDataServer
except:
    sys.exit('ERROR: Can\'t find the ECMWF Python api....\nSee https://software.ecmwf.int/wiki/display/WEBAPI/ECMWF+Web+API+Home')


def download_sfc(year, month, day, lat, lon, size, download_path, case):
    """
    Download the ERA5 surface data
    Data is saved in the `download_path` directory in format:
        `download_path/yyyymmdd_case_sfc.nc`

    Arguments:
        year, month, day -- date to download
        lat, lon -- location to download
        size -- download an area of lat+/-size, lon+/-size
        download_path -- directory to save the data
        case -- case name
    """

    server = ECMWFDataServer()
    server.retrieve({
        "class": "ea",
        "dataset": "era5",
        "date": "{0:04d}-{1:02d}-{2:02d}".format(year, month, day),
        "expver": "1",
        "levtype": "sfc",
        "param": "78.128/79.128/89.228/90.228/134.128/136.128/137.128/151.128/159.128/164.128/165.128/166.128/167.128/168.128/186.128/187.128/188.128/229.128/230.128/231.128/232.128/235.128/244.128/245.128/246.228/247.228",
        "stream": "oper",
        "grid": "0.3/0.3",
        "area": "{}/{}/{}/{}".format(lat+size, lon-size, lat-size, lon+size),
        "time": "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
        "type": "an",
        "format": "netcdf",
        "target": "{0:}{1:04d}{2:02d}{3:02d}_{4:}_sfc.nc".format(download_path, year, month, day, case)
    })


def download_model(year, month, day, lat, lon, size, download_path, case):
    """
    Download the ERA5 model level data
    Data is saved in the `download_path` directory in format:
        `download_path/yyyymmdd_case_model.nc`

    Arguments:
        year, month, day -- date to download
        lat, lon -- location to download
        size -- download an area of lat+/-size, lon+/-size
        download_path -- directory to save the data
        case -- case name
    """

    server = ECMWFDataServer()
    server.retrieve({
        "class": "ea",
        "dataset": "era5",
        "date": "{0:04d}-{1:02d}-{2:02d}".format(year, month, day),
        "expver": "1",
        "levelist": "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137",
        "levtype": "ml",
        "param": "75/76/129/130/131/132/133/135/152/246/247/248",
        "stream": "oper",
        "grid": "0.3/0.3",
        "area": "{}/{}/{}/{}".format(lat+size, lon-size, lat-size, lon+size),
        "time": "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
        "type": "an",
        "format": "netcdf",
        "target": "{0:}{1:04d}{2:02d}{3:02d}_{4:}_model.nc".format(download_path, year, month, day, case)
    })


def download_pressure(year, month, day, lat, lon, size, download_path, case):
    """
    Download the ERA5 pressure level data
    Data is saved in the `download_path` directory in format:
        `download_path/yyyymmdd_case_pressure.nc`

    Arguments:
        year, month, day -- date to download
        lat, lon -- location to download
        size -- download an area of lat+/-size, lon+/-size
        download_path -- directory to save the data
        case -- case name
    """

    server = ECMWFDataServer()
    server.retrieve({
        "class": "ea",
        "dataset": "era5",
        "date": "{0:04d}-{1:02d}-{2:02d}".format(year, month, day),
        "expver": "1",
        "levelist": "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
        "levtype": "pl",
        "param": "129.128",
        "stream": "oper",
        "grid": "0.3/0.3",
        "area": "{}/{}/{}/{}".format(lat+size, lon-size, lat-size, lon+size),
        "time": "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
        "type": "an",
        "format": "netcdf",
        "target": "{0:}{1:04d}{2:02d}{3:02d}_{4:}_pressure.nc".format(download_path, year, month, day, case)
    })


def download_ERA5(year, month, day, lat, lon, size, download_path, case):
    """
    Download the ERA5 surface/model/pressure level data in parallel
    Data is saved in the `download_path` directory in format:
        `download_path/yyyymmdd_case_level_type.nc`

    Arguments:
        year, month, day -- date to download
        lat, lon -- location to download
        size -- download an area of lat+/-size, lon+/-size
        download_path -- directory to save the data
        case -- case name
    """

    if download_path[-1] != '/':
        download_path += '/'

    # Check if download directory exists
    if not os.path.isdir(download_path):
        sys.exit('ERA5 download_path does not exist')

    # Push the requests in parallel, given the (sometimes) long queue times of Mars
    # ECMWF allows up to 3 requests in parallel.
    processes = []
    processes.append(Process(target=download_sfc,      args=(year, month, day, lat, lon, size, download_path, case)))
    processes.append(Process(target=download_model,    args=(year, month, day, lat, lon, size, download_path, case)))
    processes.append(Process(target=download_pressure, args=(year, month, day, lat, lon, size, download_path, case)))

    # Start all downloads:
    for p in processes:
        p.start()

    # Wait for them to finish:
    for p in processes:
        p.join()


if __name__ == "__main__":
    """ Test / example, only executed if script is called directly """

    lat   = 51.971
    lon   = 4.927
    year  = 2016
    month = 5
    day   = 1
    size  = 2
    case  = 'cabauw'
    path  = '/home/scratch1/meteo_data/ERA5/'

    download_ERA5(year, month, day, lat, lon, size, path, case) 
