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
import shutil
import sys
import datetime
import os

# Third party modules
import netCDF4 as nc4
import xarray as xr
import numpy as np

# LS2D modules
sys.path.append('/home/bart/meteo/models/LS2D')
from ls2d.src.messages import *


def patch_netcdf(nc_file_path):
    """
    With the introduction of the new Copernicus Data Store (CDS) in September 2024,
    the NetCDF format for ERA5 data has undergone some changes. As a result, these
    NetCDF files differ from the previous CDS format and from files currently retrieved from MARS.
    This function patches the new NetCDF files, to make them +/- identical to the old format.

    NOTE: the patched files are not 100% identical to the old format, just
    identical enough for (LS)2D to read and parse them. 
    """

    # Backup old file, and remove original.
    backup_file_path = f'{nc_file_path}.unpatched'
    shutil.copyfile(nc_file_path, backup_file_path)
    os.remove(nc_file_path)

    # Edit with Xarray. Read the copied file, so that we can overwrite the original one.
    ds = xr.open_dataset(backup_file_path, decode_times=False)

    # Check if we actually have a new NetCDF file.
    if 'valid_time' not in ds.dims:
        error('Provided NetCDF is not a new (>09/2024) CDS file!')

    # Drop `expver`; we need to save this file in classic NetCDF4 format, which
    # does not support variable length strings.
    if 'expver' in ds.variables:
        ds = ds.drop('expver')

    file_name = os.path.basename(nc_file_path)
    
    if file_name == 'model_an.nc' or file_name == 'eac4_ml.nc':
        new_ds = ds.rename({
                'model_level': 'level',
                'valid_time': 'time'})

    elif file_name == 'pressure_an.nc':
        new_ds = ds.rename({
                'pressure_level': 'level',
                'valid_time': 'time'})

        # Yeah, somehow they thought it was a good idea to reverse the pressure levels......
        new_ds = new_ds.reindex(level=new_ds.level[::-1])

    elif file_name == 'surface_an.nc' or file_name == 'eac4_sfc.nc':
        new_ds = ds.rename({
                'valid_time': 'time'})

    else:
        error('Not a valid file type!')

    # Fix time. Old format was `hours since 1900-01-01 00:00:00.0`, new format `seconds since 1970-01-01`.
    old_ref = datetime.datetime(year=1900, month=1, day=1)
    new_ref = datetime.datetime(year=1970, month=1, day=1)

    new_ds['time'] = [(new_ref + datetime.timedelta(seconds=int(s)) - old_ref).total_seconds() / 3600. for s in new_ds.time.values]
    new_ds['time'].attrs['units'] = 'hours since 1900-01-01 00:00:00.0'

    # Remove Grib attributes.
    for v in new_ds.variables:
        da = new_ds[v]

        to_rm = []
        for attr in da.attrs:
            if 'GRIB' in attr:
                to_rm.append(attr)

        for attr in to_rm:
            del da.attrs[attr]

    # Overwrite old file.
    new_ds.to_netcdf(nc_file_path, format='NETCDF4_CLASSIC')

    return new_ds   # Just for debugging...



if __name__ == '__main__':
    """
    For debugging...
    """
    #model_an = patch_netcdf('/home/scratch1/bart/LS2D_ERA5/cabauw_test/2016/08/15/model_an.nc')
    #surf_an = patch_netcdf('/home/scratch1/bart/LS2D_ERA5/cabauw_test/2016/08/15/surface_an.nc')
    #pres_an = patch_netcdf('/home/scratch1/bart/LS2D_ERA5/cabauw_test/2016/08/15/pressure_an.nc')

    model_an = patch_netcdf('/home/scratch1/bart/LS2D_CAMS/cabauw_test/2016/08/15/eac4_ml.nc')
    surf_an = patch_netcdf('/home/scratch1/bart/LS2D_CAMS/cabauw_test/2016/08/15/eac4_sfc.nc')