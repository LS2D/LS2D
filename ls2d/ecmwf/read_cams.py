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
import datetime
import sys, os

# Third party modules
import xarray as xr
import pandas as pd
import numpy as np
from scipy import interpolate

# LS2D modules
from ls2d.src.messages import *
import ls2d.src.spatial_tools as spatial
import ls2d.ecmwf.era_tools as era_tools
from ls2d.ecmwf.IFS_tools import IFS_tools


class Read_cams:
    """
    Read the ERA5 model/pressure/surface level data,
    and optionally calculate the LES/SCM forcings
    """

    def __init__(self, settings, variables):

        self.settings = settings
        self.variables = variables
        self.start = settings['start_date']
        self.end   = settings['end_date']

        # Open all required NetCDF files:
        self.read_netcdf()


    def read_netcdf(self):
        """
        Open all NetCDF files required for start->end period using Xarray.
        """

        header('Reading CAMS from {} to {}'.format(self.start, self.end))

        # Get list of required forecast and analysis times
        an_dates = era_tools.get_required_analysis(self.start, self.end, freq=3)

        # Check if output directory ends with '/'
        if self.settings['cams_path'][-1] != '/':
            self.settings['cams_path'] += '/'

        # Create lists with required files
        path = self.settings['cams_path']
        case = self.settings['case_name']

        # Populate list of NetCDF files for each CAMS file type.
        cams_files = []
        for date in an_dates:
            for cams_type in ['eac4_ml', 'eac4_sfc', 'egg4_ml']:
                if cams_type in self.variables.keys():
                    
                    if cams_type == 'egg4_ml':
                        error('Processing CAMS EGG4 data currently does not work due to an open ADS bug.')

                    cams_files.append(era_tools.era5_file_path(
                        date.year, date.month, date.day, path, case, cams_type, False))

        # Check if all files are present.
        files_missing = False
        for f in cams_files:
            if not os.path.exists(f):
                 warning(f'File {f} does not exist...')
                 files_missing = True

        if files_missing:
            error('One or more required CAMS files are missing...!')

        # Open NetCDF files with Xarray
        ds = xr.open_mfdataset(cams_files)

        # Interpolate to hourly frequency, to stay in line with ERA5.
        # This automagically selects the correct time period as a bonus.
        dates = pd.date_range(self.start, self.end, freq='H')
        self.ds_ml = ds.interp(time=dates)

        # Reverse height dimension such that height increases with increasing levels.
        self.ds_ml = self.ds_ml.reindex(level=self.ds_ml.level[::-1])


    def calc_model_levels(self):
        """
        Calculate model level pressure and heights.
        """

        # Short-cut
        ds = self.ds_ml

        dims = self.ds_ml.dims
        dim_name = ['time', 'level', 'latitude', 'longitude']

        ntime = dims['time']
        nlevel = dims['level']
        nlon = dims['longitude']
        nlat = dims['latitude']

        # Help class for vertical grid calculations IFS.
        ifs_tools = IFS_tools('L60')

        # Calculate virtual temperature (neglecting qc et al.)
        Tv = ifs_tools.calc_virtual_temp(ds.t.values, ds.q.values)

        # Calculate half level pressure and height.
        dim_size_h = [ntime, nlevel+1, nlat, nlon]
        ph = np.zeros(dim_size_h, np.float32)
        zh = np.zeros(dim_size_h, np.float32)

        for t in range(ntime):
            for j in range(nlat):
                for i in range(nlon):
                    ph[t,:,j,i] = ifs_tools.calc_half_level_pressure(float(ds.sp[t,j,i]))
                    zh[t,:,j,i] = ifs_tools.calc_half_level_Zg(ph[t,:,j,i], Tv[t,:,j,i])

        # Full level pressure and height as interpolation of the half level values
        p = 0.5 * (ph[:,1:,:,:] + ph[:,:-1:,:])
        z = 0.5 * (zh[:,1:,:,:] + zh[:,:-1:,:])

        # Assign to Dataset.
        ds['p'] = (dim_name, p)
        ds['z'] = (dim_name, z)


    def get_les_input(self, z, n_av=0):
        """
        Interpolate variables required for LES onto model grid,
        and return xarray.Dataset with all possible LES input.
        """
    
        # Find nearest index in dataset.
        clon = self.settings['central_lon']
        clat = self.settings['central_lat']
    
        ic = int(np.abs(self.ds_ml.longitude - clon).argmin())
        jc = int(np.abs(self.ds_ml.latitude  - clat).argmin())
    
        # Some debugging output
        distance = spatial.haversine(self.ds_ml.longitude[ic], self.ds_ml.latitude[jc], clon, clat)
        message('Using nearest lat/lon = {0:.2f}/{1:.2f} (requested = {2:.2f}/{3:.2f}), distance ~= {4:.1f} km'\
                .format(self.ds_ml.latitude[jc], self.ds_ml.longitude[ic], clat, clon, distance/1000.))

        # Calculate and output averaging area.
        dlon = (1+2*n_av) * abs(float(self.ds_ml.longitude[1] - self.ds_ml.longitude[0]))
        dlat = (1+2*n_av) * abs(float(self.ds_ml.latitude[0] - self.ds_ml.latitude[1]))
        message(f'Averaging CAMS over a {dlon:.2f}°×{dlat:.2f}° spatial area.')
    
        # Slice out averaging sub-domain, and calculate mean over sub-domain.
        self.ds_ml = self.ds_ml.isel(longitude=slice(ic-n_av, ic+n_av+1), latitude=slice(jc-n_av, jc+n_av+1))

        # Load dataset to get rid of Dask/chuncks.
        self.ds_ml = self.ds_ml.load()

        # Calculate model level pressure / height.
        self.calc_model_levels()

        # Calculate Spatial mean over requested area.
        self.ds_ml_mean = self.ds_ml.mean(dim=['longitude', 'latitude'])
    
        # Create new Dataset with LES model levels as main height coordinate.    
        self.ds_les = xr.Dataset(
                coords = {
                    'time': self.ds_ml.time,
                    'z': z,
                    'lay': self.ds_ml.level.values
                })
    
        dims_sfc = ['time']
        dims_lay = ['time', 'lay']
        dims_les = ['time', 'z']
    
        ntime = self.ds_les.dims['time']
        ktot = self.ds_les.dims['z']
    
        # Parse all variables in the CAMS NetCDF files. This might be more than the variables
        # specified in settings['cams_vars'], but since those variable names differ from
        # the variable names in NetCDF, it is difficult to link them without a huge lookup
        # table containing **ALL** possible CAMS variables...
        blacklist = ['level', 'time']
    
        for name, da in self.ds_ml_mean.data_vars.items():
            if name not in blacklist and da.ndim == 2:
    
                # Add "raw" data on CAMS model levels.
                self.ds_les[f'{name}_lay'] = (dims_lay, da.values)
    
                if name != 'z':
                    # Interpolate data onto LES grid.
                    out = np.empty((ntime, ktot), np.float32)
    
                    for t in range(ntime):
                        out[t,:] = interpolate.interp1d(
                                self.ds_ml_mean['z'][t,:].values, da[t,:].values, fill_value='extrapolate')(z)
    
                    self.ds_les[name] = (dims_les, out)

        # Calculate time in seconds since start of LES.
        date = self.ds_les.time.values
        self.ds_les['time_sec'] = (date - date[0]).astype(np.float32) * 1e-9
    
        return self.ds_les
