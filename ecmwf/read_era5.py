#
# This file is part of LS2D.
#
# Copyright (c) 2017-2021 Wageningen University & Research
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
import numpy as np
import netCDF4 as nc4
from scipy import interpolate

# LS2D modules
import src.spatial_tools as spatial
import src.finite_difference as fd
from src.messages import *

import ecmwf.era_tools as era_tools
from ecmwf.IFS_tools import IFS_tools

class Slice:
    def __init__(self, istart, iend, jstart, jend):
        self.istart = istart
        self.iend   = iend
        self.jstart = jstart
        self.jend   = jend

    def __call__(self, dj, di):
        return np.s_[:,:,self.jstart+dj:self.jend+dj,\
                         self.istart+di:self.iend+di]

class Read_era5:
    """
    Read the ERA5 model/pressure/surface level data,
    and optionally calculate the LES/SCM forcings
    """

    def __init__(self, settings):

        self.settings = settings
        self.start = settings['start_date']
        self.end   = settings['end_date']

        # Open all required NetCDF files:
        self.open_netcdf_files()

        # Read all the required variables:
        self.read_data()

        # Calculate derived properties needed for LES:
        self.calc_derived_data()


    def open_netcdf_files(self):
        """
        Open all NetCDF files required for start->end period
        """
        header('Reading ERA5 from {} to {}'.format(self.start, self.end))

        # Get list of required forecast and analysis times
        an_dates = era_tools.get_required_analysis(self.start, self.end)
        fc_dates = era_tools.get_required_forecast(self.start, self.end)

        # Check if output directory ends with '/'
        if self.settings['era5_path'][-1] != '/':
            self.settings['era5_path'] += '/'

        # Create lists with required files
        path = self.settings['era5_path']
        case = self.settings['case_name']
        an_sfc_files   = [era_tools.era5_file_path(
            d.year, d.month, d.day, path, case, 'surface_an',  False) for d in an_dates]
        an_model_files = [era_tools.era5_file_path(
            d.year, d.month, d.day, path, case, 'model_an',    False) for d in an_dates]
        an_pres_files  = [era_tools.era5_file_path(
            d.year, d.month, d.day, path, case, 'pressure_an', False) for d in an_dates]
        fc_model_files = [era_tools.era5_file_path(
            d.year, d.month, d.day, path, case, 'model_fc',    False) for d in fc_dates]

        # Check if all files exist, and exit if not..
        def check_files(files):
            file_missing = False
            for f in files:
                if not os.path.exists(f):
                    error('File \"{}\" does not exist...'.format(f), exit=False)
                    file_missing = True
            return file_missing

        files_missing = False
        files_missing += check_files(an_sfc_files  )
        files_missing += check_files(an_model_files)
        files_missing += check_files(an_pres_files )
        files_missing += check_files(fc_model_files)
        if files_missing:
            error('One or more required ERA5 files are missing..')

        # Open NetCDF files: MFDataset automatically merges the files / time dimensions
        self.fsa = nc4.MFDataset(an_sfc_files,   aggdim='time')
        self.fma = nc4.MFDataset(an_model_files, aggdim='time')
        self.fpa = nc4.MFDataset(an_pres_files,  aggdim='time')
        self.fmf = nc4.MFDataset(fc_model_files, aggdim='time')


    def read_data(self):
        """
        Read all the required variables from the NetCDF files
        """

        def flip(array):
            """
            Flip the height and/or latitude dimensions
            """
            if len(array.shape) == 4:
                # Reverse order of 4-dimensional field (time, height, lat, lon)
                # in height (axis=1) and lat (axis=2) direction
                return np.flip(np.flip(array, axis=1), axis=2)
            elif len(array.shape) == 3:
                # Reverse order of 3-dimensional field (time, lat, lon)
                # in lat (axis=1) direction
                return np.flip(array, axis=1)
            elif len(array.shape) == 1:
                # Reverse order of 1-dimensional field (height)
                return np.flip(array, axis=0)


        def get_variable(nc, var, s, wrap_func=None, dtype=None):
            """
            Read NetCDF variable, and flip height and latitude dimensions.
            Optionally, apply the `wrap_func` function on the data,
            and/or cast data to requested `dtype`.
            """
            data = flip(nc.variables[var][s])
            # Apply wrapper function (if provided):
            data = wrap_func(data) if wrap_func is not None else data
            # Cast to requested data type (if provided):
            data = data.astype(dtype) if dtype is not None else data

            return data


        # Full time records in analysis and forecast files
        an_time_tmp = self.fsa.variables['time'][:]
        fc_time_tmp = self.fmf.variables['time'][:]

        # Find start and end time indices
        # ERA5 time is in hours since 1900-01-01; convert `start` and `end` to same units
        date_00 = datetime.datetime(year=1900, month=1, day=1, hour=0)
        start_h_since = (self.start - date_00).total_seconds()/3600.
        end_h_since   = (self.end   - date_00).total_seconds()/3600.

        t0_an = np.abs(an_time_tmp - start_h_since).argmin()
        t1_an = np.abs(an_time_tmp - end_h_since  ).argmin()

        t0_fc = np.abs(fc_time_tmp - start_h_since).argmin()
        t1_fc = np.abs(fc_time_tmp - end_h_since  ).argmin()

        # Time slices
        t_an = np.s_[t0_an:t1_an+1]
        t_fc = np.s_[t0_fc:t1_fc+1]

        # Read spatial and time variables
        self.lats     = self.fma.variables['latitude'][::-1]
        self.lons     = self.fma.variables['longitude'][:]
        self.time     = self.fma.variables['time'][t_an]
        self.time_fc  = self.fmf.variables['time'][t_fc]

        # Shift grid from 0-360 to -180, 180 (if needed)
        if np.any(self.lons>180):
            self.lons = -360+self.lons

        self.time_sec = (self.time-self.time[0])*3600.

        # Time in datetime format
        self.datetime = [datetime.datetime(1900, 1, 1) + datetime.timedelta(hours=int(h)) for h in self.time]

        # Check if times are really synced, if not; quit, as things will go very wrong
        if self.time.size != self.time_fc.size:
            error('ERA5 analysis and forecast times are not synced (different size)')
        if not np.all(self.time == self.time_fc):
            error('ERA5 analysis and forecast times are not synced (different times)')

        # Grid and time dimensions
        self.nfull = self.fma.dimensions['level'].size
        self.nhalf = self.nfull+1
        self.nlat  = self.fma.dimensions['latitude'].size
        self.nlon  = self.fma.dimensions['longitude'].size
        self.ntime = self.time.size

        # Read the full fields, reversing (flip) the height axis from top-to-bottom
        # to bottom-to-top, and reversing the latitude dimension
        s1d  = np.s_[:         ]    # Slice for 1D fields
        s2d  = np.s_[t_an,:,:  ]    # Slice for 2D (surface) fields
        s3d  = np.s_[t_an,:,:,:]    # Slice for 3D (atmospheric) fields
        s3ds = np.s_[t_an,0,:,:]    # Slice for surface plane of 3D field

        s3d_fc = np.s_[t_fc,:,:,:]    # Slice for 3D (atmospheric) fields

        # Model level analysis data:
        self.u  = get_variable(self.fma, 'u',    s3d )  # v-component wind (m s-1)
        self.v  = get_variable(self.fma, 'v',    s3d )  # v-component wind (m s-1)
        self.w  = get_variable(self.fma, 'w',    s3d )  # Vertical velocity (Pa s-1)
        self.T  = get_variable(self.fma, 't',    s3d )  # Absolute temperature (K)
        self.q  = get_variable(self.fma, 'q',    s3d )  # Specific humidity (kg kg-1)
        self.qc = get_variable(self.fma, 'clwc', s3d )  # Specific cloud liquid water content (kg kg-1)
        self.qi = get_variable(self.fma, 'ciwc', s3d )  # Specific cloud ice content (kg kg-1)
        self.qr = get_variable(self.fma, 'crwc', s3d )  # Specific rain water content (kg kg-1)
        self.qs = get_variable(self.fma, 'cswc', s3d )  # Specific snow content (kg kg-1)
        self.o3 = get_variable(self.fma, 'o3',   s3d )  # Ozone (kg kg-1)
        self.ps = get_variable(self.fma, 'lnsp', s3ds, np.exp)  # Surface pressure (Pa)

        # Surface variables:
        self.sst =  get_variable(self.fsa, 'sst',  s2d)  # Sea surface temperature (K)
        self.Ts  =  get_variable(self.fsa, 'skt',  s2d)  # Skin temperature (K)
        self.H   = -get_variable(self.fsa, 'ishf', s2d)  # Surface sensible heat flux (W m-2)
        self.wqs = -get_variable(self.fsa, 'ie',   s2d)  # Surface kinematic moisture flux (g kg-1)
        self.z0m =  get_variable(self.fsa, 'fsr',  s2d)  # Surface roughness length (m)
        self.z0h =  get_variable(self.fsa, 'flsr', s2d, np.exp)  # Surface roughness length heat (m)

        self.soil_type     = get_variable(self.fsa, 'slt', s2d, np.round, np.int)  # Soil type (-)
        self.veg_type_low  = get_variable(self.fsa, 'tvl', s2d, np.round, np.int)  # Low vegetation type (-)
        self.veg_type_high = get_variable(self.fsa, 'tvh', s2d, np.round, np.int)  # High vegetation type (-)

        self.lai_low  = get_variable(self.fsa, 'lai_lv', s2d)  # LAI low veg (-)
        self.lai_high = get_variable(self.fsa, 'lai_hv', s2d)  # LAI high veg (-)

        # Soil variables:
        self.T_soil1 = get_variable(self.fsa, 'stl1', s2d)  # Top soil layer temperature (K)
        self.T_soil2 = get_variable(self.fsa, 'stl2', s2d)  # 2nd soil layer temperature (K)
        self.T_soil3 = get_variable(self.fsa, 'stl3', s2d)  # 3rd soil layer temperature (K)
        self.T_soil4 = get_variable(self.fsa, 'stl4', s2d)  # Bottom soil layer temperature (K)

        self.theta_soil1 = get_variable(self.fsa, 'swvl1', s2d)  # Top soil layer moisture (-)
        self.theta_soil2 = get_variable(self.fsa, 'swvl2', s2d)  # 2nd soil layer moisture (-)
        self.theta_soil3 = get_variable(self.fsa, 'swvl3', s2d)  # 3rd soil layer moisture (-)
        self.theta_soil4 = get_variable(self.fsa, 'swvl4', s2d)  # Bottom soil layer moistsure (-)

        # Pressure level data:
        self.z_p = get_variable(self.fpa, 'z', s3d) / IFS_tools.grav  # Geopotential height on pressure levels (m)
        self.p_p = get_variable(self.fpa, 'level', s1d) * 100         # Pressure levels (Pa)

        # Convert ozone from mass mixing ratio to volume mixing ratio
        self.o3 = 28.9644 / 47.9982 * self.o3
        #print('no fix o3')


    def calc_derived_data(self):
        """
        Calculate derived properties; conversion model levels to pressure/height,
        prognostic variables used by LES, etc.
        """

        self.ql  = self.qc + self.qi + self.qr + self.qs  # Total liquid/solid specific humidity (kg kg-1)
        self.qt  = self.q + self.ql                       # Total specific humidity (kg kg-1)
        self.Tv  = IFS_tools.calc_virtual_temp(
                self.T, self.q, self.qc, self.qi, self.qr, self.qs)  # Virtual temp on full levels (K)

        # Calculate half level pressure and heights
        self.ph  = np.zeros((self.ntime, self.nhalf, self.nlat, self.nlon))  # Half level pressure (Pa)
        self.zh  = np.zeros((self.ntime, self.nhalf, self.nlat, self.nlon))  # Half level geopotential height (m)

        # TO-DO: remove loops
        for t in range(self.ntime):
            for la in range(self.nlat):
                for lo in range(self.nlon):
                    self.ph[t,:,la,lo] = IFS_tools.calc_half_level_pressure(self.ps[t,la,lo])
                    self.zh[t,:,la,lo] = IFS_tools.calc_half_level_Zg(self.ph[t,:,la,lo], self.Tv[t,:,la,lo])

        # Full level pressure and height as interpolation of the half level values
        self.p  = 0.5 * (self.ph[:,1:,:,:] + self.ph[:,:-1:,:])  # Full level pressure (Pa)
        self.z  = 0.5 * (self.zh[:,1:,:,:] + self.zh[:,:-1:,:])  # Full level height (m)

        # Other derived quantities
        self.exn  = IFS_tools.calc_exner(self.p)  # Exner on full model levels (-)
        self.th   = (self.T / self.exn)  # Potential temperature (K)
        self.thl  = self.th - IFS_tools.Lv / (IFS_tools.cpd * self.exn) * self.ql  # Liquid water potential temperature (K)
        self.rho  = self.p / (IFS_tools.Rd * self.Tv)  # Density at full levels (kg m-3)
        self.wls  = -self.w / (self.rho * IFS_tools.grav)  # Vertical velocity (m s-1)
        self.U    = (self.u**2. + self.v**2)**0.5  # Absolute horizontal wind (m s-1)

        self.Tvs  = IFS_tools.calc_virtual_temp(self.Ts, self.q[:,0])  # Estimate surface Tv using lowest model q (...)
        self.rhos = self.ph[:,0] / (IFS_tools.Rd * self.Tvs)  # Surface density (kg m-3)
        self.exns = IFS_tools.calc_exner(self.ps)  # Exner at surface (-)
        self.wths = self.H / (self.rhos * IFS_tools.cpd * self.exns)  # Surface kinematic heat flux (K m s-1)

        self.fc  = 2 * 7.2921e-5 * np.sin(np.deg2rad(self.settings['central_lat']))  # Coriolis parameter

        # Store soil temperature, and moisture content, in 3D array
        self.T_soil = np.zeros((self.ntime, 4, self.nlat, self.nlon))
        self.theta_soil = np.zeros((self.ntime, 4, self.nlat, self.nlon))

        self.T_soil[:,0,:,:] = self.T_soil1[:,:,:]
        self.T_soil[:,1,:,:] = self.T_soil2[:,:,:]
        self.T_soil[:,2,:,:] = self.T_soil3[:,:,:]
        self.T_soil[:,3,:,:] = self.T_soil4[:,:,:]

        self.theta_soil[:,0,:,:] = self.theta_soil1[:,:,:]
        self.theta_soil[:,1,:,:] = self.theta_soil2[:,:,:]
        self.theta_soil[:,2,:,:] = self.theta_soil3[:,:,:]
        self.theta_soil[:,3,:,:] = self.theta_soil4[:,:,:]


    def calculate_forcings(self, n_av=0, method='4th'):
        """
        Calculate the advective tendencies, geostrophic wind, et cetera.
        """
        header('Calculating large-scale forcings')

        # Find nearest location on (regular lat/lon) grid
        self.i = np.abs(self.lons - self.settings['central_lon']).argmin()
        self.j = np.abs(self.lats - self.settings['central_lat']).argmin()

        # Some debugging output
        distance = spatial.haversine(
                self.lons[self.i], self.lats[self.j],
                self.settings['central_lon'], self.settings['central_lat'])

        message('Using nearest lat/lon = {0:.2f}/{1:.2f} (requested = {2:.2f}/{3:.2f}), distance ~= {4:.1f} km'\
                .format(self.lats[self.j], self.lons[self.i],
                        self.settings['central_lat'], self.settings['central_lon'], distance/1000.))

        # Start and end indices of averaging domain:
        istart = self.i - n_av
        iend   = self.i + n_av + 1
        jstart = self.j - n_av
        jend   = self.j + n_av + 1

        # Numpy slicing tuples for averaging domain
        center4d = np.s_[:,:,jstart:jend,istart:iend]
        center3d = np.s_[:,  jstart:jend,istart:iend]

        # Variables averaged from (time, height, lon, lat) to (time, height):
        var_4d_mean = [
                'z', 'zh', 'p', 'ph', 'T', 'thl', 'qt', 'qc', 'qi',
                'u', 'v', 'U', 'wls', 'rho', 'o3',
                'T_soil', 'theta_soil']
        for var in var_4d_mean:
            mean = getattr(self, var)[center4d].mean(axis=(2,3))
            setattr(self, '{}_mean'.format(var), mean)

        # Variables averaged from (time, lon, lat) to (time):
        var_3d_mean = [
                'ps', 'Ts', 'sst', 'wths', 'wqs', 'ps',
                'lai_low', 'lai_high', 'z0m', 'z0h']
        for var in var_3d_mean:
            mean = getattr(self, var)[center3d].mean(axis=(1,2))
            setattr(self, '{}_mean'.format(var), mean)

        # Variables selected as nearest-neighbour
        var_nn = ['soil_type', 'veg_type_low', 'veg_type_high']
        for var in var_nn:
            data = getattr(self, var)
            setattr(self, '{}_nn'.format(var), data[:, self.j, self.i])

        # Half level values temperature for radiation
        self.Th_mean = np.zeros_like(self.zh_mean)
        self.Th_mean[:,1:-1] = 0.5 * (self.T_mean[:,1:] + self.T_mean[:,:-1])

        dTdz = (self.Th_mean[:,1] - self.T_mean[:,0]) / (self.zh_mean[:,1] - self.z_mean[:,0])
        self.Th_mean[:,0] = self.T_mean[:,0] - dTdz * self.z_mean[:,0]

        dTdz = (self.T_mean[:,-1] - self.Th_mean[:,-2]) / (self.z_mean[:,-1] - self.zh_mean[:,-2])
        self.Th_mean[:,-1] = self.T_mean[:,-1] + dTdz * (self.zh_mean[:,-1] - self.z_mean[:,-1])

        # Estimate horizontal grid spacing (assumed constant in averaging domain)\
        dx = spatial.dlon(self.lons[self.i-1], self.lons[self.i+1], self.lats[self.j]) / 2.
        dy = spatial.dlat(self.lats[self.j-1], self.lats[self.j+1]) / 2.

        if (method == '2nd'):

            r_earth = 6.37e6

            lat_rad = np.deg2rad(self.lats)
            lon_rad = np.deg2rad(self.lons)
            cos_lat = np.cos(lat_rad)

            dxdi = np.zeros((self.nlat, self.nlon))
            dydj = np.zeros((self.nlat, self.nlon))

            dxdi[:,:] = r_earth * cos_lat[:,None]*np.gradient(lon_rad[None, :], axis=1)
            dydj[:,:] = r_earth * np.gradient(lat_rad[:, None], axis=0)

            def advec(var):
                dvardx = np.gradient(var, axis=3) / dxdi[None, None, :, :]
                dvardy = np.gradient(var, axis=2) / dydj[None, None, :, :]
                dtvar  = -self.u * dvardx - self.v * dvardy
                return dtvar[center4d].mean(axis=(2,3))

            # Calculate advective tendencies:
            self.dtthl_advec_mean = advec(self.thl)
            self.dtqt_advec_mean  = advec(self.qt)
            self.dtu_advec_mean   = advec(self.u)
            self.dtv_advec_mean   = advec(self.v)

            # Geostrophic wind:
            dzdx = np.gradient(self.z_p, axis=3) / dxdi[None, None, :, :]
            dzdy = np.gradient(self.z_p, axis=2) / dydj[None, None, :, :]

            ug = -IFS_tools.grav / self.fc * dzdy
            vg =  IFS_tools.grav / self.fc * dzdx

            ug_p_mean = ug[center4d].mean(axis=(2,3))
            vg_p_mean = vg[center4d].mean(axis=(2,3))

        #if (method == '2nd'):

        #    s = Slice(istart, iend, jstart, jend)

        #    # Calculate advective tendencies
        #    self.dtthl_advec_mean = (
        #        -self.u[s(0,0)] * fd.grad2c( self.thl[s(0,-1)], self.thl[s(0,+1)], dx) \
        #        -self.v[s(0,0)] * fd.grad2c( self.thl[s(-1,0)], self.thl[s(+1,0)], dy) ).mean(axis=(2,3))

        #    self.dtqt_advec_mean = (
        #        -self.u[s(0,0)] * fd.grad2c( self.qt[s(0,-1)], self.qt[s(0,+1)], dx) \
        #        -self.v[s(0,0)] * fd.grad2c( self.qt[s(-1,0)], self.qt[s(+1,0)], dy) ).mean(axis=(2,3))

        #    self.dtu_advec_mean = (
        #        -self.u[s(0,0)] * fd.grad2c( self.u[s(0,-1)], self.u[s(0,+1)], dx) \
        #        -self.v[s(0,0)] * fd.grad2c( self.u[s(-1,0)], self.u[s(+1,0)], dy) ).mean(axis=(2,3))

        #    self.dtv_advec_mean = (
        #        -self.u[s(0,0)] * fd.grad2c( self.v[s(0,-1)], self.v[s(0,+1)], dx) \
        #        -self.v[s(0,0)] * fd.grad2c( self.v[s(-1,0)], self.v[s(+1,0)], dy) ).mean(axis=(2,3))

        #    # Geostrophic wind (gradient geopotential height on constant pressure levels)
        #    vg_p_mean = (  IFS_tools.grav / self.fc * fd.grad2c(
        #        self.z_p[s(0,-1)], self.z_p[s(0,+1)], dx) ).mean(axis=(2,3))
        #    ug_p_mean = ( -IFS_tools.grav / self.fc * fd.grad2c(
        #        self.z_p[s(-1,0)], self.z_p[s(+1,0)], dy) ).mean(axis=(2,3))

        elif (method == '4th'):

            s = Slice(istart, iend, jstart, jend)

            # Calculate advective tendencies
            self.dtthl_advec_mean = (
                -self.u[s(0,0)] * fd.grad4c(
                    self.thl[s(0,-2)], self.thl[s(0,-1)], self.thl[s(0,+1)], self.thl[s(0,+2)], dx) \
                -self.v[s(0,0)] * fd.grad4c(
                    self.thl[s(-2,0)], self.thl[s(-1,0)], self.thl[s(+1,0)], self.thl[s(+2,0)], dy)
                                    ).mean(axis=(2,3))

            self.dtqt_advec_mean =  (
                -self.u[s(0,0)] * fd.grad4c(
                    self.qt[s(0,-2)], self.qt[s(0,-1)], self.qt[s(0,+1)], self.qt[s(0,+2)], dx) \
                -self.v[s(0,0)] * fd.grad4c(
                    self.qt[s(-2,0)], self.qt[s(-1,0)], self.qt[s(+1,0)], self.qt[s(+2,0)], dy)
                                    ).mean(axis=(2,3))

            self.dtu_advec_mean = (
                -self.u[s(0,0)] * fd.grad4c(
                    self.u[s(0,-2)], self.u[s(0,-1)], self.u[s(0,+1)], self.u[s(0,+2)], dx) \
                -self.v[s(0,0)] * fd.grad4c(
                    self.u[s(-2,0)], self.u[s(-1,0)], self.u[s(+1,0)], self.u[s(+2,0)], dy)
                                  ).mean(axis=(2,3))

            self.dtv_advec_mean = (
                -self.u[s(0,0)] * fd.grad4c(
                    self.v[s(0,-2)], self.v[s(0,-1)], self.v[s(0,+1)], self.v[s(0,+2)], dx) \
                -self.v[s(0,0)] * fd.grad4c(
                    self.v[s(-2,0)], self.v[s(-1,0)], self.v[s(+1,0)], self.v[s(+2,0)], dy)
                                  ).mean(axis=(2,3))

            # Geostrophic wind (gradient geopotential height on constant pressure levels)
            vg_p_mean = (
                IFS_tools.grav / self.fc * fd.grad4c(
                    self.z_p[s(0,-2)], self.z_p[s(0,-1)], self.z_p[s(0,+1)], self.z_p[s(0,+2)], dx)
                        ).mean(axis=(2,3))
            ug_p_mean = (
               -IFS_tools.grav / self.fc * fd.grad4c(
                    self.z_p[s(-2,0)], self.z_p[s(-1,0)], self.z_p[s(+1,0)], self.z_p[s(+2,0)], dy)
                        ).mean(axis=(2,3))

        # Interpolate geostrophic wind onto model grid. Use Scipy's interpolation,
        # as it can extrapolate (in case ps > 1000 hPa)
        self.ug_mean = np.zeros_like(self.p_mean)
        self.vg_mean = np.zeros_like(self.p_mean)

        for t in range(self.ntime):
            self.ug_mean[t,:] = interpolate.interp1d(
                    self.p_p, ug_p_mean[t,:], fill_value='extrapolate')(self.p_mean[t,:])
            self.vg_mean[t,:] = interpolate.interp1d(
                    self.p_p, vg_p_mean[t,:], fill_value='extrapolate')(self.p_mean[t,:])

        # Momentum tendency coriolis
        self.dtu_coriolis_mean = +self.fc * (self.v_mean - self.vg_mean)
        self.dtv_coriolis_mean = -self.fc * (self.u_mean - self.ug_mean)

        # Total momentum tendency
        self.dtu_total_mean = self.dtu_advec_mean + self.dtu_coriolis_mean
        self.dtv_total_mean = self.dtv_advec_mean + self.dtv_coriolis_mean


    def interpolate_to_fixed_height(self, z):
        """
        Interpolate list of `variables` to a requested (fixed in time) height `z`.
        """

        variables = ['thl', 'qt', 'u', 'v', 'wls', 'p',
            'dtthl_advec', 'dtqt_advec', 'dtu_advec', 'dtv_advec',
            'ug' ,'vg' ,'o3', 'z']

        def interp_z(array, z):
            out = np.empty((self.ntime, z.size))
            for t in range(self.ntime):
                out[t,:] = np.interp(z, self.z_mean[t,:], array[t,:])
            return out

        output = {}
        for var in variables:
            var_era5 = '{}_mean'.format(var)
            if hasattr(self, var_era5):
                output[var] = interp_z(getattr(self, var_era5), z)
            else:
                error('Can\'t interpolate variable \"{}\"...'.format(var))

        return output


    def save_forcings(self, file_name):
        """
        Write the large-scale forcings to a NetCDF file
        """

        def add_variable(nc, name, dims, dtype, units, long_name, data):
            var = nc.createVariable(name, dtype, dims)
            var[:] = data
            var.units = units
            var.long_name = long_name

        # Define new NetCDF file
        nc = nc4.Dataset(file_name, 'w')
        nc.createDimension('time', self.ntime)
        nc.createDimension('level', self.nfull)

        # Add variables
        add_variable(nc, 'time', ('time'), int,
                'hours since 1900-01-01 00:00:00.0', 'time', self.time)

        # Mean profiles
        add_variable(nc, 'p', ('time','level'), float, 'pa',
                'Pressure', self.p_mean)
        add_variable(nc, 'z', ('time','level'), float, 'm',
                'Height above ground level', self.z_mean)
        add_variable(nc, 'thl', ('time','level'), float, 'K',
                'Liquid water potential temperature', self.thl_mean)
        add_variable(nc, 'qt', ('time','level'), float, 'kg kg-1',
                'Total specific humidity', self.qt_mean)
        add_variable(nc, 'u', ('time','level'), float, 'm s-1',
                'Zonal wind component', self.u_mean)
        add_variable(nc, 'v', ('time','level'), float, 'm s-1',
                'Meridional wind component', self.v_mean)
        add_variable(nc, 'wls', ('time','level'), float, 'm s-1',
                'Vertical wind component', self.wls_mean)

        # Geostrophic wind
        add_variable(nc, 'ug', ('time','level'), float, 'm s-1',
                'Zonal component geostrophic wind', self.ug_mean)
        add_variable(nc, 'vg', ('time','level'), float, 'm s-1',
                'Meridional component geostrophic wind', self.vg_mean)

        # Advective tendencies
        add_variable(nc, 'dtthl_advec', ('time','level'), float, 'K s-1',
                'Advective tendency liquid water potential temperature', self.dtthl_advec_mean)
        add_variable(nc, 'dtqt_advec', ('time','level'), float, 'kg kg-1 s-1',
                'Advective tendency total specific humidity', self.dtqt_advec_mean)
        add_variable(nc, 'dtu_advec', ('time','level'), float, 'm s-2',
                'Advective tendency zonal wind component', self.dtu_advec_mean)
        add_variable(nc, 'dtv_advec', ('time','level'), float, 'm s-2',
                'Advective tendency meridional wind component', self.dtv_advec_mean)

        # Radiative tendencies
        add_variable(nc, 'dtthl_sw', ('time','level'), float, 'K s-1',
                'Shortwave radiative tendency liquid water potential temperature', self.dtthl_sw_mean)
        add_variable(nc, 'dtthl_lw', ('time','level'), float, 'K s-1',
                'Longwave radiative tendency liquid water potential temperature', self.dtthl_lw_mean)

        nc.close()
