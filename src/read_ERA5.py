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
import numpy as np
import netCDF4 as nc4
from scipy import interpolate
import os
import sys
import datetime

# Custom tools (in src subdirectory)
import spatial_tools as st
import time_tools as tt
from IFS_tools import IFS_tools
from conventions import ERA5_file_path
from messages import *


def format_h_since(hours):
    return datetime.timedelta(hours=float(hours)) + datetime.datetime(1900, 1, 1)


class Read_ERA:
    """
    Read the ERA5 model/pressure/surface level data,
    and optionally calculate the LES/SCM forcings
    """

    #def __init__(self, lat, lon, year, month, day, ERA_data, case_name, quiet=False):
    def __init__(self, start, end, lat, lon, path, case_name, quiet=False):
        """
        Read all the required fields to memory
        """

        header('Reading ERA5 from {} to {}'.format(start, end))

        # For now (?), only start and end at full hours
        start = tt.lower_to_hour(start)
        end   = tt.lower_to_hour(end)

        # Get list of required forecast and analysis times
        an_dates = tt.get_required_analysis(start, end)
        fc_dates = tt.get_required_forecast(start, end)

        # Check if output directory ends with '/'
        if path[-1] != '/':
            path += '/'

        # Create lists with required files
        an_sfc_files   = [ERA5_file_path(d.year, d.month, d.day, path, case_name, 'surface_an' ) for d in an_dates]
        an_model_files = [ERA5_file_path(d.year, d.month, d.day, path, case_name, 'model_an'   ) for d in an_dates]
        an_pres_files  = [ERA5_file_path(d.year, d.month, d.day, path, case_name, 'pressure_an') for d in an_dates]
        fc_model_files = [ERA5_file_path(d.year, d.month, d.day, path, case_name, 'model_fc'   ) for d in fc_dates]

        # Open NetCDF files: MFDataset automatically merges the files / time dimensions
        fsa = nc4.MFDataset(an_sfc_files  )
        fma = nc4.MFDataset(an_model_files)
        fpa = nc4.MFDataset(an_pres_files )
        fmf = nc4.MFDataset(fc_model_files)

        # Full time records in analysis and forecast files
        an_time_tmp = fsa.variables['time'][:]
        fc_time_tmp = fmf.variables['time'][:]

        # Find start and end time indices
        # ERA5 time is in hours since 1900-01-01; convert `start` and `end` to same units
        start_h_since = (start - datetime.datetime(1900, 1, 1)).total_seconds()/3600.
        end_h_since   = (end   - datetime.datetime(1900, 1, 1)).total_seconds()/3600.

        t0_an = np.abs(an_time_tmp - start_h_since).argmin()
        t1_an = np.abs(an_time_tmp - end_h_since  ).argmin()

        t0_fc = np.abs(fc_time_tmp - start_h_since).argmin()
        t1_fc = np.abs(fc_time_tmp - end_h_since  ).argmin()

        # Time slices
        t_an = np.s_[t0_an:t1_an+1]
        t_fc = np.s_[t0_fc:t1_fc+1]

        # Read spatial and time variables
        self.lats    = fma.variables['latitude'][:]
        self.lons    = fma.variables['longitude'][:]
        self.time    = fma.variables['time'][t_an]
        self.time_fc = fmf.variables['time'][t_fc]

        # Check if times are really synced, if not; quit, as things will go very wrong
        assert np.all(self.time == self.time_fc), 'Analysis and forecast times are not synced'

        # Grid and time dimensions
        self.nfull = fma.dimensions['level'].size
        self.nhalf = self.nfull+1
        self.nlat  = fma.dimensions['latitude'].size
        self.nlon  = fma.dimensions['longitude'].size
        self.ntime = self.time.size

        # Find nearest location on (regular lat/lon) grid
        self.i = np.abs(self.lons - lon).argmin()
        self.j = np.abs(self.lats - lat).argmin()

        # Some debugging output
        distance = st.haversine(self.lons[self.i], self.lats[self.j], lon, lat)
        message('Using nearest lat/lon = {0:.2f}/{1:.2f} (requested = {2:.2f}/{3:.2f}), distance = {4:.1f} km'\
                  .format(self.lats[self.j], self.lons[self.i], lat, lon, distance/1000.))

        # Read the full fields, reversing (flip) the height axis from top-to-bottom to bottom-to-top
        # ------------------------------
        # Model level analysis data:
        self.u   = np.flip(fma.variables['u']   [t_an, :, :, :], axis=1)  # u-component wind (m s-1)
        self.v   = np.flip(fma.variables['v']   [t_an, :, :, :], axis=1)  # v-component wind (m s-1)
        self.w   = np.flip(fma.variables['w']   [t_an, :, :, :], axis=1)  # Vertical velocity (Pa s-1)
        self.T   = np.flip(fma.variables['t']   [t_an, :, :, :], axis=1)  # Absolute temperature (K)
        self.q   = np.flip(fma.variables['q']   [t_an, :, :, :], axis=1)  # Specific humidity (kg kg-1)
        self.qc  = np.flip(fma.variables['clwc'][t_an, :, :, :], axis=1)  # Specific cloud liquid water content (kg kg-1)
        self.qi  = np.flip(fma.variables['ciwc'][t_an, :, :, :], axis=1)  # Specific cloud ice water content (kg kg-1)
        self.qr  = np.flip(fma.variables['crwc'][t_an, :, :, :], axis=1)  # Specific cloud rain water content (kg kg-1)
        self.qs  = np.flip(fma.variables['cswc'][t_an, :, :, :], axis=1)  # Specific cloud snow water content (kg kg-1)
        lnps     = fma.variables['lnsp'][t_an, 0, :, :]                   # Logaritm of surface pressure

        # Model level forecast data:

        # Surface variables:
        self.Ts  = fsa.variables['skt'] [t_an, :, :]          # Skin temperature (K)
        self.H   =-fsa.variables['ishf'][t_an, :, :]          # Surface sensible heat flux (W m-2)
        self.wqs =-fsa.variables['ie']  [t_an, :, :]          # Surface kinematic moisture flux (g kg-1)
        self.cc  = fsa.variables['tcc'] [t_an, :, :]          # Total cloud cover (-)
        self.z0m = fsa.variables['fsr'] [t_an, :, :]          # Surface roughness length (m)
        self.z0h = np.exp(fsa.variables['flsr'][t_an, :, :])  # Surface roughness length heat (m)

        # Pressure level data:
        self.z_p = np.flip(fpa.variables['z'][t_an, :, :, :], axis=1) / IFS_tools.grav  # Geopotential height on pressure levels
        self.p_p = np.flip(fpa.variables['level'][:],         axis=0) * 100.            # Pressure levels (Pa)

        # Calculate derived variables:
        # ------------------------------
        self.ql  = self.qc + self.qi + self.qr + self.qs  # Total liquid/solid specific humidity (kg kg-1)
        self.qt  = self.q + self.ql                       # Total specific humidity (kg kg-1)
        self.ps  = np.exp(lnps)                           # Non-logaritmic surface pressure... (Pa)
        self.Tv  = IFS_tools.calc_virtual_temp(self.T, self.q, self.qc, self.qi, self.qr, self.qs) # Virtual temp on full levels (K)

        # Calculate half level pressure and heights
        self.ph  = np.zeros((self.ntime, self.nhalf, self.nlat, self.nlon))   # Half level pressure (Pa)
        self.zh  = np.zeros((self.ntime, self.nhalf, self.nlat, self.nlon))   # Half level geopotential height (m)

        # TO-DO: remove loops
        for t in range(self.ntime):
            for la in range(self.nlat):
                for lo in range(self.nlon):
                    self.ph[t,:,la,lo] = IFS_tools.calc_half_level_pressure(self.ps[t,la,lo])
                    self.zh[t,:,la,lo] = IFS_tools.calc_half_level_Zg(self.ph[t,:,la,lo], self.Tv[t,:,la,lo])

        # Full level pressure and height as interpolation of the half level values
        self.p  = 0.5 * (self.ph[:,1:,:,:] + self.ph[:,:-1:,:])  # Full level pressure (Pa)
        self.z  = 0.5 * (self.zh[:,1:,:,:] + self.zh[:,:-1:,:])  # Full level height (m)

        # Some more derived quantities
        self.exn  = IFS_tools.calc_exner(self.p)                                   # Exner (-)
        self.th   = (self.T / self.exn)                                            # Potential temperature (K)
        self.thl  = self.th - IFS_tools.Lv / (IFS_tools.cpd * self.exn) * self.ql  # Liquid water potential temperature (K)
        self.rho  = self.p / (IFS_tools.Rd * self.Tv)                              # Density at full levels (kg m-3)
        self.wls  = -self.w / (self.rho * IFS_tools.grav)                          # Vertical velocity (m s-1)
        self.U    = (self.u**2. + self.v**2)**0.5                                  # Absolute horizontal wind (m s-1)

        self.Tvs  = IFS_tools.calc_virtual_temp(self.Ts, self.q[:,0])              # Estimate surface Tv using lowest model q (...)
        self.rhos = self.ph[:,0] / (IFS_tools.Rd * self.Tvs)                       # Surface density (kg m-3)
        self.exns = IFS_tools.calc_exner(self.ps)                                  # Exner at surface (-)
        self.wths = self.H / (self.rhos * IFS_tools.cpd * self.exns)               # Surface kinematic heat flux (K m s-1)

        self.fc  = 2 * 7.2921e-5 * np.sin(np.deg2rad(lat))                         # Coriolis parameter

    def calculate_forcings(self, n_av=1):
        """
        Calculate the advective tendencies, geostrophic wind, ....
        """
        header('Calculating large-scale forcings')

        # Slicing tuples to calculate means
        center4d = np.s_[:, :, self.j-n_av:self.j+n_av+1, self.i-n_av:self.i+n_av+1]
        center3d = np.s_[:,    self.j-n_av:self.j+n_av+1, self.i-n_av:self.i+n_av+1]

        # Slicing of boxes east, west, north and south of main domain
        box_size = 2*n_av+1
        east  = np.s_[:, :, self.j-n_av:self.j+n_av+1,  self.i+1:self.i+box_size+1]
        west  = np.s_[:, :, self.j-n_av:self.j+n_av+1,  self.i-box_size:self.i    ]
        north = np.s_[:, :, self.j-box_size:self.j,     self.i-n_av:self.i+n_av+1 ]
        south = np.s_[:, :, self.j+1:self.j+box_size+1, self.i-n_av:self.i+n_av+1 ]

        # 1. Main domain (location +/- n_av grid points)
        self.z_mean   = self.z   [center4d].mean(axis=(2,3))
        self.p_mean   = self.p   [center4d].mean(axis=(2,3))
        self.thl_mean = self.thl [center4d].mean(axis=(2,3))
        self.qt_mean  = self.qt  [center4d].mean(axis=(2,3))
        self.u_mean   = self.u   [center4d].mean(axis=(2,3))
        self.v_mean   = self.v   [center4d].mean(axis=(2,3))
        self.U_mean   = self.U   [center4d].mean(axis=(2,3))
        self.wls_mean = self.wls [center4d].mean(axis=(2,3))
        self.rho_mean = self.rho [center4d].mean(axis=(2,3))

        self.ps_mean  = self.ps  [center3d].mean(axis=(1,2))
        self.wth_mean = self.wths[center3d].mean(axis=(1,2))
        self.wq_mean  = self.wqs [center3d].mean(axis=(1,2))
        self.ps_mean  = self.ps  [center3d].mean(axis=(1,2))
        self.cc_mean  = self.cc  [center3d].mean(axis=(1,2))

        self.z0m_mean = self.z0m [center3d].mean(axis=(1,2))
        self.z0h_mean = self.z0h [center3d].mean(axis=(1,2))

        # 2. Calculate advective tendencies
        # Distance east-west and north_south of boxes
        distance_WE = st.dlon(self.lons[self.i-n_av-1], self.lons[self.i+n_av+1], self.lats[self.j])
        distance_NS = st.dlat(self.lats[self.j+n_av+1], self.lats[self.j-n_av-1])

        # Liquid water potential temperature
        self.thl_advec_x = -self.u_mean * (self.thl[east] .mean(axis=(2,3)) - self.thl[west ].mean(axis=(2,3))) / distance_WE
        self.thl_advec_y = -self.v_mean * (self.thl[north].mean(axis=(2,3)) - self.thl[south].mean(axis=(2,3))) / distance_NS
        self.thl_advec   = self.thl_advec_x + self.thl_advec_y

        # Total specific humidity
        self.qt_advec_x = -self.u_mean * (self.qt[east] .mean(axis=(2,3)) - self.qt[west ].mean(axis=(2,3))) / distance_WE
        self.qt_advec_y = -self.v_mean * (self.qt[north].mean(axis=(2,3)) - self.qt[south].mean(axis=(2,3))) / distance_NS
        self.qt_advec   = self.qt_advec_x + self.qt_advec_y

        # Momentum
        self.u_advec_x = -self.u_mean * (self.u[east] .mean(axis=(2,3)) - self.u[west ].mean(axis=(2,3))) / distance_WE
        self.u_advec_y = -self.v_mean * (self.u[north].mean(axis=(2,3)) - self.u[south].mean(axis=(2,3))) / distance_NS
        self.u_advec   = self.u_advec_x + self.u_advec_y

        self.v_advec_x = -self.u_mean * (self.v[east] .mean(axis=(2,3)) - self.v[west ].mean(axis=(2,3))) / distance_WE
        self.v_advec_y = -self.v_mean * (self.v[north].mean(axis=(2,3)) - self.v[south].mean(axis=(2,3))) / distance_NS
        self.v_advec   = self.v_advec_x + self.v_advec_y

        # 3. Geostrophic wind (gradient geopotential height on constant pressure levels)
        ug_p = -IFS_tools.grav / self.fc * (self.z_p[north].mean(axis=(2,3)) - self.z_p[south].mean(axis=(2,3))) / distance_NS
        vg_p =  IFS_tools.grav / self.fc * (self.z_p[east ].mean(axis=(2,3)) - self.z_p[west ].mean(axis=(2,3))) / distance_WE

        # Interpolate geostrophic wind onto model grid. Use Scipy's interpolation, as it can extrapolate (in case ps > 1000 hPa)
        self.ug = np.zeros_like(self.p_mean)
        self.vg = np.zeros_like(self.p_mean)
        for t in range(self.ntime):
            self.ug[t,:] = interpolate.interp1d(self.p_p, ug_p[t,:], fill_value='extrapolate')(self.p_mean[t,:])
            self.vg[t,:] = interpolate.interp1d(self.p_p, vg_p[t,:], fill_value='extrapolate')(self.p_mean[t,:])


if __name__ == '__main__':
    """ Test / example, only executed if script is called directly """

    import matplotlib.pyplot as pl
    pl.ion()
    pl.close('all')

    lat   = 51.971
    lon   = 4.927
    #path  = '/Users/bart/meteo/data/ERA5/LS2D/'
    path  = '/nobackup/users/stratum/ERA5/LS2D/'

    start = datetime.datetime(year=2016, month=5, day=1, hour=5, minute=5)
    end   = datetime.datetime(year=2016, month=5, day=2, hour=23, minute=45)

    e5 = Read_ERA(start, end, lat, lon, path, 'cabauw')
    e5.calculate_forcings(n_av=1)
