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
from datetime import datetime
import sys, os
import subprocess
import shutil

# Third party modules
import xarray as xr
import numpy as np

# LS2D & custom modules
sys.path.append('/home/bart/meteo/models/LS2D')
import ls2d
import microhh_tools as mht

# Constants
Rd = 287.04
Rv = 461.5
ep = Rd/Rv

#
# Download ERA5 and generate LES initialisation and forcings
#
settings = {
    'central_lat' : 51.971,
    'central_lon' : 4.927,
    'area_size'   : 1,
    'case_name'   : 'cabauw',
    'era5_path'   : '/home/scratch1/meteo_data/LS2D/',
    'start_date'  : datetime(year=2016, month=8, day=15, hour=6),
    'end_date'    : datetime(year=2016, month=8, day=15, hour=18),
    'write_log'   : True,
    'data_source' : 'CDS',
    'ntasks'      : 3
    }

# Download required ERA5 files:
ls2d.download_era5(settings)

# Read ERA5 data, and calculate derived properties (thl, etc.):
era = ls2d.Read_era5(settings)

# Calculate initial conditions and large-scale forcings for LES:
era.calculate_forcings(n_av=0, method='4th')

# Define vertical grid LES:
grid = ls2d.grid.Grid_linear_stretched(kmax=176, dz0=20, alpha=0.009)
#grid.plot()

# Interpolate ERA5 variables and forcings onto LES grid:
les_input = era.interpolate_to_fixed_height(grid.z)

#
# MicroHH specific initialisation
#
# Settings:
# ------------------
# Root fraction coefficients (see IFS documentation):
a_r = 10.739
b_r = 2.608

# Soil index in `van_genuchten_parameters.nc`:
soil_index = 2

# Nudging time scale atmosphere
tau_nudge = 10800
# ------------------

# 1. RRTMGP input: ERA5 at model levels, averaged
#    in time of the entire simulation period:
z_lay = era.z_mean .mean(axis=0)
z_lev = era.zh_mean.mean(axis=0)

p_lay = era.p_mean .mean(axis=0)
p_lev = era.ph_mean.mean(axis=0)

T_lay = era.T_mean .mean(axis=0)
T_lev = era.Th_mean.mean(axis=0)

qt_rad = era.qt_mean.mean(axis=0)
o3_rad = era.o3_mean.mean(axis=0)

# Fixed concentrations:
co2 = 348.e-6
ch4 = 1650.e-9
n2o = 306.e-9
n2  = 0.7808
o2  = 0.2095

# Profiles on LES grid:
qt_les = les_input['qt'].mean(axis=0)
o3_les = les_input['o3'].mean(axis=0)

# Conversion moisture from mass to volume mixing ratio
h2o_rad = qt_rad / (ep - ep*qt_rad)
h2o_les = qt_les / (ep - ep*qt_les)

# 2. Soil
z_soil     = np.array([-1.945, -0.64, -0.175, -0.035])
soil_index = np.ones_like(z_soil)*soil_index
root_frac  = mht.calc_root_frac(z_soil, a_r, b_r)

# 3. Nudge factor
nudge_fac = np.ones(grid.z.size) / tau_nudge

#
# Write NetCDF input file for MicroHH
#
init_profiles = {
        'z': grid.z,
        'thl': les_input['thl'][0,:],
        'qt': les_input['qt'][0,:],
        'u': les_input['u'][0,:],
        'v': les_input['v'][0,:],
        'nudgefac': nudge_fac,
        'co2': co2,
        'ch4': ch4,
        'n2o': n2o,
        'n2': n2,
        'o2': o2,
        'o3': o3_les,
        'h2o': h2o_les}

radiation  = {
        'z_lay': z_lay,
        'z_lev': z_lev,
        'p_lay': p_lay,
        'p_lev': p_lev,
        't_lay': T_lay,
        't_lev': T_lev,
        'co2': co2,
        'ch4': ch4,
        'n2o': n2o,
        'n2': n2,
        'o2': o2,
        'o3': o3_rad,
        'h2o': h2o_rad}

timedep_surface = {
        'time_surface': era.time_sec,
        'p_sbot': era.ps_mean }

timedep_ls = {
        'time_ls': era.time_sec,
        'u_geo': les_input['ug'],
        'v_geo': les_input['vg'],
        'w_ls': les_input['wls'],
        'thl_ls': les_input['dtthl_advec'],
        'qt_ls': les_input['dtqt_advec'],
        'u_ls': les_input['dtu_advec'],
        'v_ls': les_input['dtv_advec'],
        'thl_nudge': les_input['thl'],
        'qt_nudge': les_input['qt'],
        'u_nudge': les_input['u'],
        'v_nudge': les_input['v']}

soil = {
        'z': z_soil,
        'theta_soil': era.theta_soil_mean[0,::-1],
        't_soil': era.T_soil_mean[0,::-1],
        'index_soil': soil_index,
        'root_frac': root_frac}

mht.write_NetCDF_input(
        'cabauw', 'f8', init_profiles,
        timedep_surface, timedep_ls, radiation, soil)

#
# Update namelist
#
nl = mht.read_namelist('cabauw.ini.base')

nl['grid']['ktot'] = grid.kmax
nl['grid']['zsize'] = grid.zsize
nl['time']['endtime'] = era.time_sec.max()
nl['force']['fc'] = era.fc
nl['radiation']['lon'] = settings['central_lon']
nl['radiation']['lat'] = settings['central_lat']

start = settings['start_date']
datetime_utc = '{0:04d}-{1:02d}-{2:02d} {3:02d}:{4:02d}:{5:02d}'.format(
        start.year, start.month, start.day, start.hour, start.minute, start.second)
nl['time']['datetime_utc'] = datetime_utc

# Add column locations
column_x = np.array([1300,1800,2300])
column_y = np.array([1300,1800,2300])
x,y = np.meshgrid(column_x, column_y)
nl['column']['coordinates[x]'] = list(x.flatten())
nl['column']['coordinates[y]'] = list(y.flatten())

mht.write_namelist('cabauw.ini', nl)
