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
from datetime import datetime
import sys

# Third party modules
import numpy as np

# LS2D & custom modules
# The next line is only needed if you did not install (LS)2D using PyPI:
#sys.path.append('/home/bart/meteo/models/LS2D')

import ls2d
import microhh_ls2d_tools as mlt

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
    'data_source' : 'CDS'
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

# Interpolate ERA5 variables and forcings onto LES grid.
# In addition, `get_les_input` returns additional variables needed to init LES.
les_input = era.get_les_input(grid.z)

# Remove top ERA5 level, to ensure that pressure stays
# above the minimum reference pressure in RRTMGP
les_input = les_input.sel(lay=slice(0,135), lev=slice(0,136))

#
# MicroHH specific initialisation
#
# Settings:
# ------------------
# Root fraction coefficients (see IFS documentation):
a_r = 10.739
b_r = 2.608

# Nudging time scale atmosphere
tau_nudge = 10800
# ------------------

# Fixed background concentrations RRTMGP:
co2 = 348.e-6
ch4 = 1650.e-9
n2o = 306.e-9
n2  = 0.7808
o2  = 0.2095

# Mean radiation profiles on LES grid:
qt_les = les_input['qt'].mean(axis=0)
o3_les = les_input['o3'].mean(axis=0)

# Conversion moisture from mass to volume mixing ratio
h2o_les = qt_les / (ep - ep*qt_les)

## Soil
z_soil     = les_input['zs'][::-1].values
soil_index = les_input['type_soil'].values-1
soil_index = np.ones_like(z_soil)*soil_index
root_frac  = mlt.calc_root_frac(z_soil, a_r, b_r)

# Nudge factor
nudge_fac = np.ones(grid.kmax) / tau_nudge

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
        'o3': o3_les*1e-6,
        'h2o': h2o_les}

radiation  = {
        'z_lay': les_input['z_lay'].mean(axis=0),
        'z_lev': les_input['z_lev'].mean(axis=0),
        'p_lay': les_input['p_lay'].mean(axis=0),
        'p_lev': les_input['p_lev'].mean(axis=0),
        't_lay': les_input['t_lay'].mean(axis=0),
        't_lev': les_input['t_lev'].mean(axis=0),
        'o3': les_input['o3_lay'].mean(axis=0)*1e-6,
        'h2o': les_input['h2o_lay'].mean(axis=0),
        'co2': co2,
        'ch4': ch4,
        'n2o': n2o,
        'n2': n2,
        'o2': o2}

timedep_surface = {
        'time_surface': les_input['time_sec'],
        'p_sbot': les_input['ps'] }

timedep_ls = {
        'time_ls': les_input['time_sec'],
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
        'theta_soil': les_input['theta_soil'][0,::-1],
        't_soil': les_input['t_soil'][0,::-1],
        'index_soil': soil_index,
        'root_frac': root_frac}

mlt.write_netcdf_input(
        'cabauw', 'f8', init_profiles,
        timedep_surface, timedep_ls, radiation, soil)

#
# Update namelist
#
nl = mlt.read_namelist('cabauw.ini.base')

nl['grid']['ktot'] = grid.kmax
nl['grid']['zsize'] = grid.zsize
nl['time']['endtime'] = float(les_input['time_sec'][-1])
nl['force']['fc'] = les_input.attrs['fc']
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

mlt.write_namelist('cabauw.ini', nl)
