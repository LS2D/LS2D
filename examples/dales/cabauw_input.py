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
from collections import OrderedDict as odict
import sys

# Third party modules
import numpy as np

# LS2D & custom modules
sys.path.append('/home/bart/meteo/models/LS2D')
sys.path.append('/Users/bart/meteo/models/LS2D')
import ls2d
import dales_ls2d_tools as dlt

#
# Download ERA5 and generate LES initialisation and forcings
#
settings = {
    'central_lat' : 51.971,
    'central_lon' : 4.927,
    'area_size'   : 1,
    'case_name'   : 'cabauw',
    'era5_path'   : '/home/scratch1/meteo_data/LS2D/',
    #'era5_path'   : '/Users/bart/meteo/data/ERA5/LS2D',
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
era.calculate_forcings(n_av=0, method='2nd')

# Define vertical grid LES:
grid = ls2d.grid.Grid_linear_stretched(kmax=176, dz0=20, alpha=0.009)
#grid.plot()

# Interpolate ERA5 variables and forcings onto LES grid.
# In addition, `get_les_input` returns additional variables needed to init LES.
les_input = era.get_les_input(grid.z)

#
# DALES specific initialisation.
#
# Settings:
expnr = 1
tau_nudge = 10800    # Nudging time scale atmosphere
init_tke = 0.1       # Initial SGS-TKE

docstring = '(LS)2D case Cabauw: {} to {}'.format(
        settings['start_date'].isoformat(), settings['end_date'].isoformat())

#
# Write initial profiles to `prof.inp.expnr`.
#
output = odict([
        ('z (m)',        grid.z),
        ('thl (K)',      les_input.thl[0,:].values),
        ('qt (kg kg-1)', les_input.qt[0,:].values),
        ('u (m s-1)',    les_input.u[0,:].values),
        ('v (m s-1)',    les_input.v[0,:].values),
        ('tke (m2 s-2)', np.ones(grid.kmax)*init_tke)])

dlt.write_profiles(
        'prof.inp.{0:03d}'.format(expnr), output, grid.kmax, docstring)

#
# Write initial scalar profiles to `scalar.inp.expnr`.
#
zero = np.zeros(grid.kmax)
output = odict([
        ('z (m)',        grid.z),
        ('qr (kg kg-1)', zero),
        ('nr (kg kg-1)', zero)])

dlt.write_profiles(
        'scalar.inp.{0:03d}'.format(expnr), output, grid.kmax, docstring)

#
# Write large-scale forcings to `ls_flux.inp.expnr`.
#
output_sfc = odict([
        ('time',   les_input.time_sec.values),
        ('p_s',    les_input.ps.values),
        ('T_s',    np.zeros_like(les_input.time_sec)),      # Not sure if this works..
        ('qt_s',   np.zeros_like(les_input.time_sec))])

output_ls = odict([
        ('time',   les_input.time_sec.values),
        ('z',      grid.z),
        ('ug',     les_input.ug.values),
        ('vg',     les_input.vg.values),
        ('wls',    les_input.wls.values),
        ('dqtdt',  les_input.dtqt_advec.values),
        ('dthldt', les_input.dtthl_advec.values),
        ('dudt',   les_input.dtu_advec.values),
        ('dvdt',   les_input.dtv_advec.values)])

dlt.write_forcings(
        'ls_flux.inp.{0:03d}'.format(expnr), output_sfc, output_ls, docstring)

#
# Write nudging profiles to `nudge.inp.expnr`.
#
output = odict([
        ('z (m)',        grid.z),
        ('factor (-)',   np.ones_like(les_input.u.values)),
        ('u (m s-1)',    les_input.u.values),
        ('v (m s-1)',    les_input.v.values),
        ('w (m s-1)',    np.zeros_like(les_input.u.values)),
        ('thl (K)',      les_input.thl.values),
        ('qt (kg kg-1)', les_input.qt.values)])

dlt.write_time_profiles(
        'nudge.inp.{0:03d}'.format(expnr), les_input.time_sec.values, output, grid.kmax, docstring)

#
# Also create non-time dependent file (lscale.inp), required by DALES (why?)
#
zero = np.zeros_like(grid.z)

output = odict([
        ('height', grid.z),
        ('ug', zero),
        ('vg', zero),
        ('wfls', zero),
        ('dqtdxls', zero),
        ('dqtdyls', zero),
        ('dqtdtls', zero),
        ('dthldt', zero)])

dlt.write_profiles('lscale.inp.{0:03d}'.format(expnr), output, grid.kmax, docstring)

#
# Write radiation background profiles to `backrad.inp.expnr`.
#
dlt.create_backrad(
        les_input['p_lay'].mean(axis=0),
        les_input['t_lay'].mean(axis=0),
        les_input['h2o_lay'].mean(axis=0))
