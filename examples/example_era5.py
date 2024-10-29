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
import subprocess
import shutil
import sys,os

# Third party modules
import matplotlib.pyplot as pl
import numpy as np

import ls2d

settings = {
    'central_lat' : 51.97,
    'central_lon' : 4.93,
    'area_size'   : 1,
    'case_name'   : 'cabauw',
    'era5_path'   : '/home/scratch1/bart/LS2D_ERA5/',
    'era5_expver' : 1,   # 1=normal ERA5, 5=ERA5 near-realtime
    'start_date'  : datetime(year=2016, month=8, day=15, hour=6),
    'end_date'    : datetime(year=2016, month=8, day=15, hour=18),
    'write_log'   : False,
    'data_source' : 'CDS'
    }

# Download required ERA5 files:
ls2d.download_era5(settings)

# Read ERA5 data, and calculate derived properties (thl, etc.):
era = ls2d.Read_era5(settings)

# Calculate large-scale forcings:
# `n_av` is the number of ERA5 gridpoints (+/-) over which
# the ERA5 variables and forcings are averaged.
era.calculate_forcings(n_av=1, method='2nd')

# Interpolate ERA5 to fixed height grid:
z = np.arange(10, 5000, 20).astype(float)
les_input = era.get_les_input(z)

# `les_input` is an xarray.Dataset, which can easily be save to NetCDF:
les_input.to_netcdf('ls2d_era5.nc')

# Plot variables as example:
nrow = 5
ncol = 2
sp = 1

pl.figure(figsize=(8,8))
pl.subplot(nrow, ncol, sp); sp+=1
pl.plot(les_input['thl'].T, les_input['z'])
pl.xlabel(r'$\theta_l$ (K)')

pl.subplot(nrow, ncol, sp); sp+=1
pl.plot(les_input['dtthl_advec'].T*3600, les_input['z'])
pl.xlabel(r'$\partial \theta_l/\partial t$ advec (K h$^{-1}$)')

pl.subplot(nrow, ncol, sp); sp+=1
pl.plot(les_input['qt'].T*1e3, les_input['z'])
pl.xlabel(r'$q_t$ (g kg$^{-1}$)')

pl.subplot(nrow, ncol, sp); sp+=1
pl.plot(les_input['dtqt_advec'].T*3600000, les_input['z'])
pl.xlabel(r'$\partial q_t/\partial t$ advec (g kg$^{-1}$ h$^{-1}$)')

pl.subplot(nrow, ncol, sp); sp+=1
pl.plot(les_input['u'].T, les_input['z'])
pl.xlabel(r'$u$ (m s$^{-1}$)')

pl.subplot(nrow, ncol, sp); sp+=1
pl.plot(les_input['dtu_advec'].T*3600, les_input['z'])
pl.xlabel(r'$\partial u/\partial t$ advec (m s$^{-1}$ h$^{-1}$)')

pl.subplot(nrow, ncol, sp); sp+=1
pl.plot(les_input['v'].T, les_input['z'])
pl.xlabel(r'$v$ (m s$^{-1}$)')

pl.subplot(nrow, ncol, sp); sp+=1
pl.plot(les_input['dtv_advec'].T*3600, les_input['z'])
pl.xlabel(r'$\partial v/\partial t$ advec (m s$^{-1}$ h$^{-1}$)')

pl.subplot(nrow, ncol, sp); sp+=1
pl.plot(les_input['ug'].T, les_input['z'])
pl.xlabel(r'$u_g$ (m s$^{-1}$)')

pl.subplot(nrow, ncol, sp); sp+=1
pl.plot(les_input['vg'].T, les_input['z'])
pl.xlabel(r'$v_g$ (m s$^{-1}$)')

pl.tight_layout()
