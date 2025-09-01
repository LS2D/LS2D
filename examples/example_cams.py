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

import matplotlib.pyplot as pl
import numpy as np

import ls2d

env = {
    'system': 'arch',
    'era5_path': '/home/scratch1/bart/LS2D_ERA5/',
    'cams_path': '/home/scratch1/bart/LS2D_CAMS/',
    'cdsapirc': '/home/bart/.cdsapirc'}

cams_eac4_vars = {
    'eac4_ml': [
        'dust_aerosol_0.03-0.55um_mixing_ratio',
        'dust_aerosol_0.55-0.9um_mixing_ratio',
        'dust_aerosol_0.9-20um_mixing_ratio',
        'hydrophilic_black_carbon_aerosol_mixing_ratio',
        'hydrophilic_organic_matter_aerosol_mixing_ratio',
        'hydrophobic_black_carbon_aerosol_mixing_ratio',
        'hydrophobic_organic_matter_aerosol_mixing_ratio',
        'sea_salt_aerosol_0.03-0.5um_mixing_ratio',
        'sea_salt_aerosol_0.5-5um_mixing_ratio',
        'sea_salt_aerosol_5-20um_mixing_ratio',
        'specific_humidity',
        'sulphate_aerosol_mixing_ratio',
        'temperature'],
    'eac4_sfc': [
        'surface_pressure']
    }

# NOTE: EGG4 provides several variables defined on model levels,
#       but some of them (surface pressure, geopotential, etc.)
#       are only available at model level 1. These cannot
#       be downloaded together with standard model-level variables.
#       To retrieve them, use the `egg4_sl` variable group instead.
cams_egg4_vars = {
    'egg4_ml': [
        'carbon_dioxide',
        'methane',
        'temperature',
        'specific_humidity'],
    'egg4_sl': [
        'logarithm_of_surface_pressure']
    }


# Dictionary with (LS)2D settings
settings = {
    'central_lon'   : 4.92,
    'central_lat'   : 51.97,
    'start_date'    : datetime(year=2016, month=8, day=15, hour=6),
    'end_date'      : datetime(year=2016, month=8, day=15, hour=18),
    'area_size'     : 2,
    'case_name'     : 'cabauw',
    'cams_path'     : env['cams_path'],
    'cdsapirc'      : env['cdsapirc'],
    'era5_expver'   : 1,
    'write_log'     : False,
    'data_source'   : 'CDS',
    'ntasks'        : 1
    }

#settings = {
#    'central_lon'   : -57.7,
#    'central_lat'   : 13.3,
#    'start_date'    : datetime(year=2022, month=1, day=2, hour=6),
#    'end_date'      : datetime(year=2022, month=1, day=2, hour=18),
#    'area_size'     : 2,
#    'case_name'     : 'barbados',
#    'cams_path'     : env['cams_path'],
#    'cdsapirc'      : env['cdsapirc'],
#    'era5_expver'   : 1,
#    'write_log'     : False,
#    'data_source'   : 'CDS',
#    'ntasks'        : 1
#    }


# Download data.
ls2d.download_cams(settings, variables=cams_eac4_vars, grid=0.25)
ls2d.download_cams(settings, variables=cams_egg4_vars, grid=0.25)

# Read, average over 3x3 grid points, and interpolate on LES grid.
cams_eac = ls2d.Read_cams(settings, variables=cams_eac4_vars)
cams_egg = ls2d.Read_cams(settings, variables=cams_egg4_vars)

z = np.arange(10, 5000, 20).astype(float)
les_eac = cams_eac.get_les_input(z, n_av=1)
les_egg = cams_egg.get_les_input(z, n_av=1)

les_eac.to_netcdf('ls2d_eac4_cams.nc')
les_egg.to_netcdf('ls2d_egg4_cams.nc')

# Quick plot.
for les_input in [les_eac, les_egg]:

    pl.figure(figsize=(10,8), layout='constrained')
    ncol = 5
    nrow = 6
    sp = 1

    for name, da in les_input.data_vars.items():
        pl.subplot(nrow, ncol, sp); sp+=1
        if 'lay' in name:
            pl.plot(da[0,:], les_input.z_lay[0,:])
        else:
            pl.plot(da[0,:], les_input.z)
        pl.xlabel(name)