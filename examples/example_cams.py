#
# This file is part of LS2D.
#
# Copyright (c) 2017-2023 Wageningen University & Research
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

# Only necessary if (LS)2D is not installed.
sys.path.append('/home/bart/meteo/models/LS2D')

import ls2d

env = {
        'system': 'arch',
        'era5_path': '/home/scratch1/meteo_data/LS2D/',
        'cdsapirc': '/home/bart/.cdsapirc'}

cams_vars = {
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
            'surface_pressure'],
        'egg4_ml': [
            'carbon_dioxide',
            'methane']
        }

# Dictionary with (LS)2D settings
settings = {
    'central_lon'   : 4.92,
    'central_lat'   : 51.97,
    'start_date'    : datetime(year=2016, month=8, day=15, hour=6),
    'end_date'      : datetime(year=2016, month=8, day=15, hour=18),
    'area_size'     : 1,
    'case_name'     : 'cabauw',
    'era5_path'     : env['era5_path'],
    'cdsapirc'      : env['cdsapirc'],
    'cams_vars'     : cams_vars,
    'era5_expver'   : 1,
    'write_log'     : False,
    'data_source'   : 'CDS',
    'ntasks'        : 1
    }

ls2d.download_cams(settings)
