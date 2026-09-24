#
# This file is part of LS2D.
#
# Copyright (c) 2017-2026 Wageningen University & Research
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

"""
Definition of the generic 3D dataset, which every backend (ERA5, ...) must provide.

Conventions:
    - SI units.
    - `level` and `level_half` from surface to top of atmosphere.
    - `pressure_level` [Pa] from surface to top of atmosphere.
    - `latitude` south to north, `longitude` west to east (regular grid, in degrees).
    - `soil_layer` from top to bottom, with coordinate `z_soil` [m, negative].
    - Attributes `central_lat` and `central_lon`: location of the column.
    - All model-specific conversions (thermodynamics, fluxes, land surface) are done by the backend,
      so that the generic code does not depend on any model constants.
"""

# Third party modules
import numpy as np

# LS2D modules
from ls2d.core.logger import logger

_d3 = ('time', 'level', 'latitude', 'longitude')
_d3h = ('time', 'level_half', 'latitude', 'longitude')
_d3p = ('time', 'pressure_level', 'latitude', 'longitude')
_d3s = ('time', 'soil_layer', 'latitude', 'longitude')
_d3s_static = ('soil_layer', 'latitude', 'longitude')
_d2 = ('time', 'latitude', 'longitude')

# name: (dims, units, long name)
required = {
    'p': (_d3, 'Pa', 'full level pressure'),
    'z': (_d3, 'm', 'full level height above surface'),
    'T': (_d3, 'K', 'absolute temperature'),
    'thl': (_d3, 'K', 'liquid water potential temperature'),
    'qt': (_d3, 'kg kg-1', 'total specific humidity'),
    'u': (_d3, 'm s-1', 'zonal wind'),
    'v': (_d3, 'm s-1', 'meridional wind'),
    'w': (_d3, 'm s-1', 'vertical wind'),
    'h2o': (_d3, 'mol mol-1', 'moisture volume mixing ratio'),
    'ph': (_d3h, 'Pa', 'half level pressure'),
    'zh': (_d3h, 'm', 'half level height above surface'),
    'phi_p': (_d3p, 'm2 s-2', 'geopotential on pressure levels'),
    'ps': (_d2, 'Pa', 'surface pressure'),
    'ts': (_d2, 'K', 'surface (skin) temperature'),
    'wth': (_d2, 'K m s-1', 'surface kinematic sensible heat flux (positive upward)'),
    'wq': (_d2, 'kg kg-1 m s-1', 'surface kinematic moisture flux (positive upward)'),
}

optional = {
    'o3': (_d3, 'mol mol-1', 'ozone volume mixing ratio'),
    'sst': (_d2, 'K', 'sea surface temperature'),
    'z0m': (_d2, 'm', 'roughness length momentum'),
    'z0h': (_d2, 'm', 'roughness length scalars'),
    'lai_low': (_d2, '-', 'leaf area index low vegetation'),
    'lai_high': (_d2, '-', 'leaf area index high vegetation'),
    'c_low_veg': (_d2, '-', 'fraction low vegetation'),
    'c_high_veg': (_d2, '-', 'fraction high vegetation'),
    'type_soil': (_d2, '-', 'soil type'),
    'type_low_veg': (_d2, '-', 'low vegetation type'),
    'type_high_veg': (_d2, '-', 'high vegetation type'),
    't_soil': (_d3s, 'K', 'soil temperature'),
    'theta_soil': (_d3s, 'm3 m-3', 'soil moisture content'),
    'root_frac_low_veg': (_d3s_static, '-', 'root fraction low vegetation'),
    'root_frac_high_veg': (_d3s_static, '-', 'root fraction high vegetation'),
}


def validate(ds):
    """
    Check if `ds` follows the generic 3D dataset definition.
    """

    for name, (dims, _, _) in required.items():
        if name not in ds:
            msg = f'Generic dataset: missing required variable "{name}"'
            logger.error(msg)
            raise ValueError(msg)
        if ds[name].dims != dims:
            msg = f'Generic dataset: "{name}" has dims {ds[name].dims}, expected {dims}'
            logger.error(msg)
            raise ValueError(msg)

    for name, (dims, _, _) in optional.items():
        if name in ds and ds[name].dims != dims:
            msg = f'Generic dataset: "{name}" has dims {ds[name].dims}, expected {dims}'
            logger.error(msg)
            raise ValueError(msg)

    for attr in ('central_lat', 'central_lon'):
        if attr not in ds.attrs:
            msg = f'Generic dataset: missing attribute "{attr}"'
            logger.error(msg)
            raise ValueError(msg)

    if np.any(np.diff(ds.latitude.values) <= 0):
        msg = 'Generic dataset: latitude should be ascending (south to north)'
        logger.error(msg)
        raise ValueError(msg)
    if np.any(np.diff(ds.pressure_level.values) >= 0):
        msg = 'Generic dataset: pressure_level should be descending (surface to top)'
        logger.error(msg)
        raise ValueError(msg)
    if np.any(ds.p.isel(level=0) < ds.p.isel(level=-1)):
        msg = 'Generic dataset: levels should go from surface to top'
        logger.error(msg)
        raise ValueError(msg)
