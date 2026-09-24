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

# Third party modules
import numpy as np

# Depth soil layers (m), top to bottom.
z_soil = np.array([-0.035, -0.175, -0.64, -1.945])

# Root-fraction coefficients per vegetation type, from ECMWF documentation.
ar = np.array([
    5.558, 10.739, 6.706, 7.066, 5.99, 7.344, 8.235, 4.372, 8.992, 5.558,
    4.372, -1.0, 7.344, -1.0, -1.0, 6.326, 6.326, 4.453, 4.453, -1.0,
])
br = np.array([
    2.614, 2.608, 2.175, 1.953, 1.955, 1.303, 1.627, 0.978, 8.992, 2.614,
    0.978, -1.0, 1.303, -1.0, -1.0, 1.567, 1.567, 1.631, 1.631, -1.0,
])


def root_fraction(veg_type):
    """
    Root fraction per soil layer (top to bottom) for HTESSEL vegetation type(s) `veg_type` (Fortran indexing).
    Returns array with shape (soil layer, *veg_type.shape).
    """

    # Soil layer interfaces, top to bottom, with the full levels halfway.
    zi = np.zeros(z_soil.size + 1)
    for k in range(z_soil.size):
        zi[k + 1] = 2 * z_soil[k] - zi[k]

    i = np.asarray(veg_type) - 1
    a, b = ar[i][None], br[i][None]
    shape = (-1,) + (1,) * i.ndim
    zt, zb = zi[:-1].reshape(shape), zi[1:].reshape(shape)

    rf = 0.5 * (np.exp(a * zt) + np.exp(b * zt) - np.exp(a * zb) - np.exp(b * zb))

    # Bottom layer gets the remainder.
    rf[-1] = 1 - rf[:-1].sum(axis=0)
    return rf
