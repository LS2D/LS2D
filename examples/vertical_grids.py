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

# Python modules.
import sys

# Third party modules.
import matplotlib.pyplot as pl
import numpy as np

# LS2D modules.
sys.path.append('/home/bart/meteo/models/LS2D')
import ls2d

# Equidistant grid.
grid1 = ls2d.grid.Grid_equidist(kmax=28, dz0=20)

# Linearly stretched grid.
grid2 = ls2d.grid.Grid_linear_stretched(kmax=24, dz0=15, alpha=0.03)

# Manually stretched grid.
heights = [0,100,10000]
alpha = [1.02, 1.05]
grid3 = ls2d.grid.Grid_stretched_manual(kmax=24, dz0=12.5, heights=heights, factors=alpha)

# Smoothly stretched grids (from Chiel's DNSs).
grid4 = ls2d.grid.Grid_stretched(kmax=32, dz0=10, nloc1=15, nbuf1=5, dz1=25)

# Plot!
pl.figure()
pl.plot(grid1.dz, grid1.z, '-o', ms=5, label='Grid_equidist')
pl.plot(grid2.dz, grid2.z, '-o', ms=5, label='Grid_linear_stretched')
pl.plot(grid3.dz, grid3.z, '-o', ms=5, label='Grid_stretched_manual')
pl.plot(grid4.dz, grid4.z, '-o', ms=5, label='Grid_stretched')
pl.xlabel(r'$\Delta$ z (m)')
pl.ylabel(r'$z$ (m)')
pl.legend()
