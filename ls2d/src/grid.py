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

import matplotlib.pyplot as pl
import numpy as np
import sys

#
# Vertical grids
#
class _Grid:
    def __init__(self, kmax, dz0):
        self.kmax = kmax
        self.dz0  = dz0

        self.z = np.zeros(kmax)
        self.zh = np.zeros(kmax+1)
        self.dz = np.zeros(kmax)
        self.zsize = None

    def plot(self, logx=False, logy=False):

        fig=pl.figure()
        ax=pl.subplot(121)
        pl.title(r'$z_\mathrm{{size}}$ = {0:.1f} m'.format(self.zsize), loc='left')
        pl.plot(self.dz, self.z, 'k-x')
        pl.xlabel(r'$\Delta z$ (m)')
        pl.ylabel(r'$z$ (m)')
        pl.grid()

        ax=pl.subplot(122)
        pl.plot(np.arange(self.kmax)+1, self.z, 'k-x')
        pl.xlabel(r'Level (-)')
        pl.ylabel(r'$z$ (m)')
        pl.grid()

        if logx:
            for ax in fig.axes:
                ax.set_xscale('log')
        if logy:
            for ax in fig.axes:
                ax.set_yscale('log')

        pl.tight_layout()


class Grid_equidist(_Grid):
    def __init__(self, kmax, dz0):
        _Grid.__init__(self, kmax, dz0)

        self.zsize = kmax * dz0
        self.z[:]  = np.arange(dz0/2, self.zsize, dz0)
        self.zh[:] = np.arange(0, self.zsize+0.1, dz0)
        self.dz[:] = dz0


class Grid_stretched(_Grid):
    def __init__(self, kmax, dz0, nloc1, nbuf1, dz1, nloc2=None, nbuf2=None, dz2=None):
        _Grid.__init__(self, kmax, dz0)

        double_stretched = nloc2 is not None and nbuf2 is not None and dz2 is not None

        dn = 1./kmax
        n = np.linspace(dn, 1.-dn, kmax)

        nloc1 *= dn
        nbuf1 *= dn

        if double_stretched:
            nloc2 *= dn
            nbuf2 *= dn

        dzdn1 = dz0/dn
        dzdn2 = dz1/dn

        if double_stretched:
            dzdn3  = dz2/dn
            dzdn = dzdn1 + 0.5*(dzdn2-dzdn1)*(1. + np.tanh((n-nloc1)/nbuf1)) \
                         + 0.5*(dzdn3-dzdn2)*(1. + np.tanh((n-nloc2)/nbuf2))
        else:
            dzdn = dzdn1 + 0.5*(dzdn2-dzdn1)*(1. + np.tanh((n-nloc1)/nbuf1))

        self.dz[:] = dzdn*dn

        stretch = np.zeros(self.dz.size)
        self.z[0]  = 0.5*self.dz[0]
        stretch[0] = 1.

        for k in range(1, self.kmax):
              self.z [k] = self.z[k-1] + 0.5*(self.dz[k-1]+self.dz[k])
              stretch[k] = self.dz[k]/self.dz[k-1]

        self.zsize = self.z[kmax-1] + 0.5*self.dz[kmax-1]

        self.zh[1:-1] = self.z[1:] - self.z[:-1]
        self.zh[-1] = self.zsize



class Grid_linear_stretched(_Grid):
    def __init__(self, kmax, dz0, alpha):
        _Grid.__init__(self, kmax, dz0)

        self.dz[:] = dz0 * (1 + alpha)**np.arange(kmax)
        self.zh = np.zeros(kmax+1)
        self.zh[1:]= np.cumsum(self.dz)
        self.z[:] = 0.5 * (self.zh[1:] + self.zh[:-1])
        self.zsize = self.zh[-1]


class Grid_stretched_manual(_Grid):
    def __init__(self, kmax, dz0, heights, factors):
        _Grid.__init__(self, kmax, dz0)

        self.z[0]  = dz0/2.
        self.dz[0] = dz0

        def index(z, goal):
            return np.where(z-goal>0)[0][0]-1

        for k in range(1, kmax):
            self.dz[k] = self.dz[k-1] * factors[index(heights, self.z[k-1])]
            self.z[k] = self.z[k-1] + self.dz[k]

        self.zsize = self.z[kmax-1] + 0.5*self.dz[kmax-1]

        self.zh[1:-1] = self.z[1:] - self.z[:-1]
        self.zh[-1] = self.zsize
