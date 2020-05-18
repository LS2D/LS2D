import matplotlib.pyplot as pl
import numpy as np

#
# Vertical grids
#
class Grid:
    def __init__(self, kmax, dz0):
        self.kmax = kmax
        self.dz0  = dz0

        self.z = np.zeros(kmax)
        self.dz = np.zeros(kmax)
        self.zsize = None

    def plot(self):
        pl.figure()
        pl.title(r'$z_\mathrm{{size}}$ = {0:.1f} m'.format(self.zsize), loc='left')
        pl.plot(self.dz, self.z, '-x')
        pl.xlabel(r'$\Delta z$ (m)')
        pl.ylabel(r'$z$ (m)')


class Grid_equidist(Grid):
    def __init__(self, kmax, dz0):
        Grid.__init__(self, kmax, dz0)

        self.zsize = kmax * dz0
        self.z[:]  = np.arange(dz0/2, self.zsize, dz0)
        self.dz[:] = dz0


class Grid_stretched(Grid):
    def __init__(self, kmax, dz0, nloc1, nbuf1, dz1, nloc2=None, nbuf2=None, dz2=None):
        Grid.__init__(self, kmax, dz0)

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

        if double_stretched:
            dzdn = dzdn1 + 0.5*(dzdn2-dzdn1)*(1. + np.tanh((n-nloc1)/nbuf1)) + 0.5*(dzdn3-dzdn2)*(1. + np.tanh((n-nloc2)/nbuf2))
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


class Grid_linear_stretched(Grid):
    def __init__(self, kmax, dz0, alpha):
        Grid.__init__(self, kmax, dz0)

        self.dz[:] = dz0 * (1 + alpha)**np.arange(kmax)
        zh         = np.zeros(kmax+1)
        zh[1:]     = np.cumsum(self.dz)
        self.z[:]  = 0.5 * (zh[1:] + zh[:-1])
        self.zsize = zh[-1]


if __name__ == '__main__':
    """
    Just for testing...
    """
    import matplotlib.pyplot as pl
    pl.close('all'); pl.ion()

    ktot = 224
    dz0  = 20

    equidist   = Grid_equidist(ktot, dz0)
    linear     = Grid_linear_stretched(ktot, dz0, 0.01)
    stretched1 = Grid_stretched(ktot, dz0, 90, 20, 150)
    stretched2 = Grid_stretched(ktot, dz0, 100, 20, 100, 210, 10, 500)

    pl.figure()
    pl.plot(equidist.dz, equidist.z, '-x', label='equidistant')
    pl.plot(linear.dz, linear.z, '-x', label='linear')
    pl.plot(stretched1.dz, stretched1.z, '-x', label='stretched-single')
    pl.plot(stretched2.dz, stretched2.z, '-x', label='stretched-double')
    pl.legend()
    pl.xlabel('dz (m)')
    pl.ylabel('z (m)')
