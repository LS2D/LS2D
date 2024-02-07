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
import pkg_resources
import os

# Third party modules
import numpy as np

class IFS_tools:
    """
    Various tools to calculate e.g. properties of the vertical IFS/ERA grid,
    or the thermodynamics as used by IFS. Wrapped in a class, to prevent
    mixing up differences in methods/constants/.. between IFS and other models
    """
    def __init__(self, grid_def='L137'):

        # Constants (see IFS part IV, chapter 12)
        self.grav  = 9.80665
        self.Rd    = 287.0597
        self.Rv    = 461.5250
        self.eps   = self.Rv/self.Rd-1.
        self.cpd   = 1004.7090
        self.Lv    = 2.5008e6

        # Read the table with the vertical grid properties/parameters
        # From: https://www.ecmwf.int/en/forecasts/documentation-and-support/137-model-levels
        if grid_def == 'L137':
            path = pkg_resources.resource_filename(__name__, 'L137_grid.txt')
        elif grid_def == 'L60':
            path = pkg_resources.resource_filename(__name__, 'L60_grid.txt')
        else:
            sys.exit(f'Unknow grid definition \"{grid_def}\".')
        f = np.loadtxt(path)

        # Half and full level number
        self.nh  = f[:,0]
        self.nf  = 0.5 * (self.nh[1:] + self.nh[:-1])

        # Reverse all arrays (::-1) from top-to-bottom to bottom-to-top
        self.a   = f[ :,1][::-1]       # a-coefficient pressure calculation
        self.b   = f[ :,2][::-1]       # b-coefficient pressure calculation

        # Values below are only valid for reference atmosphere with ps=101325 Pa
        self.ph  = f[ :,3][::-1]*100   # Reference half level pressure (Pa)
        self.pf  = f[1:,4][::-1]*100   # Referente full level pressure (Pa)
        self.gpa = f[1:,5][::-1]       # Reference geopotential height (m)
        self.gma = f[1:,6][::-1]       # Reference geometric altitude (m)
        self.T   = f[1:,7][::-1]       # Reference temperature (K)
        self.rho = f[1:,8][::-1]       # Reference density (kg m-3)

    def calc_half_level_pressure(self, ps):
        """
        Calculate half level pressure
        See IFS part III, eq. 2.11
        Equation: p = a + b * ps
        Top value is set to a small non-zero number to prevent div-by-0's
        Keyword arguments: 
            ps -- surface pressure (Pa)
        """

        ph = self.a + self.b * ps
        ph[-1] = 0.34     # Chosen to match IFS values for standard atmosphere
        return ph

    def calc_full_level_pressure(self, ps):
        """
        Calculate full level pressure as a linear interpolation of the half level pressure
        Equation: p = a + b * ps
        See IFS part III, eq. 2.11
        Keyword arguments: 
            ps -- surface pressure (Pa)
        """

        p = self.calc_half_level_pressure(ps)
        return 0.5 * (p[1:] + p[:-1])

    def calc_half_level_Zg(self, ph, Tv):
        """
        Calculate half level geopotential height
        Equation: sums dZg = -Rd / g * Tv * ln(p+ / p-)
        See IFS part III, eq. 2.20-2.21
        Keyword arguments: 
            ph -- half level pressure (Pa)
            Tv -- full level virtual temperature (K) 
    `    """

        pfrac = ph[1:] / ph[:-1]
        dZg   = -self.Rd * Tv * np.log(pfrac) / self.grav
        Zg    = np.cumsum(dZg) 
        Zg    = np.insert(Zg, 0, 0) 

        return Zg

    def calc_full_level_Zg(self, ph, Tv):
        """
        Calculate full level geopotential height
        Equation: sums dZg = -Rd / g * Tv * ln(p+ / p-)
        See IFS part III, eq. 2.20-2.21
        Keyword arguments: 
            ph -- half level pressure (Pa)
            Tv -- full level virtual temperature (K) 
        """

        Zg = self.calc_half_level_Zg(ph, Tv)
        return 0.5 * (Zg[1:] + Zg[:-1])

    def calc_virtual_temp(self, T, qv, ql=0, qi=0, qr=0, qs=0):
        """
        Calculate the virtual temperature
        Equation: Tv = T * ([Rv/Rd-1]*qv - ql - qi - qr - qs)
        See IFS part IV, eq. 12.6
        Keyword arguments:
            T -- absolute temperature (K)
            q* -- specific humidities (kg kg-1):
                qv = vapor
                ql = liquid (optional)
                qi = ice    (optional)
                qr = rain   (optional)
                qs = snow   (optional)
        """

        return T * (1+self.eps*qv - ql - qi - qr - qs)

    def calc_exner(self, p):
        return (p/1e5)**(self.Rd/self.cpd)

    def validate(self):
        """
        Compare calculated grid properties vs. the tabulated
        ones in the ECMWF table
        """

        ps = 101325     # Surface pressure (Pa)

        # Calculate pressure at full and half levels
        pf = self.calc_full_level_pressure(ps)
        ph = self.calc_half_level_pressure(ps)

        rhof = pf / (self.Rd * self.T)  # assumes T=Tv

        # Calculate geopotential height on full and half levels
        Zgf = self.calc_full_level_Zg(ph, self.T)
        Zgh = self.calc_half_level_Zg(ph, self.T)

        # Compare calculated values with reference data from ECMWF table
        pl.close('all')
        pl.figure(figsize=[14,8])

        pl.subplot(231)
        pl.plot(self.nh, self.ph, label='Table half')
        pl.plot(self.nh, ph, label='Calculated half', marker='x', linestyle='')
        pl.plot(self.nf+20, self.pf, label='Table full (+offs)')
        pl.plot(self.nf+20, pf, label='Calculated full (+offs)', marker='x', linestyle='')
        pl.legend(frameon=False, loc='best')
        pl.xlabel('Level (-)')
        pl.ylabel('p (Pa)')

        pl.subplot(234)
        pl.title('rel. diff. p', loc='left')
        pl.semilogy(self.nh, np.abs((self.ph-ph)/self.ph), label='half')
        pl.semilogy(self.nf, np.abs((self.pf-pf)/self.pf), label='full')
        pl.legend(frameon=False, loc='best')
        pl.xlabel('Level (-)')
        pl.ylabel('(-)')

        pl.subplot(232)
        pl.plot(self.nf, self.rho, label='Table')
        pl.plot(self.nf, rhof, label='Calculated', marker='x', linestyle='')
        pl.legend(frameon=False, loc='best')
        pl.xlabel('Level (-)')
        pl.ylabel('rho (kg m-3)')

        pl.subplot(235)
        pl.title('rel. diff. rho', loc='left')
        pl.semilogy(self.nf, np.abs((self.rho-rhof)/self.rho))
        pl.xlabel('Level (-)')
        pl.ylabel('(-)')

        pl.subplot(233)
        pl.plot(self.nf, self.gpa, label='Table')
        pl.plot(self.nf, Zgf, label='Calculated', marker='x', linestyle='')
        pl.legend(frameon=False, loc='best')
        pl.xlabel('Level (-)')
        pl.ylabel('Geop. height (m)')

        pl.subplot(236)
        pl.title('rel. diff. z', loc='left')
        pl.semilogy(self.nf, np.abs((self.gpa-Zgf)/self.gpa))
        pl.xlabel('Level (-)')
        pl.ylabel('(-)')

        pl.tight_layout()
        pl.savefig('IFS_tools_validation.pdf')


if __name__ == "__main__":
    """ Test / example, only executed if script is called directly """

    import matplotlib.pyplot as pl
    pl.ion()
    pl.close('all')

    ifs_tools = IFS_tools('L60')

    # Validate the IFS_tools conversion from surface pressure to model level pressure, height, etc.
    ifs_tools.validate()
