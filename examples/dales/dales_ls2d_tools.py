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

from collections import OrderedDict as odict
import matplotlib.pyplot as pl
import netCDF4 as nc4
import numpy as np
import datetime

# ---------------------------
# "Private" help functions
# ---------------------------
def _get_or_default(dict, name, shape, default_value):
    if name in dict:
        return dict[name]
    else:
        print(' - No input for \"{}\", defaulting values at zero'.format(name))
        return default_value * np.ones(shape)

# ---------------------------
# Function to write the DALES input
# ---------------------------
def write_profiles(file_name, variables, nlev, docstring=''):
    """
    Write the prof.inp.xxx input profiles for DALES
    """

    print(' - Saving {}'.format(file_name))

    f = open(file_name, 'w')

    # Write header (description file)
    if docstring == '':
        f.write('DALES\n')
    else:
        f.write('{}\n'.format(docstring))

    # Write header (column names)
    for var in variables.keys():
        f.write('{0:^17s} '.format(var))
    f.write('\n')

    # Write data
    for k in range(nlev):
        for var in variables.keys():
            f.write('{0:+1.10E} '.format(variables[var][k]))
        f.write('\n')

    f.close()


def write_time_profiles(file_name, time, variables, nlev, docstring=''):
    """
    Write time varying input profiles for DALES
    """

    print(' - Saving {}'.format(file_name))

    f = open(file_name, 'w')

    # Write header (description file)
    if docstring == '':
        f.write('DALES\n')
    else:
        f.write('{}\n'.format(docstring))

    # Write time dependent profiles
    for t in range(time.size):

        f.write('\n')

        # Write header (column names)
        for var in variables.keys():
            f.write('{0:^17s} '.format(var))
        f.write('\n')

        # Write time
        f.write('# {0:1.8E}\n'.format(time[t]))

        # Write data
        for k in range(nlev):
            for var in variables.keys():
                if len(variables[var].shape) == 1:
                    f.write('{0:+1.10E} '.format(variables[var][k]))
                else:
                    f.write('{0:+1.10E} '.format(variables[var][t,k]))
            f.write('\n')

    f.close()


def write_dummy_forcings(file_name, n_scalars, z, docstring):
    """
    Write dummy forcings
    """

    print('Saving {}'.format(file_name))

    f = open(file_name, 'w')

    # Write header (description file)
    if docstring == '':
        f.write('DALES\n\n')
    else:
        f.write('{}\n\n'.format(docstring))

    # Surface fluxes (zero)
    f.write('{0:^15s} '.format('time'))
    for i in range(n_scalars):
        f.write('{0:>10s}{1:<8d}'.format('sv', i+1))
    f.write('\n')

    for time in [0,1e6]:
        f.write('{0:+1.10E} '.format(time))
        for i in range(n_scalars):
            f.write('{0:+1.10E} '.format(0))
        f.write('\n')

    # Atmospheric forcings
    f.write('\n')

    for time in [0,1e6]:
        f.write('# {0:+1.10E}\n'.format(time))
        for k in range(z.size):
            f.write('{0:+1.10E} '.format(z[k]))
            for i in range(n_scalars):
                f.write('{0:+1.10E} '.format(0))
            f.write('\n')

    f.close()


def write_forcings(file_name, timedep_sfc, timedep_atm, docstring=''):
    """
    Write the ls_flux.inp.xxx files
    """

    print(' - Saving {}'.format(file_name))


    f = open(file_name, 'w')

    # Always write something; DALES expects three line header
    if docstring == '':
        f.write('DALES time dependent input\n')
    else:
        f.write('{}\n'.format(docstring))

    # Write surface variables
    f.write('{0:^15s} {1:^15s} {2:^15s} {3:^15s} {4:^15s} {5:^15s}\n'\
        .format('time', 'wthl_s', 'wqt_s', 'T_s', 'qt_s', 'p_s'))
    f.write('{0:^15s} {1:^15s} {2:^15s} {3:^15s} {4:^15s} {5:^15s}\n'\
        .format('(s)', '(K m s-1)', '(kg kg-1 m s-1)', '(K)', '(kg kg-1)', '(Pa)'))

    if timedep_sfc is None:
        # Write a large initial time, so DALES will disable the surface timedep
        f.write('{0:+1.8E} {1:+1.8E} {2:+1.8E} {3:+1.8E} {4:+1.8E} {5:+1.8E}\n'\
            .format(1e16, -1, -1, -1, -1, -1))
    else:
        nt    = timedep_sfc['time'].size
        time  = _get_or_default(timedep_sfc, 'time',  nt, 0)
        wthls = _get_or_default(timedep_sfc, 'wthl_s',nt, 0)
        wqts  = _get_or_default(timedep_sfc, 'wqt_s', nt, 0)
        Ts    = _get_or_default(timedep_sfc, 'T_s',   nt, 0)
        qts   = _get_or_default(timedep_sfc, 'qt_s',  nt, 0)
        ps    = _get_or_default(timedep_sfc, 'p_s',   nt, 0)

        for t in range(nt):
            f.write('{0:+1.8E} {1:+1.8E} {2:+1.8E} {3:+1.8E} {4:+1.8E} {5:+1.8E}\n'\
                .format(time[t], wthls[t], wqts[t], Ts[t], qts[t], ps[t]))

    if timedep_atm is not None:
        time = timedep_atm['time']
        z    = timedep_atm['z']
        nt   = time.size
        nlev = z.size

        ug   = _get_or_default(timedep_atm, 'ug',     [nt,nlev], 0)
        vg   = _get_or_default(timedep_atm, 'vg',     [nt,nlev], 0)
        wls  = _get_or_default(timedep_atm, 'wls',    [nt,nlev], 0)
        dxq  = _get_or_default(timedep_atm, 'dqtdx',  [nt,nlev], 0)
        dyq  = _get_or_default(timedep_atm, 'dqtdy',  [nt,nlev], 0)
        dtq  = _get_or_default(timedep_atm, 'dqtdt',  [nt,nlev], 0)
        dtth = _get_or_default(timedep_atm, 'dthldt', [nt,nlev], 0)
        dtu  = _get_or_default(timedep_atm, 'dudt',   [nt,nlev], 0)
        dtv  = _get_or_default(timedep_atm, 'dvdt',   [nt,nlev], 0)

        # Write atmospheric data
        for t in range(nt):
            f.write('\n')
            # Write header:
            f.write('{0:^19s} {1:^19s} {2:^19s} {3:^19s} {4:^19s} {5:^19s} {6:^19s} {7:^19s} {8:^19s} {9:^19s}\n'\
                .format('z (m)', 'u_g (m s-1)', 'v_g (m s-1)', 'w_ls (m s-1)',
                        'dqtdx (kg kg-1 m-1)', 'dqtdy (kg kg m-1)', 'dqtdt (kg kg-1 s-1)', 'dthldt (K s-1)', 'dudt (m s-2)', 'dvdt (m s-2)'))

            # Write current time:
            f.write('# {0:1.8E}\n'.format(time[t]))

            # Write profiles:
            for k in range(nlev):
                f.write('{0:+1.12E} {1:+1.12E} {2:+1.12E} {3:+1.12E} {4:+1.12E} {5:+1.12E} {6:+1.12E} {7:+1.12E} {8:+1.12E} {9:+1.12E}\n'\
                    .format(z[k], ug[t,k], vg[t,k], wls[t,k], dxq[t,k], dyq[t,k], dtq[t,k], dtth[t,k], dtu[t,k], dtv[t,k]))

    f.close()


def create_backrad(p, T, q, expnr=1):
    """
    Create the background profiles for RRTMG
    """

    print(' - Saving backrad.inp.{0:03d}.nc'.format(expnr))

    nc_file = nc4.Dataset('backrad.inp.{0:03d}.nc'.format(expnr), 'w')
    dims = nc_file.createDimension('lev', p.size)

    p_var = nc_file.createVariable('lev', 'f4', ('lev'))
    T_var = nc_file.createVariable('T',   'f4', ('lev'))
    q_var = nc_file.createVariable('q',   'f4', ('lev'))

    p_var[:] = p
    T_var[:] = T
    q_var[:] = q

    nc_file.close()
