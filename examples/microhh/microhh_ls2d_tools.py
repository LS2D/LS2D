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
from collections import OrderedDict
import struct as st
import re
import sys

# Third party modules
import netCDF4 as nc4
import xarray as xr
import numpy as np


# General help functions
def _int_or_float_or_str(value):
    """ Helper function: convert a string to int/float/str """
    try:
        if ('.' in value):
            return float(value)
        else:
            return int(float(value))
    except:
        return value.rstrip()


def _convert_value(value):
    """ Helper function: convert namelist value or list """
    if ',' in value:
        value = value.split(',')
        return [_int_or_float_or_str(val) for val in value]
    else:
        return _int_or_float_or_str(value)


def read_namelist(namelist_file):
    """
    Read a .ini namelist into a nested dictionary (dict[group][variable])
    """

    ini = OrderedDict()
    with open(namelist_file, 'r') as f:
        for line in f:
            lstrip = line.strip()
            if (len(lstrip) > 0 and lstrip[0] != "#"):
                if lstrip[0] == '[' and lstrip[-1] == ']':
                    group_name = lstrip[1:-1]
                    ini[group_name] = OrderedDict()
                elif ("=" in line):
                    var_name = lstrip.split('=')[0]
                    value = _convert_value(lstrip.split('=')[1])
                    ini[group_name][var_name] = value
    return ini


def write_namelist(namelist_file, namelist_dict):
    """
    Write a .ini namelist from a (nested) dictionary
    """

    with open(namelist_file, 'w') as f:
        for group, items in namelist_dict.items():
            f.write('[{}]\n'.format(group))
            for variable,value in items.items():
                if isinstance(value, list):
                    value = ','.join([str(elem) for elem in value])
                elif isinstance(value, bool):
                    value = 1 if value==True else 0
                f.write('{}={}\n'.format(variable, value))
            f.write('\n')


def write_netcdf_input(
        case_name, float_type, init_profiles, tdep_surface=None,
        tdep_ls=None, radiation=None, soil=None):
    """
    Function for writing the MicroHH2 NetCDF input
    """

    def add_variable(nc_group, name, dims, data, float_type):
        """
        Add variable to NetCDF file (or group), and write data
        """
        if dims is None:
            var = nc_group.createVariable(name, float_type)
            var[:] = data
        else:
            var = nc_group.createVariable(name, float_type, dims)
            var[:] = data[:]

    def is_array(data):
        """
        Check if value if array or scalar
        """
        if isinstance(data, np.ndarray) or isinstance(data, xr.DataArray):
            return True
        return False

    # Define new NetCDF file
    nc_file = nc4.Dataset('{}_input.nc'.format(case_name), mode='w', datamodel='NETCDF4')

    # Create height dimension, and set height coordinate
    nc_file.createDimension('z', init_profiles['z'].size)
    add_variable(nc_file, 'z', ('z'), init_profiles['z'], float_type)

    # Create a group called "init" for the initial profiles.
    nc_group_init = nc_file.createGroup('init')

    # Set the initial profiles
    for name, data in init_profiles.items():
        # Switch between vector and scalar values
        dims = 'z' if is_array(data) else None
        add_variable(nc_group_init, name, dims, data, float_type)

    # Create a group called "timedep" for the time dependent input
    if tdep_surface is not None or tdep_ls is not None:
        nc_group_timedep = nc_file.createGroup('timedep')

    # Write the time dependent surface values
    if tdep_surface is not None:
        nc_group_timedep.createDimension('time_surface', tdep_surface['time_surface'].size)

        for name, data in tdep_surface.items():
            add_variable(nc_group_timedep, name, ('time_surface'), data, float_type)

    # Write the time dependent atmospheric values
    if tdep_ls is not None:
        nc_group_timedep.createDimension('time_ls', tdep_ls['time_ls'].size)

        for name, data in tdep_ls.items():
            dims = ('time_ls') if name == 'time_ls' else ('time_ls', 'z')
            add_variable(nc_group_timedep, name, dims, data, float_type)

    if radiation is not None:
        nc_group_rad = nc_file.createGroup('radiation')

        nc_group_rad.createDimension("lay", radiation['p_lay'].size)
        nc_group_rad.createDimension("lev", radiation['p_lev'].size)

        for name, data in radiation.items():
            # Switch between vector and scalar values
            if not is_array(data):
                dims = None
            else:
                dims = ('lay') if data.size == radiation['p_lay'].size else ('lev')

            add_variable(nc_group_rad, name, dims, data, float_type)

    if soil is not None:
        nc_group_soil = nc_file.createGroup('soil')
        nc_group_soil.createDimension("z", soil['z'].size)

        for name, data in soil.items():
            add_variable(nc_group_soil, name, 'z', data, float_type)

    nc_file.close()


def calc_root_frac(z, a_r, b_r):
    """
    Calculate root fraction from `a_r` and `b_r` coeffients.
    See IFS documentation..
    """

    # Calculate half level soil depths
    zh = np.zeros(z.size+1)
    for k in range(zh.size-2, -1, -1):
        zh[k] = zh[k+1] - 2*(zh[k+1] - z[k])

    # Calculate root fraction
    root_frac = np.zeros_like(z)
    for k in range(1, root_frac.size):
        root_frac[k] = 0.5 * (np.exp(a_r * zh[k+1]) + \
                              np.exp(b_r * zh[k+1]) - \
                              np.exp(a_r * zh[k  ]) - \
                              np.exp(b_r * zh[k  ]));

    root_frac[0] = 1.-root_frac.sum()

    return root_frac


def check_grid_decomposition(itot, jtot, ktot, npx, npy):
    """
    Check whether grid / MPI decomposition is valid
    """

    print('Checking grid: itot={}, jtot={}, ktot={}, npx={}, npy={}'.format(
        itot, jtot, ktot, npx, npy))

    err = False
    if itot%npx != 0:
        print('ERROR: itot%npx != 0')
        err = True

    if itot%npy != 0:
        print('ERROR: itot%npy != 0')
        err = True

    if jtot%npx != 0 and npy > 1:
        print('ERROR: jtot%npx != 0')
        err = True

    if jtot%npy != 0:
        print('ERROR: jtot%npy != 0')
        err = True

    if ktot%npx != 0:
        print('ERROR: ktot%npx != 0')
        err = True

    if err:
        sys.exit('Invalid grid configuration!')
    else:
        print('Grid okay!')
