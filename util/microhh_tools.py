import netCDF4 as nc4
import numpy   as np
import struct  as st
import glob
import re
from collections import OrderedDict

# -------------------------
# General help functions
# -------------------------

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


def _find_namelist_file():
    """ Helper function: automatically find the .ini file in the current directory """
    namelist_file = glob.glob('*.ini')
    if len(namelist_file) == 0:
        raise RuntimeError('Can\'t find any .ini files in the current directory!')
    if len(namelist_file) > 1:
        raise RuntimeError('There are multiple .ini files: {}'.format(namelist_file))
    else:
        return namelist_file[0]


def _process_endian(endian):
    if endian not in ['little', 'big']:
        raise ValueError('endian has to be \"little\" or \"big\"!')
    endian = '<' if endian == 'little' else '>'
    return endian


# -------------------------
# Classes and functions to read and write MicroHH things
# -------------------------
def read_namelist(namelist_file=None):
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
                f.write('{}={}\n'.format(variable, value))
            f.write('\n')


class Read_grid:
    """ Read the grid file from MicroHH.
        If no file name is provided, grid.0000000 from the current directory is read """
    def __init__(self, itot, jtot, ktot, zsize, filename=None, endian='little'):
        self.en  = _process_endian(endian)
        filename = 'grid.0000000' if filename is None else filename

        self.fin = open(filename, 'rb')
        self.x  = self.read(itot)
        self.xh = self.read(itot)
        self.y  = self.read(jtot)
        self.yh = self.read(jtot)
        self.z  = self.read(ktot)
        self.zh = self.read(ktot)
        self.fin.close()
        del self.fin

        self.itot = self.x.size
        self.jtot = self.y.size
        self.ktot = self.z.size

        self.dx = self.x[1]-self.x[0] if itot > 1 else self.xh[1]
        self.dy = self.y[1]-self.y[0] if jtot > 1 else self.yh[1]

        self.xsize = self.itot * self.dx
        self.ysize = self.jtot * self.dy
        self.zsize = zsize

    def read(self, n):
        return np.array(st.unpack('{0}{1}d'.format(self.en, n), self.fin.read(n*8)))


def read_restart_file(path, itot, jtot, ktot, endian='little'):
    """ Read a MicroHH restart file into a 3D (or 2D if ktot=1) numpy array
        The returned array has the dimensions ordered as [z,y,x] """

    en = _process_endian(endian)

    f  = open(path, 'rb')
    if (ktot > 1):
        field = np.zeros((ktot, jtot, itot))
        for k in range(ktot):
            raw = f.read(itot*jtot*8)
            tmp = np.array(st.unpack('{0}{1}d'.format(en, itot*jtot), raw))
            field[k,:,:] = tmp.reshape((jtot, itot))[:,:]
        f.close()
    else:
        raw = f.read(itot*jtot*8)
        tmp = np.array(st.unpack('{0}{1}d'.format(en, itot*jtot), raw))
        field = tmp.reshape((jtot, itot))

    return field


def write_restart_file(data, itot, jtot, ktot, path, per_slice=True, endian='little'):
    """ Write a restart file in the format requires by MicroHH.
        The input array should be indexed as [z,y,x] """

    en = _process_endian(endian)

    if(per_slice):
        # Write level by level (less memory hungry.....)
        fout  = open(path, "wb")
        for k in range(ktot):
            tmp  = data[k,:,:].reshape(itot*jtot)
            tmp2 = st.pack('{0}{1}d'.format(en, tmp.size), *tmp)
            fout.write(tmp2)
        fout.close()
    else:
        # Write entire field at once (memory hungry....)
        tmp  = data.reshape(data.size)
        tmp2 = st.pack('{0}{1}d'.format(en, tmp.size), *tmp)
        fout = open(path, "wb")
        fout.write(tmp2)
        fout.close()


def write_NetCDF_input(
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
        dims = None if not isinstance(data, np.ndarray) else 'z'
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
            if not isinstance(data, np.ndarray):
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
