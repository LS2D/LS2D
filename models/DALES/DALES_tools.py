from collections import OrderedDict as odict

import matplotlib.pyplot as pl
import netCDF4 as nc4
import numpy as np

import datetime
import shutil
import re
import os

# ---------------------------
# "Private" help functions
# ---------------------------
def _get_or_default(dict, name, shape, default_value):
    if name in dict:
        return dict[name]
    else:
        #print('No input profile for \"{}\", defaulting values at zero'.format(name))
        return default_value * np.ones(shape)


def _int_or_float_or_str(value):
    """ Helper function: convert a string to int/float/str """
    try:
        if 'true' in value:
            return True
        elif 'false' in value:
            return False
        elif ('.' in value):
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


def _find_next_dividable_number(number, n):
    if number % n == 0:
        return number
    else:
        return ((number // n)+1)*n


def interp_z(z_input, z_output, variable):
    """
    Interpolate (linear) `variable` from input grid (`z_input`) to output grid (`z_output`)
    """
    return np.interp(z_output, z_input, variable)


def interp_z_time(z_input, z_output, variable):
    """
    Interpolate time varying `variable` from input grid (`z_input`) to output grid (`z_output`)
    """
    data = np.zeros((variable.shape[0], z_output.size))
    for t in range(variable.shape[0]):
        data[t,:] = interp_z(z_input[t,:], z_output, variable[t,:])
    return data


# ---------------------------
# DALES constants (from modglobal.f90)
# ---------------------------
constants = dict(rd = 287.04,
                 cp = 1004.,
                 lv = 2.53e6,
                 p0 = 1e5)

# ---------------------------
# Vertical grids
# ---------------------------
class Grid:
    def __init__(self, kmax, dz0):
        self.kmax = kmax
        self.dz0  = dz0

        self.z = np.zeros(kmax)
        self.dz = np.zeros(kmax)
        self.zsize = None

    def plot(self):
        pl.figure()
        pl.title('zsize = {0:.1f} m'.format(self.zsize), loc='left')
        pl.plot(self.dz, self.z, '-x')
        pl.xlabel('dz (m)')
        pl.ylabel('z (m)')


class Grid_equidist(Grid):
    def __init__(self, kmax, dz0):
        Grid.__init__(self, kmax, dz0)

        self.zsize = kmax * dz0
        self.z[:]  = np.arange(dz0/2, self.zsize, dz0)
        self.dz[:] = dz0


class Grid_stretched(Grid):
    def __init__(self, kmax, dz0, nloc1, nbuf1, dz1):
        Grid.__init__(self, kmax, dz0)

        dn         = 1./kmax
        n          = np.linspace(dn, 1.-dn, kmax)
        nloc1     *= dn
        nbuf1     *= dn
        dzdn1      = dz0/dn
        dzdn2      = dz1/dn

        dzdn       = dzdn1 + 0.5*(dzdn2-dzdn1)*(1. + np.tanh((n-nloc1)/nbuf1))
        self.dz[:] = dzdn*dn

        stretch    = np.zeros(self.dz.size)

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

# ---------------------------
# Function to read DALES in/output
# ---------------------------
class Read_namelist:
    def __init__(self, namelist_file):
        self.namelist_file = namelist_file

        self.groups = {}   # Dictionary holding all the data
        with open(namelist_file) as f:
            for line in f:
                lstrip = line.strip()
                if (len(lstrip) > 0 and lstrip[0] != "#"):
                    if lstrip[0] == '&':
                        curr_group_name = lstrip[1:].lower()
                        self.groups[curr_group_name] = {}
                    elif ("=" in line):
                        var_name = lstrip.split('=')[0].strip()
                        value = _convert_value(lstrip.split('=')[1])
                        self.groups[curr_group_name][var_name] = value

    def __getitem__(self, name):
        if name in self.groups.keys():
            return self.groups[name]
        else:
            raise RuntimeError('Can\'t find group \"{}\" in .ini file'.format(name))

    def __repr__(self):
        return 'Available groups in {}:\n{}'.format(self.namelist_file, ', '.join(self.groups.keys()))


def replace_namelist_value(namelist, variable, new_value, group=None):
    # Backup namelist
    shutil.copy(namelist, '{}.bak'.format(namelist))

    # Read the entire namelist to memory
    with open(namelist, 'r') as source:
        lines = source.readlines()

    # Loop over lines, and replace matches
    curr_group = None
    with open(namelist, 'w') as source:
        for line in lines:
            lstrip = line.strip()
            if len(lstrip) > 0 and lstrip[0] == '&':
                curr_group = lstrip[1:].lower()

            if group is None or curr_group == group.lower():
                source.write(re.sub(r'({}[ |=]).*'.format(variable), r'\1 = {}'.format(new_value), line))
            else:
                source.write(line)

# ---------------------------
# Function to write the DALES input
# ---------------------------
def write_profiles(file_name, variables, nlev, docstring=''):
    """
    Write the prof.inp.xxx input profiles for DALES
    """

    print('Saving {}'.format(file_name))

    f = open(file_name, 'w')

    # Write header (description file)
    if docstring is '':
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

    print('Saving {}'.format(file_name))

    f = open(file_name, 'w')

    # Write header (description file)
    if docstring is '':
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
    if docstring is '':
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

    print('Saving {}'.format(file_name))


    f = open(file_name, 'w')

    # Always write something; DALES expects three line header
    if docstring is '':
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


# ---------------------------
# Function to convert/write the Harmonie LES forcings
# ---------------------------
def get_file_list(path, starttime, endtime):
    """
    Get list of required DDH NetCDF files to force
    LES from `starttime` to `endtime`
    """

    # For now limited to runs starting at a complete hour, to prevent having
    # to interpolate the inital conditions
    if starttime.minute != 0 or starttime.second != 0:
        raise RuntimeError('Can only create forcings starting at a complete hour!')

    # If experiment starts at start of cycle (t=0,3,6,..) we also need the previous cycle..
    if starttime.hour % 3 == 0:
        starttime -= datetime.timedelta(hours=3)

    # Number of cycles to convert
    n_cycles = int((endtime-starttime).total_seconds() / 3600. / 3.) + 1

    # Create list of cycles to convert
    files   = []
    success = True
    for t in range(n_cycles):
        date    = starttime + t * datetime.timedelta(hours=3)
        in_file = '{0:}/{1:04d}/{2:02d}/{3:02d}/{4:02d}/LES_forcing_{1:04d}{2:02d}{3:02d}{4:02d}.nc'.\
            format(path, date.year, date.month, date.day, date.hour)

        files.append(in_file)

        # Check if file actually exists..
        if not os.path.exists(in_file):
            print('ERROR: Can not find input file {}'.format(in_file))
            success = False

    if not success:
        raise RuntimeError('One or more required files could not be found...')
    else:
        return files


def get_start_end_indices(start, end, time):
    """
    Get indices in `time` that correspond to the requested `start` and `end` times
    """

    t0 = np.abs(np.datetime64(start) - time).argmin()
    t1 = np.abs(np.datetime64(end)   - time).argmin() + 1

    print('Using {} (index {}) to {} (index {})'.format(time[t0], t0, time[t1-1], t1-1))

    return t0, t1


def create_initial_profiles(nc_data, grid, t0, t1, iloc, docstring, expnr=1):
    """
    Interpolate Harmonie data onto LES grid,
    and create the `prof.inp` and `scalar.inp` files
    """

    p   = interp_z( nc_data['z'][t0, iloc, :], grid.z, nc_data['p' ][t0, iloc, :] )
    T   = interp_z( nc_data['z'][t0, iloc, :], grid.z, nc_data['T' ][t0, iloc, :] )
    qt  = interp_z( nc_data['z'][t0, iloc, :], grid.z, nc_data['q' ][t0, iloc, :] )
    ql  = interp_z( nc_data['z'][t0, iloc, :], grid.z, nc_data['ql'][t0, iloc, :] )
    u   = interp_z( nc_data['z'][t0, iloc, :], grid.z, nc_data['u' ][t0, iloc, :] )
    v   = interp_z( nc_data['z'][t0, iloc, :], grid.z, nc_data['v' ][t0, iloc, :] )
    tke = np.ones(grid.kmax) * 0.1

    # Conversions from Harmonie -> DALES
    exner  = (p / constants['p0'])**(constants['rd']/constants['cp'])
    theta  = T / exner
    thetal = theta - constants['lv'] / (constants['cp'] * exner) * ql

    # Write to prof.inp.001
    output = odict([('z (m)',        grid.z),
                    ('thl (K)',      thetal),
                    ('qt (kg kg-1)', qt),
                    ('u (m s-1)',    u),
                    ('v (m s-1)',    v),
                    ('tke (m2 s-2)', tke)])

    write_profiles('prof.inp.{0:03d}'.format(expnr), output, grid.kmax, docstring)

    # Initial scalar profiles (for microphysics) are zero
    zero = np.zeros(grid.kmax)
    output = odict([('z (m)',        grid.z),
                    ('qr (kg kg-1)', zero),
                    ('nr (kg kg-1)', zero)])

    write_profiles('scalar.inp.{0:03d}'.format(expnr), output, grid.kmax, docstring)


def create_ls_forcings(nc_data, grid, t0, t1, iloc, docstring, n_accumulate=1, expnr=1, harmonie_rad=False):
    """
    Create all the (partially time dependent) large scale forcings
    """

    # Conversion of Harmonie prognostic variables (T) to LES (thl)
    exner       = (nc_data['p'][:,iloc,:] / constants['p0'])**(constants['rd']/constants['cp'])
    dtthl_dyn_f = nc_data['dtT_dyn'][:,iloc,:] / exner - constants['lv'] / (constants['cp'] * exner) * nc_data['dtql_dyn'][:,iloc,:]
    dtthl_rad_f = nc_data['dtT_rad'][:,iloc,:] / exner - constants['lv'] / (constants['cp'] * exner) * nc_data['dtql_dyn'][:,iloc,:]

    if (n_accumulate > 1):
        # Aaargg.. exceptions, exceptions...
        if (n_accumulate == 2):
            pad = 0     # "Valid" time of accumulated tendencies falls exactly at t==0
        else:
            pad = 1     # First valid time of "    " is after t=0; add padding in front to fix t==0 later...

        # Pick end time (t1) such that t1-t0 is dividable by the number of accumulation steps
        n  = _find_next_dividable_number(t1-t0, n_accumulate)
        nt = int(n / n_accumulate)
        nz = nc_data.dims['level']

        t1 = t0 + n

        z         = np.zeros((nt+pad, nz))
        dtthl_dyn = np.zeros((nt+pad, nz))
        dtthl_rad = np.zeros((nt+pad, nz))
        dtu_dyn   = np.zeros((nt+pad, nz))
        dtv_dyn   = np.zeros((nt+pad, nz))
        dtq_dyn   = np.zeros((nt+pad, nz))

        # Average `n_accumulate` time steps
        z        [pad:,:] = nc_data['z'       ][t0:t1, iloc, :].values.reshape((-1, n_accumulate, nc_data.dims['level'])).mean(axis=1)
        dtthl_dyn[pad:,:] = dtthl_dyn_f        [t0:t1,       :].values.reshape((-1, n_accumulate, nc_data.dims['level'])).mean(axis=1)
        dtthl_rad[pad:,:] = dtthl_rad_f        [t0:t1,       :].values.reshape((-1, n_accumulate, nc_data.dims['level'])).mean(axis=1)
        dtu_dyn  [pad:,:] = nc_data['dtu_dyn' ][t0:t1, iloc, :].values.reshape((-1, n_accumulate, nc_data.dims['level'])).mean(axis=1)
        dtv_dyn  [pad:,:] = nc_data['dtv_dyn' ][t0:t1, iloc, :].values.reshape((-1, n_accumulate, nc_data.dims['level'])).mean(axis=1)
        dtq_dyn  [pad:,:] = nc_data['dtqv_dyn'][t0:t1, iloc, :].values.reshape((-1, n_accumulate, nc_data.dims['level'])).mean(axis=1)
    else:
        z         = nc_data['z'       ][t0:t1, iloc, :].values
        dtthl_dyn = dtthl_dyn_f        [t0:t1,       :].values
        dtthl_rad = dtthl_rad_f        [t0:t1,       :].values
        dtu_dyn   = nc_data['dtu_dyn' ][t0:t1, iloc, :].values
        dtv_dyn   = nc_data['dtv_dyn' ][t0:t1, iloc, :].values
        dtq_dyn   = nc_data['dtqv_dyn'][t0:t1, iloc, :].values

    # Time in seconds since start of run
    time_sec = (nc_data.time[t0:t1  ] - nc_data.time[t0]).values / 1e9

    # Forcings are accumulated, so "valid" time is half a dt earlier..
    dt = time_sec[1] - time_sec[0]      # This is a bit inaccurate....
    time_sec_ls = time_sec - 0.5 * dt

    # Calculate valid time of accumulated tendencies
    if (n_accumulate > 1):
        if (pad == 0):
            time_sec_ls = time_sec_ls.reshape((-1, n_accumulate)).mean(axis=1)
        else:
            tmp = time_sec_ls.copy()
            time_sec_ls = np.zeros(nt+1)
            time_sec_ls[1:] = tmp.reshape((-1, n_accumulate)).mean(axis=1)
    else:
        time_sec_ls[0] = 0

    # Fix t==0 forcings (if necessary)
    if (n_accumulate == 1 or n_accumulate > 2):
        dtthl_dyn[0,:] = 0.5 * (dtthl_dyn_f        [t0-1,     :] + dtthl_dyn_f        [t0,     :])
        dtthl_rad[0,:] = 0.5 * (dtthl_rad_f        [t0-1,     :] + dtthl_rad_f        [t0,     :])
        dtu_dyn  [0,:] = 0.5 * (nc_data['dtu_dyn' ][t0-1,iloc,:] + nc_data['dtu_dyn' ][t0,iloc,:])
        dtv_dyn  [0,:] = 0.5 * (nc_data['dtv_dyn' ][t0-1,iloc,:] + nc_data['dtv_dyn' ][t0,iloc,:])
        dtq_dyn  [0,:] = 0.5 * (nc_data['dtqv_dyn'][t0-1,iloc,:] + nc_data['dtqv_dyn'][t0,iloc,:])

    # Interpolate onto LES grid
    dtthl_dyn = interp_z_time(z, grid.z, dtthl_dyn)
    dtthl_rad = interp_z_time(z, grid.z, dtthl_rad)
    dtu_dyn   = interp_z_time(z, grid.z, dtu_dyn  )
    dtv_dyn   = interp_z_time(z, grid.z, dtv_dyn  )
    dtq_dyn   = interp_z_time(z, grid.z, dtq_dyn  )
    zero_a    = np.zeros_like(dtthl_dyn)

    if (harmonie_rad):
        print('Adding radiative tendency from Harmonie..')
        dtthl_dyn += dtthl_rad

    # Surface forcings
    time_sec_sfc = time_sec[::n_accumulate]
    ps = nc_data['p_s'][t0:t1:n_accumulate, iloc].values
    Ts = nc_data['T_s'][t0:t1:n_accumulate, iloc].values
    qs = nc_data['q_s'][t0:t1:n_accumulate, iloc].values
    zero_s = np.zeros_like(Ts)

    # Write to ls_flux.inp.expnr
    output_sfc = odict([('time',   time_sec_sfc),
                        ('p_s',    ps          ),
                        ('T_s',    Ts          ),
                        ('qt_s',   qs          )])

    output_ls  = odict([('time',   time_sec_ls),
                        ('z',      grid.z     ),
                        ('dqtdt',  dtq_dyn    ),
                        ('dthldt', dtthl_dyn  ),
                        ('dudt',   dtu_dyn    ),
                        ('dvdt',   dtv_dyn    )])

    write_forcings('ls_flux.inp.{0:03d}'.format(expnr), output_sfc, output_ls, docstring)


    # Dummy forcings for the microphysics scalars
    write_dummy_forcings('ls_fluxsv.inp.{0:03d}'.format(expnr), 2, grid.z, docstring)

    # Also create non-time dependent file (lscale.inp), required by DALES (why?)
    zero = np.zeros_like(grid.z)
    output_ls2  = odict([('height', grid.z), ('ug', zero), ('vg', zero), ('wfls', zero), \
                         ('dqtdxls', zero), ('dqtdyls', zero), ('dqtdtls', zero), ('dthldt', zero)])
    write_profiles('lscale.inp.{0:03d}'.format(expnr), output_ls2, grid.kmax, docstring)

    # Return data from debugging....
    return output_sfc, output_ls


def create_nudging_profiles(nc_data, grid, nudgefactor, t0, t1, iloc, docstring, interval=1, expnr=1):
    """
    Create the nudging profiles
    """

    # Time in seconds since start of run
    time_sec    = (nc_data.time[t0:t1  ] - nc_data.time[t0]).values / 1e9

    # Vertical profiles
    T  = interp_z_time( nc_data['z'][t0:t1, iloc, :], grid.z, nc_data['T' ][t0:t1, iloc, :] )
    u  = interp_z_time( nc_data['z'][t0:t1, iloc, :], grid.z, nc_data['u' ][t0:t1, iloc, :] )
    v  = interp_z_time( nc_data['z'][t0:t1, iloc, :], grid.z, nc_data['v' ][t0:t1, iloc, :] )
    p  = interp_z_time( nc_data['z'][t0:t1, iloc, :], grid.z, nc_data['p' ][t0:t1, iloc, :] )
    qt = interp_z_time( nc_data['z'][t0:t1, iloc, :], grid.z, nc_data['q' ][t0:t1, iloc, :] )
    ql = interp_z_time( nc_data['z'][t0:t1, iloc, :], grid.z, nc_data['ql'][t0:t1, iloc, :] )
    zero = np.zeros_like(T)

    # Nudging factor (0-1) with height; is multiplied with nudging time from namelist
    nudgefac = np.zeros_like(T)
    nudgefac[:,:] = nudgefactor

    # Conversions from Harmonie -> DALES
    exner  = (p / constants['p0'])**(constants['rd']/constants['cp'])
    theta  = T / exner
    thetal = theta - constants['lv'] / (constants['cp'] * exner) * ql

    # Write to nudge.inp
    output = odict([('z (m)',        grid.z),
                    ('factor (-)',   nudgefac[::interval,:]),
                    ('u (m s-1)',    u       [::interval,:]),
                    ('v (m s-1)',    v       [::interval,:]),
                    ('w (m s-1)',    zero    [::interval,:]),
                    ('thl (K)',      thetal  [::interval,:]),
                    ('qt (kg kg-1)', qt      [::interval,:])])

    write_time_profiles('nudge.inp.{0:03d}'.format(expnr), time_sec[::interval], output, grid.kmax, docstring)


def create_backrad(nc_data, t0, iloc, expnr=1):
    """
    Create the background profiles for RRTMG
    """

    print('Saving backrad.inp.{0:03d}.nc'.format(expnr))

    nc_file = nc4.Dataset('backrad.inp.{0:03d}.nc'.format(expnr), 'w')
    dims = nc_file.createDimension('lev', nc_data.dims['level'])

    p = nc_file.createVariable('lev', 'f4', ('lev'))
    T = nc_file.createVariable('T',   'f4', ('lev'))
    q = nc_file.createVariable('q',   'f4', ('lev'))

    p[:] = nc_data['p'][t0, iloc, :]
    T[:] = nc_data['T'][t0, iloc, :]
    q[:] = nc_data['q'][t0, iloc, :]

    nc_file.close()






if __name__ == '__main__':
    #
    # Only executed (for debugging/examples) if script is called directly
    #

    pl.close('all')
    pl.ion()

    if True:
        """
        Demo of the different vertical grids
        """
        ktot = 192
        dz0  = 20

        equidist  = Grid_equidist(ktot, dz0)
        linear    = Grid_linear_stretched(ktot, dz0, 0.01)
        stretched = Grid_stretched(ktot, dz0, 110, 30, 150)

        pl.figure()
        pl.plot(equidist.dz, equidist.z, '-x', label='equidistant')
        pl.plot(linear.dz, linear.z, '-x', label='linear')
        pl.plot(stretched.dz, stretched.z, '-x', label='stretched')
        pl.legend()
        pl.xlabel('dz (m)')
        pl.ylabel('z (m)')

    if False:
        """
        Read DALES namelist
        """
        nl = Read_namelist('namoptions.001')

        print(nl)
        print('variables in domain:', nl['domain'])
        print('itot:', nl['domain']['itot'])

    if False:
        """
        Write initial vertical profiles
        Profiles which aren't specified are initialized at zero.
        """
        nlev  = 64
        zsize = 3200
        dz    = zsize/nlev
        z     = np.arange(0.5*dz, zsize, dz)

        th  = 290 + z * 0.006
        qt  = np.zeros(nlev)
        u   = np.ones(nlev) * 1
        v   = np.ones(nlev) * 2
        tke = np.ones(nlev) * 2

        data = odict([('z', z), ('thl', th), ('qt', qt), ('u', u), ('v', v), ('tke',tke)])
        write_profiles('prof.inp.001', data, docstring='Example of initial DALES profiles')

    if False:
        """
        Write the time dependent surface variables,
        and large scale forcings.
        """
        nlev  = 64
        zsize = 3200
        dz    = zsize/nlev
        z     = np.arange(0.5*dz, zsize, dz)

        time     = np.arange(0,7200.01,300)
        thls     = np.ones(time.size)*300
        data_sfc = odict([('time', time), ('thl_s', thls)])

        time_ls  = np.arange(0,7200.01, 1800)
        ug       = np.ones((time_ls.size, nlev))*5
        vg       = np.ones((time_ls.size, nlev))*-5
        data_ls  = odict([('time', time_ls), ('z', z), ('u_g', ug), ('v_g', vg)])

        write_forcings('ls_flux.inp.001', data_sfc, data_ls, 'Example of DALES forcings')

    if True:
        """
        Write the time dependent profiles (e.g. nudging),
        """

        nlev  = 32
        zsize = 3200
        dz    = zsize/nlev
        z     = np.arange(0.5*dz, zsize, dz)
        time  = np.arange(0,7200.01,1800)

        var1 = np.arange(time.size*nlev).reshape(time.size,-1)
        var2 = np.arange(time.size*nlev).reshape(time.size,-1)+1
        var3 = np.arange(time.size*nlev).reshape(time.size,-1)+2

        data = odict([('z (m)', z), ('var1', var1), ('var2', var2), ('var3', var3)])

        write_time_profiles('nudge.inp.001', time, data, nlev) 



