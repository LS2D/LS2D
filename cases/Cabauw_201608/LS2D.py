import matplotlib.pyplot as pl
import datetime
import numpy as np
import sys
import os
import shutil

# Add `src` subdirectory of LS2D to Python path
abs_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append('{}/../../src/'.format(abs_path))
sys.path.append('{}/..'.format(abs_path))

# Import the LS2D/MicroHH specific scripts
from download_ERA5 import download_ERA5
from read_ERA5     import Read_ERA
from messages      import header, message, error

import microhh_tools as mht
from grid import Grid_stretched

# Path to the ERA5 files
era5_path = '/archive/bstratum/ERA5/'

# Working directory; individual cases are placed in `YYYYMMDD_tHH` subdirectory
workdir = '/home/bstratum/scratch/cabauw_aug2016/'

# Path to MicroHH binary
microhh_bin = '/home/bstratum/models/microhh/build_dp_cpumpi/microhh'

# Path to RRTMGP repository, for radiation coefficient files.
rrtmgp_path = '/home/bstratum/models/rte-rrtmgp/'

float_type  = 'f8'    # MicroHH float type ('f4', 'f8')
link_files = False    # Switch between linking or copying files

# Time of day to simulate
start_hour = 5
end_hour = 19

# Days in Aug 2016:
start_day = 4
end_day = 5#18

for day in range(start_day, end_day):

    start = datetime.datetime(year=2016, month=8, day=day, hour=start_hour)
    end   = datetime.datetime(year=2016, month=8, day=day, hour=end_hour)

    # Dictionary with settings
    settings = {
        'central_lat' : 51.971,
        'central_lon' : 4.927,
        'area_size'   : 1,
        'case_name'   : 'cabauw',
        'base_path'   : era5_path,      # Arch
        'start_date'  : start,
        'end_date'    : end,
        'write_log'   : False,
        'data_source' : 'MARS',
        'ntasks'      : 1
        }

    header('Creating LES input')

    #
    # Download the ERA5 data (or check whether it is available local).
    #
    download_ERA5(settings)

    #
    # Read ERA5 data, and calculate LES forcings, using +/-`n_av` grid point averages in ERA5.
    #
    e5 = Read_ERA(settings)
    e5.calculate_forcings(n_av=0, method='4th')

    #
    # Read MicroHH namelist and create stretched vertical grid
    #
    grid = Grid_stretched(kmax=224, dz0=20, nloc1=100, nbuf1=20, dz1=100, nloc2=210, nbuf2=10, dz2=500)
    #grid.plot()

    #
    # Create nudge factor, controlling where nudging is aplied, and time scale
    #
    tau_nudge = 10800        # Nudge time scale (s)
    nudge_fac = np.ones(grid.z.size) / tau_nudge

    #
    # Interpolate ERA5 onto LES grid
    #
    def interp_time(z, ze, arr):
        out = np.empty((arr.shape[0], z.size))
        for i in range(arr.shape[0]):
            out[i,:] = np.interp(z, ze[i,:], arr[i,:])
        return out

    thl   = interp_time(grid.z, e5.z_mean, e5.thl_mean   )
    qt    = interp_time(grid.z, e5.z_mean, e5.qt_mean    )
    u     = interp_time(grid.z, e5.z_mean, e5.u_mean     )
    v     = interp_time(grid.z, e5.z_mean, e5.v_mean     )
    w     = interp_time(grid.z, e5.z_mean, e5.wls_mean   )
    p     = interp_time(grid.z, e5.z_mean, e5.p_mean     )
    thlls = interp_time(grid.z, e5.z_mean, e5.dtthl_advec)
    qtls  = interp_time(grid.z, e5.z_mean, e5.dtqt_advec )
    uls   = interp_time(grid.z, e5.z_mean, e5.dtu_advec  )
    vls   = interp_time(grid.z, e5.z_mean, e5.dtv_advec  )
    ug    = interp_time(grid.z, e5.z_mean, e5.ug         )
    vg    = interp_time(grid.z, e5.z_mean, e5.vg         )
    o3    = interp_time(grid.z, e5.z_mean, e5.o3_mean    )

    # Surface / soil
    z_soil = np.array([-1.945, -0.64, -0.175, -0.035])
    soil_index = np.ones_like(z_soil)*2

    #
    # Radiation profiles for RRTMGP
    #
    co2 = 348.e-6
    ch4 = 1650.e-9
    n2o = 306.e-9
    n2  = 0.7808
    o2  = 0.2095

    g1  = 3.6478
    g2  = 0.83209
    g3  = 11.3515

    # Background profiles on pressure levels
    z_lay = e5.z_mean [0,:]
    z_lev = e5.zh_mean[0,:]

    p_lay = e5.p_mean [0,:]
    p_lev = e5.ph_mean[0,:]

    T_lay = e5.T_mean [0,:]
    T_lev = e5.Th_mean[0,:]

    h2o_rad = e5.qt_mean[0,:]
    co2_rad = np.ones(e5.nfull) * co2
    ch4_rad = np.ones(e5.nfull) * ch4
    n2o_rad = np.ones(e5.nfull) * n2o
    n2_rad  = np.ones(e5.nfull) * n2
    o2_rad  = np.ones(e5.nfull) * o2
    o3_rad  = e5.o3_mean[0,:]

    # Profiles on LES grid
    h2o_atmo = qt[0,:]
    co2_atmo = np.ones(grid.kmax) * co2
    ch4_atmo = np.ones(grid.kmax) * ch4
    n2o_atmo = np.ones(grid.kmax) * n2o
    n2_atmo  = np.ones(grid.kmax) * n2
    o2_atmo  = np.ones(grid.kmax) * o2
    o3_atmo  = o3[0,:]

    #
    # Write MicroHH input
    #
    message('Writing forcings as LES input')

    #
    # Update namelist variables
    #
    nl_file = '{}.ini'.format(settings['case_name'])
    nl_backup = '{}.ini.bak'.format(settings['case_name'])

    shutil.copyfile(nl_file, nl_backup)
    nl = mht.read_namelist(nl_file)

    nl['grid']['ktot'] = grid.kmax
    nl['grid']['zsize'] = grid.zsize
    nl['time']['endtime'] = e5.time_sec.max()
    nl['force']['fc'] = e5.fc
    nl['radiation']['lon'] = settings['central_lon']
    nl['radiation']['lat'] = settings['central_lat']

    datetime_utc = '{0:04d}-{1:02d}-{2:02d} {3:02d}:{4:02d}:{5:02d}'.format(
            start.year, start.month, start.day, start.hour, start.minute, start.second)
    nl['time']['datetime_utc'] = datetime_utc

    mht.write_namelist(nl_file, nl)

    #
    # Write NetCDF file
    #
    init_profiles = {
            'z': grid.z, 'thl': thl[0,:], 'qt': qt[0,:], 'u': u[0,:],
            'v': v[0,:], 'nudgefac': nudge_fac, 'co2': co2_atmo, 'ch4': ch4_atmo,
            'n2o': n2o_atmo, 'n2': n2_atmo, 'o2': o2_atmo, 'o3': o3_atmo, 'h2o': h2o_atmo}

    radiation  = {
            'z_lay': z_lay, 'z_lev': z_lev, 'p_lay': p_lay, 'p_lev': p_lev,
            't_lay': T_lay, 't_lev': T_lev, 'co2': co2_rad, 'ch4': ch4_rad,
            'n2o': n2o_rad, 'n2': n2_rad, 'o2': o2_rad, 'o3': o3_rad, 'h2o': h2o_rad}

    tdep_surface = {
            'time_surface': e5.time_sec, 'thl_sbot': e5.wths_mean,
            'qt_sbot': e5.wqs_mean, 'p_sbot': e5.ps_mean }

    tdep_ls = {
            'time_ls': e5.time_sec, 'u_geo': ug, 'v_geo': vg, 'w_ls': w,
            'thl_ls': thlls, 'qt_ls': qtls, 'u_ls': uls, 'v_ls': vls,
            'thl_nudge': thl, 'qt_nudge': qt, 'u_nudge': u, 'v_nudge': v}

    soil = {'z': z_soil, 'theta': e5.phisoil_mean[0,::-1], 't': e5.Tsoil_mean[0,::-1], 'index': soil_index}

    mht.write_NetCDF_input('cabauw', float_type, init_profiles, tdep_surface, tdep_ls, radiation, soil)

    #
    # Copy/move/link files to working directory.
    #
    path = '{0}/{1:04d}{2:02d}{3:02d}_t{4:02d}'.format(workdir, start.year, start.month, start.day, start.hour)
    if os.path.exists(path):
        error('Work directory {} already exists!!'.format(path))
    else:
        os.makedirs(path)

    to_copy = ['cabauw.ini', '../van_genuchten_parameters.nc']
    to_move = ['cabauw_input.nc']
    to_link = {
            'microhh': microhh_bin,
            'coefficients_lw.nc':
                '{}/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc'.format(rrtmgp_path),
            'coefficients_sw.nc':
                '{}/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc'.format(rrtmgp_path),
            'cloud_coefficients_lw.nc':
                '{}/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc'.format(rrtmgp_path),
            'cloud_coefficients_sw.nc':
            '{}/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc'.format(rrtmgp_path)}

    for f in to_copy:
        shutil.copy(f, path)
    for f in to_move:
        shutil.move(f, path)

    for dst,src in to_link.items():
        if link_files:
            os.symlink(src, '{}/{}'.format(path,dst))
        else:
            shutil.copy(src, '{}/{}'.format(path,dst))

    # Restore namelist file
    shutil.copyfile(nl_backup, nl_file)
