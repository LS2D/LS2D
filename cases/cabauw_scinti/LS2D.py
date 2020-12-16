import matplotlib.pyplot as pl
import xarray as xr
import numpy as np
import sys
import os
import shutil
import subprocess
from datetime import datetime, timedelta
from scipy.interpolate import interp1d

# Add `src` subdirectory of LS2D to Python path
LS2D_root = '/home/bart/meteo/models/LS2D/'
sys.path.append(f'{LS2D_root}/src/')
sys.path.append(f'{LS2D_root}/util/')

# Import the LS2D/MicroHH specific scripts
from download_ERA5 import download_ERA5
from read_ERA5 import Read_ERA
from messages import header, message, error

import microhh_tools as mht
from grid import Grid_stretched, Grid_linear_stretched, Grid_stretched_manual, check_grid_decomposition

pl.close('all'); pl.ion()

# Constants
T0  = 273.15
cp  = 1004.7090
Rd  = 287.0597
Rv  = 461.5250
eps = Rd/Rv

def refine_dz(dz_in, fac):
    new_dz = dz_in / fac
    ktot_out = int(dz_in.size*fac)
    dz_out = np.zeros(ktot_out)
    for k in range(ktot_out):
        dz_out[k] = new_dz[k//fac]
    return dz_out

def interpolate_soil(variable_in, z_source, z_goal, kind='linear'):
    return interp1d(z_source, variable_in, fill_value='extrapolate', kind=kind)(z_goal)

if __name__ == '__main__':

    # Switch between different systems:
    env_cartesius = {
            'system': 'cartesius',
            'era5_path': '/archive/bstratum/ERA5/',
            'work_path': '/home/bstratum/scratch/cabauw_aug2016_small24h/',
            'microhh_bin': '/home/bstratum/models/microhh/build_dp_cpumpi/microhh',
            'rrtmgp_path': '/home/bstratum/models/rte-rrtmgp/',
            'auto_submit': True,
            'set_lfs_stripe': True,
            'link_files': False}

    env_arch = {
            'system': 'arch',
            'era5_path': '/home/scratch1/meteo_data/LS2D/',
            'work_path': '/home/bart/meteo/models/LS2D/cases/cabauw_scinti/',
            'microhh_bin': '/home/bart/meteo/models/microhh/build_dp_cpumpi/microhh',
            'rrtmgp_path': '/home/bart/meteo/models/rte-rrtmgp/',
            'auto_submit': False,
            'set_lfs_stripe': False,
            'link_files': False}

    # Switch between different systems:
    env = env_arch

    float_type  = 'f8'    # MicroHH float type ('f4', 'f8')

    start = datetime(year=2020, month=9, day=14, hour=6,  minute=0)
    end   = datetime(year=2020, month=9, day=14, hour=16, minute=0)

    hetero = True

    # Dictionary with settings
    settings = {
        'central_lat' : 51.971,
        'central_lon' : 4.927,
        'area_size'   : 1,
        'case_name'   : 'cabauw_nrt',
        'base_path'   : env['era5_path'],
        'start_date'  : start,
        'end_date'    : end,
        'write_log'   : False,
        'data_source' : 'MARS',
        'ntasks'      : 1
        }

    header('Creating LES input')

    # Download the ERA5 data (or check whether it is available local).
    download_ERA5(settings)

    # Read ERA5 data, and calculate LES forcings, using +/-`n_av` grid point averages in ERA5.
    e5 = Read_ERA(settings)
    e5.calculate_forcings(n_av=0, method='4th')

    # Read MicroHH namelist and create stretched vertical grid
    #grid = Grid_linear_stretched(kmax=512, dz0=2, alpha=.01)
    #grid = Grid_linear_stretched(kmax=128, dz0=20, alpha=0.009)
    #grid = Grid_linear_stretched(kmax=128, dz0=10, alpha=0.012)
    #grid = Grid_linear_stretched(kmax=160, dz0=4, alpha=0.012)

    #stretch_heights = np.array([0, 200, 2000, 5000, 11000, 50000])
    #stretch_factors = np.array([1.025, 1.011, 1.006, 1.02, 1.08])

    #stretch_heights = np.array([0, 600, 50000])
    #stretch_factors = np.array([1.03, 1.06])
    #grid = Grid_stretched_manual(64, 10., stretch_heights, stretch_factors)

    #stretch_heights = np.array([0, 600, 50000])
    #stretch_factors = np.array([1.005, 1.03])
    #grid = Grid_stretched_manual(192, 4., stretch_heights, stretch_factors)

    #stretch_heights = np.array([0, 600, 500000])
    #stretch_factors = np.array([1.00, 1.05])
    #grid = Grid_stretched_manual(384, 2., stretch_heights, stretch_factors)

    stretch_heights = np.array([0, 600, 500000])
    stretch_factors = np.array([1.00, 1.04])
    grid = Grid_stretched_manual(288, 3., stretch_heights, stretch_factors)

    grid.plot()

    # Interpolate ERA5 variables and forcings onto LES grid
    variables = [
            'thl', 'qt', 'u', 'v', 'wls', 'p',
            'dtthl_advec', 'dtqt_advec', 'dtu_advec', 'dtv_advec',
            'ug' ,'vg' ,'o3', 'z']
    e5_at_z = e5.interpolate_to_fixed_height(variables, grid.z)

    # Create nudge factor, controlling where nudging is aplied, and time scale
    print('NOTE: using reduced tau_nudge....')
    #tau_nudge = 10800        # Nudge time scale (s)
    tau_nudge = 21600        # Nudge time scale (s)
    nudge_fac = np.ones(grid.z.size) / tau_nudge

    # Surface / soil
    # 1. IFS/ERA grid
    dz_soil_0 = np.array([0.07, 0.21, 0.72, 1.89])
    zh_soil_0 = np.cumsum(dz_soil_0)
    zh_soil_0 = -np.insert(zh_soil_0, 0, 0)[::-1]
    z_soil_0  = 0.5*(zh_soil_0[1:] + zh_soil_0[:-1])

    # 2. New grid
    #dz_soil = refine_dz(dz_soil_0, 5)
    #zh_soil = np.cumsum(dz_soil)
    #zh_soil = -np.insert(zh_soil, 0, 0)[::-1]
    #z_soil  = 0.5*(zh_soil[1:] + zh_soil[:-1])
    zh_soil = zh_soil_0
    z_soil = z_soil_0

    soil_index_era = e5.soil_type_nn[0]-1  # -1 = Fortran->C++ indexing
    soil_index = np.ones_like(z_soil)*soil_index_era

    # Calculate root fraction
    # Only used for homogeneous land-surface
    ar = 10.739
    br = 2.608
    root_frac = np.zeros_like(z_soil)
    for k in range(1, root_frac.size):
        root_frac[k] = 0.5 * (np.exp(ar * zh_soil[k+1]) + \
                              np.exp(br * zh_soil[k+1]) - \
                              np.exp(ar * zh_soil[k  ]) - \
                              np.exp(br * zh_soil[k  ]))
    root_frac[0] = 1-np.sum(root_frac)

    # Radiation profiles for RRTMGP
    co2 = 348.e-6
    ch4 = 1650.e-9
    n2o = 306.e-9
    n2  = 0.7808
    o2  = 0.2095

    # Background profiles on pressure levels
    z_lay = e5.z_mean [0,:]
    z_lev = e5.zh_mean[0,:]

    p_lay = e5.p_mean [0,:]
    p_lev = e5.ph_mean[0,:]

    T_lay = e5.T_mean [0,:]
    T_lev = e5.Th_mean[0,:]

    #h2o_rad = e5.qt_mean[0,:]
    h2o_rad = e5.qt_mean[0,:]  / (eps - eps*e5.qt_mean[0,:]) # VMR
    co2_rad = np.ones(e5.nfull) * co2
    ch4_rad = np.ones(e5.nfull) * ch4
    n2o_rad = np.ones(e5.nfull) * n2o
    n2_rad  = np.ones(e5.nfull) * n2
    o2_rad  = np.ones(e5.nfull) * o2
    o3_rad  = e5.o3_mean[0,:]

    # Profiles on LES grid
    h2o_atmo = e5_at_z['qt'][0,:]
    co2_atmo = np.ones(grid.kmax) * co2
    ch4_atmo = np.ones(grid.kmax) * ch4
    n2o_atmo = np.ones(grid.kmax) * n2o
    n2_atmo  = np.ones(grid.kmax) * n2
    o2_atmo  = np.ones(grid.kmax) * o2
    o3_atmo  = e5_at_z['o3'][0,:]

    # Write MicroHH input
    message('Writing forcings as LES input')

    # Copy base .init file
    nl_file = '{}.ini'.format(settings['case_name'])
    shutil.copyfile('{}.base'.format(nl_file), nl_file)

    # Read & update namelist
    nl = mht.read_namelist(nl_file)

    # Grid info
    itot = nl['grid']['itot']
    jtot = nl['grid']['jtot']
    ktot = grid.kmax

    xsize = nl['grid']['xsize']
    ysize = nl['grid']['ysize']

    dx = xsize/itot
    dy = ysize/jtot

    npx = nl['master']['npx']
    npy = nl['master']['npy']

    # Check grid decomposition
    if not check_grid_decomposition(itot, jtot, ktot, npx, npy):
        error('Invalid grid configuration/decomposition')

    nl['grid']['ktot'] = grid.kmax
    nl['grid']['zsize'] = grid.zsize
    nl['time']['endtime'] = e5.time_sec.max()
    nl['force']['fc'] = e5.fc
    nl['radiation']['lon'] = settings['central_lon']
    nl['radiation']['lat'] = settings['central_lat']

    datetime_utc = '{0:04d}-{1:02d}-{2:02d} {3:02d}:{4:02d}:{5:02d}'.format(
            start.year, start.month, start.day, start.hour, start.minute, start.second)
    nl['time']['datetime_utc'] = datetime_utc

    ## Add column locations
    #if hetero:
    #    # Locations = Cabauw, Lek, Lopik, mais W, mais E
    #    col_i = np.array([44,161,172,20,106])
    #    col_j = np.array([118,54,158,124,127])
    #else:
    #    # Idealised grid for homogeneous run
    #    ic = itot//2
    #    jc = jtot//2
    #    dd = np.array([-750//dx, 0, 750//dx]) 
    #    col_i = np.meshgrid(ic+dd, jc+dd)[0].flatten()
    #    col_j = np.meshgrid(ic+dd, jc+dd)[1].flatten()

    #col_x = (col_i+0.5)*dx
    #col_y = (col_j+0.5)*dx

    #nl['column']['coordinates[x]'] = list(col_x)
    #nl['column']['coordinates[y]'] = list(col_y)

    mht.write_namelist(nl_file, nl)

    #
    # Write NetCDF file
    #
    init_profiles = {
            'z': grid.z, 'thl': e5_at_z['thl'][0,:], 'qt': e5_at_z['qt'][0,:],
            'u': e5_at_z['u'][0,:], 'v': e5_at_z['v'][0,:], 'nudgefac': nudge_fac,
            'co2': co2_atmo, 'ch4': ch4_atmo, 'n2o': n2o_atmo, 'n2': n2_atmo,
            'o2': o2_atmo, 'o3': o3_atmo, 'h2o': h2o_atmo}

    radiation  = {
            'z_lay': z_lay, 'z_lev': z_lev, 'p_lay': p_lay, 'p_lev': p_lev,
            't_lay': T_lay, 't_lev': T_lev, 'co2': co2_rad, 'ch4': ch4_rad,
            'n2o': n2o_rad, 'n2': n2_rad, 'o2': o2_rad, 'o3': o3_rad, 'h2o': h2o_rad}

    tdep_surface = {
            'time_surface': e5.time_sec, 'thl_sbot': e5.wths_mean,
            'qt_sbot': e5.wqs_mean, 'p_sbot': e5.ps_mean }

    print('NOTE: using tuned wls....')
    w_ls_fac = 0.4

    tdep_ls = {
            'time_ls': e5.time_sec, 'u_geo': e5_at_z['ug'], 'v_geo': e5_at_z['vg'],
            'w_ls': w_ls_fac*e5_at_z['wls'], 'thl_ls': e5_at_z['dtthl_advec'], 'qt_ls': e5_at_z['dtqt_advec'],
            'u_ls': e5_at_z['dtu_advec'], 'v_ls': e5_at_z['dtv_advec'],
            'thl_nudge': e5_at_z['thl'], 'qt_nudge': e5_at_z['qt'],
            'u_nudge': e5_at_z['u'], 'v_nudge': e5_at_z['v']}

    # Tuning .........................
    print('NOTE: using tuned theta soil....')
    vg = xr.open_dataset('../../data/van_genuchten_parameters.nc')
    theta_fc = float(vg.theta_fc[soil_index_era])
    theta_wp = float(vg.theta_wp[soil_index_era])
    theta_rel = (e5.theta_soil_mean[0,::-1]-theta_wp)/(theta_fc-theta_wp)
    print('Old rel. soil moisture content: ', theta_rel)
    theta_rel -= 0.2
    theta = theta_rel * (theta_fc-theta_wp) + theta_wp
    print('New rel. soil moisture content: ', theta_rel)

    #
    # Interpolate soil variables
    #
    theta_i = interpolate_soil(theta, z_soil_0, z_soil)
    t_i     = interpolate_soil(e5.T_soil_mean[0,::-1], z_soil_0, z_soil)

    soil = {'z': z_soil, 'theta_soil': theta_i,
            't_soil': t_i, 'index_soil': soil_index,
            'root_frac': root_frac}

    mht.write_NetCDF_input(
            settings['case_name'], float_type, init_profiles, tdep_surface, tdep_ls, radiation, soil)

    #
    # Move/copy/link files to work directory
    #
    path = '{0}/{1:04d}{2:02d}{3:02d}_t{4:02d}'.format(
            env['work_path'], start.year, start.month, start.day, start.hour)
    if os.path.exists(path):
        error('Work directory {} already exists!!'.format(path))
    else:
        os.makedirs(path)

    lsm_inp_path = '/home/bart/meteo/models/LASSIE/spatial_data/les_input/microhh_cabauw_scinti/'
    uhh_py_path = '/home/bart/meteo/models/microhh/python/'

    if hetero:
        to_copy = [
                '{}/data/van_genuchten_parameters.nc'.format(LS2D_root),
                '{}/util/slurm.py'.format(LS2D_root),
                '{}/c_veg.0000000'.format(lsm_inp_path),
                '{}/gD.0000000'.format(lsm_inp_path),
                '{}/index_soil.0000000'.format(lsm_inp_path),
                '{}/lai.0000000'.format(lsm_inp_path),
                '{}/lambda_skin.0000000'.format(lsm_inp_path),
                '{}/root_frac.0000000'.format(lsm_inp_path),
                '{}/rs_soil_min.0000000'.format(lsm_inp_path),
                '{}/rs_veg_min.0000000'.format(lsm_inp_path),
                '{}/theta_soil.0000000'.format(lsm_inp_path),
                '{}/t_soil.0000000'.format(lsm_inp_path),
                '{}/water_mask.0000000'.format(lsm_inp_path),
                '{}/z0m.0000000'.format(lsm_inp_path),
                '{}/z0h.0000000'.format(lsm_inp_path)]
    else:
        to_copy = [
                '{}/data/van_genuchten_parameters.nc'.format(LS2D_root),
                '{}/util/slurm.py'.format(LS2D_root)]

    to_move = [
            '{}.ini'.format(settings['case_name']),
            '{}_input.nc'.format(settings['case_name'])]

    to_link = {
            'microhh': env['microhh_bin'],
            'coefficients_lw.nc':
                '{}/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc'.format(env['rrtmgp_path']),
            'coefficients_sw.nc':
                '{}/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc'.format(env['rrtmgp_path']),
            'cloud_coefficients_lw.nc':
                '{}/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc'.format(env['rrtmgp_path']),
            'cloud_coefficients_sw.nc':
            '{}/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc'.format(env['rrtmgp_path'])}

    for f in to_copy:
        shutil.copy(f, path)
    for f in to_move:
        shutil.move(f, path)

    for dst,src in to_link.items():
        if env['link_files']:
            os.symlink(src, '{}/{}'.format(path,dst))
        else:
            shutil.copy(src, '{}/{}'.format(path,dst))
