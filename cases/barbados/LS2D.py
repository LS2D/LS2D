import xarray as xr
import numpy as np
import sys
import os
import shutil
import subprocess
from datetime import datetime, timedelta

# Add `src` subdirectory of LS2D to Python path
LS2D_root = '/home/bart/meteo/models/LS2D/'
#LS2D_root = '/home/bstratum/models/LS2D/'
sys.path.append(f'{LS2D_root}/src/')
sys.path.append(f'{LS2D_root}/util/')

# Import the LS2D/MicroHH specific scripts
from download_ERA5 import download_ERA5
from read_ERA5     import Read_ERA
from messages      import header, message, error

import microhh_tools as mht
from grid import Grid_stretched_manual, check_grid_decomposition
from slurm import submit_case

# Constants Arden-Buck equation
a_ab = 611.21
b_ab = 18.678
c_ab = 234.5
d_ab = 257.14

# Constants
T0  = 273.15
cp  = 1004.7090
Rd  = 287.0597
Rv  = 461.5250
eps = Rd/Rv

def execute(task):
   subprocess.call(task, shell=True, executable='/bin/bash')

def calc_esat(T):
    Tc = T-T0
    return a_ab * np.exp((b_ab - (Tc / c_ab)) * (Tc / (d_ab + Tc)))

def calc_qsat(T, p):
    esat = calc_esat(T)
    return eps*esat / (p-(1.-eps)*esat)

def calc_exner(p):
    return (p/1e5)**(Rd/cp)

float_type = 'f8'    # MicroHH float type ('f4', 'f8')

# Switch between different systems:
env_cartesius = {
        'system': 'cartesius',
        'era5_path': '/archive/bstratum/ERA5/',
        'work_path': '/home/bstratum/scratch/barbados/',
        'microhh_bin': '/home/bstratum/models/microhh/build_dp_cpumpi/microhh',
        'rrtmgp_path': '/home/bstratum/models/rte-rrtmgp/',
        'auto_submit': True,
        'set_lfs_stripe': True,
        'link_files': False}

env_arch = {
        'system': 'arch',
        'era5_path': '/home/scratch1/meteo_data/LS2D/',
        'work_path': '/home/scratch1/meteo_data/LS2D/MicroHH_results/barbados/rfmip_species/',
        'microhh_bin': '/home/bart/meteo/models/microhh/build_dp_cpu/microhh',
        'rrtmgp_path': '/home/bart/meteo/models/rte-rrtmgp/',
        'auto_submit': True,
        'set_lfs_stripe': False,
        'link_files': True}

env = env_arch

# Period of experiments
start_hour = 10
end_hour   = 22
run_time   = int((end_hour-start_hour)*3600)

start = datetime(year=2020, month=1, day=20, hour=start_hour)
end   = datetime(year=2020, month=2, day=20, hour=start_hour)

use_rfmip_species = True

# Location experiments Geiske
#lat = 13.3
#lon = -57.6

# Location BCO
lat = 13.163
lon = -59.429

if use_rfmip_species:
    #
    # Read RFMIP chemical species
    #
    rfm = xr.open_dataset(
            '{}/data/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc'.format(LS2D_root))
    rfm['lon'] = rfm['lon'] - 180

    # Find nearest site:
    ii = int(np.sqrt((rfm.lon-lon)**2 + (rfm.lat-lat)**2).argmin())

    # Select current day experiment and nearest site:
    rfm = rfm.sel(expt=0, site=ii)

date = start
while date < end:

    # Dictionary with settings
    settings = {
        'central_lat' : lat,
        'central_lon' : lon,
        'area_size'   : 1,
        'case_name'   : 'barbados',
        'base_path'   : env['era5_path'],
        'start_date'  : date,
        'end_date'    : date+timedelta(seconds=run_time),
        'write_log'   : False,
        'data_source' : 'MARS',
        'ntasks'      : 1
        }

    # Slurm settings
    #max_time_per_job = 3600
    #wallclocklimit = 3600

    max_time_per_job = run_time
    wallclocklimit = 999999
    partition = 'short'

    # Controls of the nudging to ERA5
    no_nudge_near_surface = True
    z0_nudge = 2000
    z1_nudge = 3000

    # Create stretched vertical grid
    heights = np.array([0, 200, 2000, 5000, 15000, 50000])
    factors = np.array([1.025, 1.011, 1.006, 1.022, 1.07])
    grid = Grid_stretched_manual(252, 10., heights, factors)

    # Checks on grid & domain decomposition
    nl = mht.read_namelist('{}.ini'.format(settings['case_name']))
    check_grid_decomposition(
            nl['grid']['itot'], nl['grid']['jtot'],
            grid.kmax, nl['master']['npx'], nl['master']['npy'])

    header('Creating LES input')

    # Download the ERA5 data (or check whether it is available local).
    download_ERA5(settings)

    # Read ERA5 data, and calculate LES forcings, using +/-`n_av` grid point averages in ERA5.
    e5 = Read_ERA(settings)
    e5.calculate_forcings(n_av=0, method='4th')

    # Calculate sea surface properties
    thl_bot = e5.sst_mean / calc_exner(e5.ps_mean)
    qt_bot  = calc_qsat(e5.sst_mean, e5.ps_mean)

    # Interpolate ERA5 variables and forcings onto LES grid
    variables = [
            'thl', 'qt', 'u', 'v', 'wls', 'p',
            'dtthl_advec', 'dtqt_advec', 'dtu_advec', 'dtv_advec',
            'ug' ,'vg' ,'o3', 'z']
    e5_at_z = e5.interpolate_to_fixed_height(variables, grid.z)

    # Create nudge factor, controlling where nudging is aplied, and time scale
    #tau_nudge = 10800  # Nudge time scale (s)
    tau_nudge = 30  # Nudge time scale (s)

    if no_nudge_near_surface:
        # Disable nudging in the ~ABL
        low  = grid.z < z0_nudge
        high = grid.z > z1_nudge
        mid  = ~(low+high)

        nudge_fac = np.zeros(grid.kmax)
        nudge_fac[low ] = 0.
        nudge_fac[high] = 1.
        nudge_fac[mid ] = (grid.z[mid]-z0_nudge)/(z1_nudge-z0_nudge)

        nudge_fac /= tau_nudge
    else:
        # Constant nudging factor with height
        nudge_fac = np.ones(grid.kmax) / tau_nudge

    #
    # Radiation profiles for RRTMGP
    #
    if use_rfmip_species:
        def get_value(ds, var):
            """
            Get variable `var` from `ds`, and scale with its units
            """
            base_var = float(ds[var])
            unit = float(ds[var].units)
            return base_var*unit

        co2     = get_value(rfm, 'carbon_dioxide_GM')
        ch4     = get_value(rfm, 'methane_GM')
        n2o     = get_value(rfm, 'nitrous_oxide_GM')
        n2      = get_value(rfm, 'nitrogen_GM')
        o2      = get_value(rfm, 'oxygen_GM')
        co      = get_value(rfm, 'carbon_monoxide_GM')
        ccl4    = get_value(rfm, 'carbon_tetrachloride_GM')
        cfc11   = get_value(rfm, 'cfc11_GM')
        cfc12   = get_value(rfm, 'cfc12_GM')
        hcfc22  = get_value(rfm, 'hcfc22_GM')
        hfc143a = get_value(rfm, 'hfc143a_GM')
        hfc125  = get_value(rfm, 'hfc125_GM')
        hfc23   = get_value(rfm, 'hfc23_GM')
        hfc32   = get_value(rfm, 'hfc32_GM')
        hfc134a = get_value(rfm, 'hfc134a_GM')
        cf4     = get_value(rfm, 'cf4_GM')
        no2     = 0.

    else:
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

    h2o_rad = e5.qt_mean[0,:]  / (eps - eps*e5.qt_mean[0,:]) # VMR
    o3_rad  = e5.o3_mean[0,:]

    # Profiles on LES grid
    h2o_atmo = e5_at_z['qt'][0,:]
    o3_atmo  = e5_at_z['o3'][0,:]

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
    nl['time']['endtime'] = run_time
    nl['force']['fc'] = e5.fc

    nl['radiation']['lon'] = settings['central_lon']
    nl['radiation']['lat'] = settings['central_lat']

    nl['boundary']['z0m'] = e5.z0m_mean.mean()
    nl['boundary']['z0h'] = e5.z0h_mean.mean()

    datetime_utc = '{0:04d}-{1:02d}-{2:02d} {3:02d}:{4:02d}:{5:02d}'.format(
            date.year, date.month, date.day, date.hour, date.minute, date.second)
    nl['time']['datetime_utc'] = datetime_utc

    nl['cross']['xy'] = '0,{0:.1f}'.format(grid.zsize)

    # Add column locations
    column_x = np.array([1300,1800,2300])
    column_y = np.array([1300,1800,2300])
    x,y = np.meshgrid(column_x, column_y)
    nl['column']['coordinates[x]'] = list(x.flatten())
    nl['column']['coordinates[y]'] = list(y.flatten())

    mht.write_namelist(nl_file, nl)

    #
    # Write NetCDF file
    #
    if use_rfmip_species:
        species = {
                'co': co,
                'ccl4': ccl4,
                'cfc11': cfc11,
                'cfc12': cfc12,
                'cfc22': hcfc22,
                'hfc143a': hfc143a,
                'hfc125': hfc125,
                'hfc23': hfc23,
                'hfc32': hfc32,
                'hfc134a': hfc134a,
                'cf4': cf4,
                'no2': no2
                }

    init_profiles = {
            'z': grid.z,
            'thl': e5_at_z['thl'][0,:],
            'qt': e5_at_z['qt'][0,:],
            'u': e5_at_z['u'][0,:],
            'v': e5_at_z['v'][0,:],
            'nudgefac': nudge_fac,
            'h2o': h2o_atmo,
            'o3': o3_atmo,
            'co2': co2,
            'n2o': n2o,
            'ch4': ch4,
            'o2': o2,
            'n2': n2
            }

    radiation  = {
            'z_lay': z_lay,
            'z_lev': z_lev,
            'p_lay': p_lay,
            'p_lev': p_lev,
            't_lay': T_lay,
            't_lev': T_lev,
            'co2': co2,
            'ch4': ch4,
            'n2o': n2o,
            'n2': n2,
            'o2': o2,
            'o3': o3_rad,
            'h2o': h2o_rad}

    if use_rfmip_species:
        init_profiles.update(species)
        radiation    .update(species)

    # Flux BC:
    #tdep_surface = {
    #        'time_surface': e5.time_sec,
    #        'thl_sbot': e5.wths_mean,
    #        'qt_sbot': e5.wqs_mean,
    #        'p_sbot': e5.ps_mean }

    # Dirichlet BC:
    tdep_surface = {
            'time_surface': e5.time_sec,
            'thl_sbot': thl_bot,
            'qt_sbot': qt_bot,
            'p_sbot': e5.ps_mean }

    tdep_ls = {
            'time_ls': e5.time_sec,
            'u_geo': e5_at_z['ug'],
            'v_geo': e5_at_z['vg'],
            'w_ls': e5_at_z['wls'],
            'thl_ls': e5_at_z['dtthl_advec'],
            'qt_ls': e5_at_z['dtqt_advec'],
            'u_ls': e5_at_z['dtu_advec'],
            'v_ls': e5_at_z['dtv_advec'],
            'thl_nudge': e5_at_z['thl'],
            'qt_nudge': e5_at_z['qt'],
            'u_nudge': e5_at_z['u'],
            'v_nudge': e5_at_z['v']}

    mht.write_NetCDF_input(settings['case_name'], float_type, init_profiles, tdep_surface, tdep_ls, radiation)

    #
    # Copy/move/link files to working directory.
    #
    path = '{0}/{1:04d}{2:02d}{3:02d}_t{4:02d}'.format(
            env['work_path'], date.year, date.month, date.day, date.hour)
    if os.path.exists(path):
        error('Work directory {} already exists!!'.format(path))
    else:
        os.makedirs(path)

    if env['set_lfs_stripe']:
        execute('lfs setstripe -c 50 {}'.format(path))

    to_copy = ['{}.ini'.format(settings['case_name']), '{}/util/slurm.py'.format(LS2D_root)]

    to_move = ['{}_input.nc'.format(settings['case_name'])]

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

    if env['auto_submit']:
        # Submit case
        nproc = nl['master']['npx']*nl['master']['npy']
        job_name = 'mhh{0:02d}{1:02d}'.format(date.month, date.day)

        submit_case(
            settings['case_name'], run_time, max_time_per_job, wallclocklimit,
            nproc, partition, path, job_name,
            'run_restart.slurm', env['auto_submit'], bypass_slurm=True)

    # Restore namelist file
    shutil.copyfile(nl_backup, nl_file)

    # Increment time
    date += timedelta(hours=24)
