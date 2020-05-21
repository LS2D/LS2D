import matplotlib.pyplot as pl
import datetime
import numpy as np
import sys
import os
import shutil
import subprocess

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

def execute(task):
   subprocess.call(task, shell=True, executable='/bin/bash')

env_cartesius = {
        'era5_path': '/archive/bstratum/ERA5/',
        'work_path': '/home/bstratum/scratch/cabauw_aug2016_small24h/',
        'microhh_bin': '/home/bstratum/models/microhh/build_dp_cpumpi/microhh',
        'rrtmgp_path': '/home/bstratum/models/rte-rrtmgp/'}

env_arch = {
        'era5_path': '/home/scratch1/meteo_data/LS2D/',
        'work_path': '.',
        'microhh_bin': '/home/bart/meteo/models/microhh/build_dp_cpu/microhh',
        'rrtmgp_path': '/home/bart/meteo/models/rte-rrtmgp/'}

# Switch between different systems:
env = env_arch

float_type  = 'f8'   # MicroHH float type ('f4', 'f8')
auto_submit = False  # Submit the case to load balancer
link_files = True    # Switch between linking or copying files

# Time of day to simulate
start_hour = 0
run_time = 24*3600

#start_hour = 8
#run_time = 4*3600

# Days in Aug 2016:
start_day = 4
end_day = 8#18

column_x = np.array([1300,1800,2300])
column_y = np.array([1300,1800,2300])

for day in range(start_day, end_day):

    start = datetime.datetime(year=2016, month=8, day=day, hour=start_hour)
    end   = start + datetime.timedelta(seconds=run_time)

    # Dictionary with settings
    settings = {
        'central_lat' : 51.971,
        'central_lon' : 4.927,
        'area_size'   : 1,
        'case_name'   : 'cabauw',
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

    # Create stretched vertical grid
    grid = Grid_stretched(kmax=228, dz0=20, nloc1=100, nbuf1=20, dz1=100, nloc2=210, nbuf2=10, dz2=500)

    # Interpolate ERA5 variables and forcings onto LES grid
    variables = [
            'thl', 'qt', 'u', 'v', 'wls', 'p',
            'dtthl_advec', 'dtqt_advec', 'dtu_advec', 'dtv_advec',
            'ug' ,'vg' ,'o3', 'z']
    e5_z = e5.interpolate_to_fixed_height(variables, grid.z)

    # Create nudge factor, controlling where nudging is aplied, and time scale
    tau_nudge = 10800  # Nudge time scale (s)
    nudge_fac = np.ones(grid.kmax) / tau_nudge

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
    h2o_atmo = e5_z['qt'][0,:]
    co2_atmo = np.ones(grid.kmax) * co2
    ch4_atmo = np.ones(grid.kmax) * ch4
    n2o_atmo = np.ones(grid.kmax) * n2o
    n2_atmo  = np.ones(grid.kmax) * n2
    o2_atmo  = np.ones(grid.kmax) * o2
    o3_atmo  = e5_z['o3'][0,:]

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

    # Add column locations
    x,y = np.meshgrid(column_x, column_y)
    nl['column']['coordinates[x]'] = list(x.flatten())
    nl['column']['coordinates[y]'] = list(y.flatten())

    mht.write_namelist(nl_file, nl)

    #
    # Write NetCDF file
    #
    init_profiles = {
            'z': grid.z, 'thl': e5_z['thl'][0,:], 'qt': e5_z['qt'][0,:], 'u': e5_z['u'][0,:],
            'v': e5_z['v'][0,:], 'nudgefac': nudge_fac, 'co2': co2_atmo, 'ch4': ch4_atmo,
            'n2o': n2o_atmo, 'n2': n2_atmo, 'o2': o2_atmo, 'o3': o3_atmo, 'h2o': h2o_atmo}

    radiation  = {
            'z_lay': z_lay, 'z_lev': z_lev, 'p_lay': p_lay, 'p_lev': p_lev,
            't_lay': T_lay, 't_lev': T_lev, 'co2': co2_rad, 'ch4': ch4_rad,
            'n2o': n2o_rad, 'n2': n2_rad, 'o2': o2_rad, 'o3': o3_rad, 'h2o': h2o_rad}

    tdep_surface = {
            'time_surface': e5.time_sec, 'thl_sbot': e5.wths_mean,
            'qt_sbot': e5.wqs_mean, 'p_sbot': e5.ps_mean }

    tdep_ls = {
            'time_ls': e5.time_sec, 'u_geo': e5_z['ug'], 'v_geo': e5_z['vg'],
            'w_ls': e5_z['wls'], 'thl_ls': e5_z['dtthl_advec'], 'qt_ls': e5_z['dtqt_advec'],
            'u_ls': e5_z['dtu_advec'], 'v_ls': e5_z['dtv_advec'],
            'thl_nudge': e5_z['thl'], 'qt_nudge': e5_z['qt'],
            'u_nudge': e5_z['u'], 'v_nudge': e5_z['v']}

    soil = {'z': z_soil, 'theta': e5.theta_soil_mean[0,::-1], 't': e5.T_soil_mean[0,::-1], 'index': soil_index}

    mht.write_NetCDF_input('cabauw', float_type, init_profiles, tdep_surface, tdep_ls, radiation, soil)

    #
    # Copy/move/link files to working directory.
    #
    path = '{0}/{1:04d}{2:02d}{3:02d}_t{4:02d}'.format(
            env['work_path'], start.year, start.month, start.day, start.hour)
    if os.path.exists(path):
        error('Work directory {} already exists!!'.format(path))
    else:
        os.makedirs(path)

    #
    # Create slurm runscript
    #
    def create_runscript(path, workdir, ntasks, wc_time, job_name):
        with open(path, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#SBATCH -p normal\n')
            f.write('#SBATCH -n {}\n'.format(ntasks))
            f.write('#SBATCH -t {}\n'.format(wc_time))
            f.write('#SBATCH --job-name={}\n'.format(job_name))
            f.write('#SBATCH --output={}/mhh-%j.out\n'.format(workdir))
            f.write('#SBATCH --error={}/mhh-%j.err\n'.format(workdir))
            f.write('#SBATCH --constraint=haswell\n\n')

            f.write('module purge\n')
            f.write('module load surfsara\n')
            f.write('module load compilerwrappers\n')
            f.write('module load 2019\n')
            f.write('module load CMake\n')
            f.write('module load intel/2018b\n')
            f.write('module load netCDF/4.6.1-intel-2018b\n')
            f.write('module load FFTW/3.3.8-intel-2018b\n\n')

            f.write('cd {}\n\n'.format(workdir))

            f.write('srun ./microhh init cabauw\n')
            f.write('srun ./microhh run cabauw\n')

    create_runscript(
            'run.slurm', path,
            nl['master']['npx']*nl['master']['npy'],
            '24:00:00', 'mhh{0:02d}{1:02d}'.format(start.month, start.day))

    to_copy = ['cabauw.ini', '../van_genuchten_parameters.nc']
    to_move = ['cabauw_input.nc', 'run.slurm']
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
        if link_files:
            os.symlink(src, '{}/{}'.format(path,dst))
        else:
            shutil.copy(src, '{}/{}'.format(path,dst))

    if auto_submit:
        execute('sbatch {}/run.slurm'.format(path))

    # Restore namelist file
    shutil.copyfile(nl_backup, nl_file)
