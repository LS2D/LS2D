import datetime
import numpy as np
import netCDF4 as nc
from scipy.special import erf
import sys
import os
import shutil

# Add `src` subdirectory of LS2D to Python path
sys.path.append('{}/../../src/'.format(os.path.dirname(os.path.abspath(__file__))))

# Import the LS2D specific scripts
from download_ERA5 import download_ERA5
from read_ERA5     import Read_ERA
from messages      import header, message, error

# Import MicroHH specific tools
import microhh_tools as mht

if (__name__ == '__main__'):

    # Start, end, and size of individual runs (dt)
    start = datetime.datetime(year=2018, month=8, day=11,  hour=5)
    end   = datetime.datetime(year=2018, month=8, day=11,  hour=19)
    dt    = datetime.timedelta(hours=24)

    # Working directory; individual cases are placed in yyyymmdd subdirectory
    #workdir = '/scratch-shared/bstratum/cabauw_aug2016/'
    workdir = '.'

    # Path to MicroHH binary
    microhh_bin = '/Users/bart/meteo/models/microhh2/build_dp_cpu/microhh'

    # Dictionary with settings
    settings = {
        'central_lat' : 51.971,
        'central_lon' : 4.927,
        'area_size'   : 1,
        'case_name'   : 'cabauw',
        #'base_path'   : '/nobackup/users/stratum/ERA5/LS2D/',  # KNMI
        #'base_path'   : '/Users/bart/meteo/data/LS2D/',   # Macbook
        #'base_path'   : '/home/scratch1/meteo_data/LS2D/',      # Arch
        #'base_path'   : '/home/bstratum/data/LS2D/',
        'base_path'   : '/Users/bart/meteo/data/ERA5/LS2D/',
        'start_date'  : start,
        'end_date'    : end,
        'write_log'   : False,
        'data_source' : 'MARS',
        'ntasks'      : 1
        }

    # MicroHH data type {f4, f8}
    float_type  = 'f8'

    header('Creating LES input')

    # Download the ERA5 data (or check whether it is available local)
    download_ERA5(settings)

    # Read ERA5 data, and calculate LES forcings, using +/-n_av grid point averages in ERA5
    e5 = Read_ERA(settings)
    e5.calculate_forcings(n_av=1)

    # Read MicroHH namelist and create stretched vertical grid
    nl   = mht.Read_namelist()
    grid = mht.Stretched_grid(kmax=128, nloc1=80, nbuf1=20, dz1=25, dz2=300)    # Same as DALES testbed
    #grid = mht.Stretched_grid(kmax=512, nloc1=150, nbuf1=60, dz1=10, dz2=100)    # Same as DALES testbed
    grid.plot()

    # Create nudge factor, controlling where nudging is aplied, and time scale
    tau_nudge = 10800        # Nudge time scale (s)
    #z0_nudge  = 1500         # Starting height of nudging (m)
    #dz_nudge  = 500          # Transition thickness
    #nudge_fac = 0.5 + 0.5*erf((grid.z-z0_nudge)/(0.25*dz_nudge))  # Nudge factor (0-1)
    #nudge_fac /= tau_nudge   # Nudge factor (1/s)
    nudge_fac = np.ones(grid.z.size) / tau_nudge

    # Interpolate ERA5 onto LES grid
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

    #
    # Radiation input
    #
    co2 = 348.e-6
    ch4 = 1650.e-9
    n2o = 306.e-9
    n2  = 0.7808
    o2  = 0.2095

    g1  = 3.6478
    g2  = 0.83209
    g3  = 11.3515

    # Background profiles
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

    p_hpa = p_lay/100.
    o3_rad = g1 * p_hpa**g2 * np.exp(-p_hpa/g3) * 1e-6

    # Profiles in LES domain
    h2o_atmo = qt[0,:]
    co2_atmo = np.ones(grid.kmax) * co2
    ch4_atmo = np.ones(grid.kmax) * ch4
    n2o_atmo = np.ones(grid.kmax) * n2o
    n2_atmo  = np.ones(grid.kmax) * n2
    o2_atmo  = np.ones(grid.kmax) * o2

    p_hpa = p[0,:]/100.
    o3_atmo = g1 * p_hpa**g2 * np.exp(-p_hpa/g3) * 1e-6

    # Surface / soil
    z_soil = np.array([-1.945, -0.64, -0.175, -0.035])
    soil_index = np.ones_like(z_soil)*2

    # ----------------------
    # Write MicroHH input
    # ----------------------
    message('Writing forcings as LES input')

    # 1. Update namelist variables
    mht.replace_namelist_value('ktot', grid.kmax, 'grid')
    mht.replace_namelist_value('zsize', grid.zsize, 'grid')
    mht.replace_namelist_value('endtime', e5.time_sec.max(), 'time')
    mht.replace_namelist_value('fc', e5.fc, 'force')

    # 2. Write NetCDF file
    init_profiles = {'z': grid.z, 'thl': thl[0,:], 'qt': qt[0,:], 'u': u[0,:], 'v': v[0,:], 'nudgefac': nudge_fac,
                     'co2': co2_atmo, 'ch4': ch4_atmo, 'n2o': n2o_atmo, 'n2': n2_atmo, 'o2': o2_atmo, 'o3': o3_atmo, 'h2o': h2o_atmo}
    radiation  = {'z_lay': z_lay, 'z_lev': z_lev, 'p_lay': p_lay, 'p_lev': p_lev, 't_lay': T_lay, 't_lev': T_lev,
                  'co2': co2_rad, 'ch4': ch4_rad, 'n2o': n2o_rad, 'n2': n2_rad, 'o2': o2_rad, 'o3': o3_rad, 'h2o': h2o_rad}
    tdep_surface = {'time_surface': e5.time_sec, 'thl_sbot': e5.wth_mean, 'qt_sbot': e5.wq_mean, 'p_sbot': e5.ps_mean }
    tdep_ls      = {'time_ls': e5.time_sec, 'u_geo': ug, 'v_geo': vg, 'w_ls': w,
                    'thl_ls': thlls, 'qt_ls': qtls, 'u_ls': uls, 'v_ls': vls,
                    'thl_nudge': thl, 'qt_nudge': qt, 'u_nudge': u, 'v_nudge': v}
    soil = {'z': z_soil, 'theta': e5.phisoil_mean[0,::-1], 't': e5.Tsoil_mean[0,::-1], 'index': soil_index}

    mht.write_NetCDF_input('testbed', float_type, init_profiles, tdep_surface, tdep_ls, radiation, soil)

    # Copy files to working directory
    path = '{0}/{1:04d}{2:02d}{3:02d}_t{4:02d}'.format(workdir, start.year, start.month, start.day, start.hour)
    if os.path.exists(path):
        error('Work directory {} already exists!!'.format(path))
    else:
        os.makedirs(path)

    to_copy = ['testbed.ini', microhh_bin]
    to_move = ['testbed_input.nc']
    to_link = {
            'coefficients_lw.nc': 'rrtmgp-data-lw-g256-2018-12-04.nc',
            'coefficients_sw.nc': 'rrtmgp-data-sw-g224-2018-12-04.nc',
            'cloud_coefficients_lw.nc': 'rrtmgp-cloud-optics-coeffs-lw.nc',
            'cloud_coefficients_sw.nc': 'rrtmgp-cloud-optics-coeffs-sw.nc'}

    for f in to_copy:
        shutil.copy(f, path)
    for f in to_move:
        shutil.move(f, path)
    for dst,src in to_link.items():
        shutil.copy(src, '{}/{}'.format(path,dst))



