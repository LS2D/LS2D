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
    end   = datetime.datetime(year=2018, month=8, day=11,  hour=18)
    dt    = datetime.timedelta(hours=24)

    # Working directory; individual cases are placed in yyyymmdd subdirectory
    #workdir = '/scratch-shared/bstratum/cabauw_aug2016/'
    workdir = '.'

    # Path to MicroHH binary
    microhh_bin = '/Users/bart/meteo/models/microhh2/build_serial/microhh'

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
    grid = mht.Stretched_grid(kmax=140, nloc1=80, nbuf1=20, dz1=25, dz2=250)    # Same as DALES testbed
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
    thlls = interp_time(grid.z, e5.z_mean, e5.dtthl_advec)
    qtls  = interp_time(grid.z, e5.z_mean, e5.dtqt_advec )
    uls   = interp_time(grid.z, e5.z_mean, e5.dtu_advec  )
    vls   = interp_time(grid.z, e5.z_mean, e5.dtv_advec  )
    ug    = interp_time(grid.z, e5.z_mean, e5.ug         )
    vg    = interp_time(grid.z, e5.z_mean, e5.vg         )

    # Add radiative tendency
    #thlls += interp_time(grid.z, e5.z_mean, e5.dtthl_sw_mean)
    #thlls += interp_time(grid.z, e5.z_mean, e5.dtthl_lw_mean)

    # ----------------------
    # Write MicroHH input
    # ----------------------
    message('Writing forcings as LES input')

    # 1. Update namelist variables
    mht.replace_namelist_value('zsize',   grid.zsize)
    mht.replace_namelist_value('endtime', e5.time_sec.max())
    mht.replace_namelist_value('fc',      e5.fc)

    # 2. Write NetCDF file
    init_profiles = {'z': grid.z, 'thl': thl[0,:], 'qt': qt[0,:], 'u': u[0,:], 'v': v[0,:], 'nudgefac': nudge_fac}
    tdep_surface  = {'time_surface': e5.time_sec, 'thl_sbot': e5.wth_mean, 'qt_sbot': e5.wq_mean, 'p_sbot': e5.ps_mean }
    tdep_ls       = {'time_ls': e5.time_sec, 'u_geo': ug, 'v_geo': vg, 'w_ls': w,
                     'thl_ls': thlls, 'qt_ls': qtls, 'u_ls': uls, 'v_ls': vls,
                     'thl_nudge': thl, 'qt_nudge': qt, 'u_nudge': u, 'v_nudge': v}

    mht.write_NetCDF_input('testbed', float_type, init_profiles, tdep_surface, tdep_ls)

    # Copy files to working directory
    dir = '{0}/{1:04d}{2:02d}{3:02d}_t{4:02d}'.format(workdir, start.year, start.month, start.day, start.hour)
    if os.path.exists(dir):
        error('Work directory {} already exists!!'.format(dir))
    else:
        os.makedirs(dir)

    to_copy = ['testbed.ini', microhh_bin]
    to_move = ['testbed_input.nc']
    for f in to_copy:
        shutil.copy(f, dir)
    for f in to_move:
        shutil.move(f, dir)
