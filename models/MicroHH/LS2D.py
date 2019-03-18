import datetime
import numpy as np
from scipy.special import erf
import sys
import os

# Add `src` subdirectory of LS2D to Python path
sys.path.append('{}/../../src/'.format(os.path.dirname(os.path.abspath(__file__))))

# Import the LS2D specific scripts
from download_ERA5 import download_ERA5
from read_ERA5     import Read_ERA
from messages      import header, message

# Import MicroHH specific tools
import microhh_tools as mht

if (__name__ == '__main__'):

    # Dictionary with settings
    settings = {
        'central_lat' : 51.971,
        'central_lon' : 4.927,
        'area_size'   : 1,
        'case_name'   : 'cabauw',
        #'base_path'   : '/nobackup/users/stratum/ERA5/LS2D/',  # KNMI
        #'base_path'   : '/Users/bart/meteo/data/ERA5/LS2D/',   # Macbook
        'base_path'   : '/home/scratch1/meteo_data/LS2D/',      # Arch
        'start_date'  : datetime.datetime(year=2018, month=8, day=11, hour=4),
        'end_date'    : datetime.datetime(year=2018, month=8, day=11, hour=20),
        'write_log'   : True
        }

    microhh_version = 1     # {1,2}

    header('Creating LES input')

    # Download the ERA5 data (or check whether it is available local)
    download_ERA5(settings)

    # Read ERA5 data, and calculate LES forcings, using +/-n_av grid point averages in ERA5
    e5 = Read_ERA(settings)
    e5.calculate_forcings(n_av=1)

    # Read MicroHH namelist and create stretched vertical grid

    nl   = mht.Read_namelist()
    grid = mht.Stretched_grid(nl['grid']['ktot'], 130, 30, 20, 150)     # 128 hrv grid
    #grid.plot()

    # Create nudge factor, controlling where nudging is aplied, and time scale
    tau_nudge = 7200         # Nudge time scale (s)
    z0_nudge  = 1500         # Starting height of nudging (m)
    dz_nudge  = 500          # Transition thickness
    nudge_fac = 0.5 + 0.5*erf((grid.z-z0_nudge)/(0.25*dz_nudge))  # Nudge factor (0-1)
    nudge_fac /= tau_nudge   # Nudge factor (1/s)

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
    thlls += interp_time(grid.z, e5.z_mean, e5.dtthl_sw_mean)
    thlls += interp_time(grid.z, e5.z_mean, e5.dtthl_lw_mean)

    # ----------------------
    # Write MicroHH input
    # ----------------------
    message('Writing forcings as LES input')

    # 1. Update namelist variables
    mht.replace_namelist_value('zsize',   grid.zsize)
    mht.replace_namelist_value('endtime', e5.time_sec.max())
    mht.replace_namelist_value('fc',      e5.fc)
    mht.replace_namelist_value('utrans',  u.mean())
    mht.replace_namelist_value('vtrans',  v.mean())
    mht.replace_namelist_value('z0m',     e5.z0m_mean[0])
    mht.replace_namelist_value('z0h',     e5.z0h_mean[0])

    # 2. Initial profiles
    variables = {'z':grid.z, 'thl':thl[0,:], 'qt':qt[0,:], 'u':u[0,:], 'v':v[0,:], 'nudgefac':nudge_fac}
    mht.write_output('testbed.prof', variables, grid.z.size)

    # 3. Time varying nudging profiles
    # Switch between old file names 1.x, and new ones in 2.0
    if microhh_version == 1:
        names = {'thl_nudge': 'thlnudge.timeprof', 'qt_nudge': 'qtnudge.timeprof',
                 'u_nudge': 'unudge.timeprof', 'v_nudge': 'vnudge.timeprof',
                 'thl_ls': 'thlls.timeprof', 'qt_ls': 'qtls.timeprof',
                 'u_ls': 'uls.timeprof', 'v_ls': 'vls.timeprof',
                 'ug': 'ug.timeprof', 'vg': 'vg.timeprof', 'wls': 'wls.timeprof'}

    elif microhh_version == 2:
        names = {'thl_nudge': 'thl_nudge.time', 'qt_nudge': 'qt_nudge.time',
                 'u_nudge': 'u_nudge.time', 'v_nudge': 'v_nudge.time',
                 'thl_ls': 'thl_ls.time', 'qt_ls': 'qt_ls.time',
                 'u_ls': 'u_ls.time', 'v_ls': 'v_ls.time',
                 'ug': 'u_geo.time', 'vg': 'v_geo.time', 'wls': 'w_ls.time'}

    mht.write_time_profs(names['thl_nudge'], grid.z, e5.time_sec, thl)
    mht.write_time_profs(names['qt_nudge'],  grid.z, e5.time_sec, qt )
    mht.write_time_profs(names['u_nudge'],   grid.z, e5.time_sec, u  )
    mht.write_time_profs(names['v_nudge'],   grid.z, e5.time_sec, v  )

    # 4. Time varying large scale advective tendencies
    mht.write_time_profs(names['thl_ls'], grid.z, e5.time_sec, thlls)
    mht.write_time_profs(names['qt_ls'],  grid.z, e5.time_sec, qtls )
    mht.write_time_profs(names['u_ls'],   grid.z, e5.time_sec, uls  )
    mht.write_time_profs(names['v_ls'],   grid.z, e5.time_sec, vls  )

    # 5. Time varying geostrophic wind and subsidence
    mht.write_time_profs(names['ug'],  grid.z, e5.time_sec, ug)
    mht.write_time_profs(names['vg'],  grid.z, e5.time_sec, vg)
    mht.write_time_profs(names['wls'], grid.z, e5.time_sec, w)

    # 7. Time varying surface variables
    if microhh_version == 1:
        variables = {'t': e5.time_sec, 'sbot[thl]': e5.wth_mean, 'sbot[qt]': e5.wq_mean, 'pbot': e5.ps_mean}
        mht.write_output('testbed.time', variables, e5.time_sec.size)

    elif microhh_version == 2:
        mht.write_output('thl_sbot.time', {'time': e5.time_sec, 'thl_sbot': e5.wth_mean}, e5.time_sec.size)
        mht.write_output('qt_sbot.time',  {'time': e5.time_sec, 'qt_sbot': e5.wq_mean}, e5.time_sec.size)
        mht.write_output('p_sbot.time',   {'time': e5.time_sec, 'p_sbot': e5.ps_mean}, e5.time_sec.size)
