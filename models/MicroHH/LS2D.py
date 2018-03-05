import datetime
import numpy as np
from scipy.special import erf

# Add `src` subdirectory of LS2D to Python path
import sys; sys.path.append('/usr/people/stratum/meteo/models/LS2D/src/')

# Import the LS2D specific scripts
from download_ERA5 import download_ERA5_period
from read_ERA5     import Read_ERA
from messages      import *

# Import MicroHH specific tools
import microhh_tools as mht

def interp_time(z, ze, arr):
    out = np.empty((arr.shape[0], z.size))
    for i in range(arr.shape[0]):
        out[i,:] = np.interp(z, ze[i,:], arr[i,:])
    return out

if (__name__ == '__main__'):
    # ------------------------
    # Settings
    # ------------------------
    central_lat  = 51.971
    central_lon  = 4.927
    area_size    = 2   # ERA5 area size (lat+/-size, lon+/-size degrees)
    case_name    = 'cabauw'
    ERA5_path    = '/nobackup/users/stratum/ERA5/LS2D/'

    # Start and end date/time of experiment. For now, limited to full hours
    start_date   = datetime.datetime(year=2016, month=5, day=4, hour=5)
    end_date     = datetime.datetime(year=2016, month=5, day=4, hour=18)
    # ------------------------
    # End settings
    # ------------------------

    # Download the ERA5 data (or check whether it is available local)
    download_ERA5_period(start_date, end_date, central_lat, central_lon, area_size, ERA5_path, case_name)

    # Read ERA5 data
    e5 = Read_ERA(start_date, end_date, central_lat, central_lon, ERA5_path, case_name)

    # Calculate LES forcings, using +/-n_av grid point averages in ERA5
    e5.calculate_forcings(n_av=1)

    # Calculate LES time (seconds since start of experiment)
    LES_time = (e5.time - e5.time[0])*3600.

    # Read MicroHH namelist and create stretched vertical grid
    header('Creating LES input')

    nl   = mht.Read_namelist()
    grid = mht.Stretched_grid(nl['grid']['ktot'], 90, 20, 20, 250)     # 128 hrv grid

    # Create nudge factor, controlling where nudging is aplied, and time scale
    tau_nudge = 7200         # Nudge time scale (s)
    z0_nudge  = 1500         # Starting height of nudging (m)
    dz_nudge  = 500          # Transition thickness
    nudge_fac = 0.5 + 0.5*erf((grid.z-z0_nudge)/(0.25*dz_nudge))  # Nudge factor (0-1)
    nudge_fac /= tau_nudge   # Nudge factor (1/s)

    # Interpolate ERA5 onto LES grid
    thl   = interp_time(grid.z, e5.z_mean, e5.thl_mean )
    qt    = interp_time(grid.z, e5.z_mean, e5.qt_mean  )
    u     = interp_time(grid.z, e5.z_mean, e5.u_mean   )
    v     = interp_time(grid.z, e5.z_mean, e5.v_mean   )
    w     = interp_time(grid.z, e5.z_mean, e5.wls_mean )
    thlls = interp_time(grid.z, e5.z_mean, e5.thl_advec)
    qtls  = interp_time(grid.z, e5.z_mean, e5.qt_advec )
    uls   = interp_time(grid.z, e5.z_mean, e5.u_advec  )
    vls   = interp_time(grid.z, e5.z_mean, e5.v_advec  )
    ug    = interp_time(grid.z, e5.z_mean, e5.ug       )
    vg    = interp_time(grid.z, e5.z_mean, e5.vg       )

    # ----------------------
    # Write MicroHH input
    # ----------------------
    # 1. Initial profiles
    variables = {'z':grid.z, 'thl':thl[0,:], 'qt':qt[0,:], 'u':u[0,:], 'v':v[0,:], 'nudgefac':nudge_fac}
    mht.write_output('testbed.prof', variables, grid.z.size)

    # 2. Time varying nudging profiles
    mht.write_time_profs('thlnudge.timeprof', grid.z, LES_time, thl)
    mht.write_time_profs('qtnudge.timeprof',  grid.z, LES_time, qt )
    mht.write_time_profs('unudge.timeprof',   grid.z, LES_time, u  )
    mht.write_time_profs('vnudge.timeprof',   grid.z, LES_time, v  )

    # 3. Time varying large scale advective tendencies
    mht.write_time_profs('thlls.timeprof', grid.z, LES_time, thlls)
    mht.write_time_profs('qtls.timeprof',  grid.z, LES_time, qtls )
    mht.write_time_profs('uls.timeprof',   grid.z, LES_time, uls  )
    mht.write_time_profs('vls.timeprof',   grid.z, LES_time, vls  )

    # 4. Time varying geostrophic wind
    mht.write_time_profs('ug.timeprof', grid.z, LES_time, ug)
    mht.write_time_profs('vg.timeprof', grid.z, LES_time, vg)

    # 5. Time varying subsidence
    mht.write_time_profs('wls.timeprof', grid.z, LES_time, w)

    # 6. Time varying surface variables
    variables = {'t':LES_time, 'sbot[thl]':e5.wth_mean, 'sbot[qt]':e5.wq_mean, 'pbot':e5.ps_mean}
    mht.write_output('testbed.time', variables, LES_time.size)

    # 7. Update namelist variables
    mht.replace_namelist_value('zsize',   grid.zsize)
    mht.replace_namelist_value('endtime', LES_time.max())
    mht.replace_namelist_value('fc',      e5.fc)
    mht.replace_namelist_value('utrans',  u.mean())
    mht.replace_namelist_value('vtrans',  v.mean())
    mht.replace_namelist_value('z0m',     e5.z0m_mean[0])
    mht.replace_namelist_value('z0h',     e5.z0h_mean[0])
