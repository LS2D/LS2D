import datetime
import numpy as np
import netCDF4 as nc
from scipy.special import erf
import sys
import os
import shutil

# Add `src` subdirectory of LS2D to Python path
sys.path.append('{}/../../src/'.format(os.path.dirname(os.path.abspath(__file__))))
sys.path.append('{}/..'.format(os.path.dirname(os.path.abspath(__file__))))

# Import the LS2D specific scripts
from download_ERA5 import download_ERA5
from read_ERA5     import Read_ERA
from messages      import header, message, error

# Import MicroHH specific tools
import microhh_tools as mht

# Start and end of individual runs
start = datetime.datetime(year=2019, month=10, day=17,  hour=6)
end   = datetime.datetime(year=2019, month=10, day=17,  hour=13)

# Working directory; individual cases are placed in yyyymmdd subdirectory
workdir = '.'

# Dictionary with settings
settings = {
    'central_lat' : 44.9579,
    'central_lon' : 25.7794,
    'area_size'   : 1,
    'case_name'   : 'romania',
    'base_path'   : '/home/scratch1/meteo_data/LS2D/',
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
grid = mht.Stretched_grid(kmax=128, nloc1=90, nbuf1=20, dz1=20, dz2=50)
grid.plot()

# Create nudge factor, controlling where nudging is aplied, and time scale
tau_nudge = 10800        # Nudge time scale (s)
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
init_profiles = {
        'z': grid.z, 'thl': thl[0,:], 'qt': qt[0,:],
        'u': u[0,:], 'v': v[0,:], 'nudgefac': nudge_fac }
tdep_surface = {
        'time_surface': e5.time_sec, 'thl_sbot': e5.wth_mean,
        'qt_sbot': e5.wq_mean, 'p_sbot': e5.ps_mean }
tdep_ls = {
        'time_ls': e5.time_sec, 'u_geo': ug, 'v_geo': vg, 'w_ls': w,
        'thl_ls': thlls, 'qt_ls': qtls, 'u_ls': uls, 'v_ls': vls,
        'thl_nudge': thl, 'qt_nudge': qt, 'u_nudge': u, 'v_nudge': v}

mht.write_NetCDF_input('memo', float_type, init_profiles, tdep_surface, tdep_ls, radiation=None, soil=None)

# Copy files to working directory
path = '{0}/{1:04d}{2:02d}{3:02d}_t{4:02d}'.format(workdir, start.year, start.month, start.day, start.hour)
if os.path.exists(path):
    error('Work directory {} already exists!!'.format(path))
else:
    os.makedirs(path)

to_copy = ['memo.ini']
to_move = ['memo_input.nc']

for f in to_copy:
    shutil.copy(f, path)
for f in to_move:
    shutil.move(f, path)
