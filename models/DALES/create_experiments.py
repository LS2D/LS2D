import datetime
import numpy as np
import netCDF4 as nc
from scipy.special import erf
import sys
import os
import shutil
from collections import OrderedDict as odict

# Add `src` subdirectory of LS2D to Python path
sys.path.append('{}/../../src/'.format(os.path.dirname(os.path.abspath(__file__))))

# Import the LS2D specific scripts
from download_ERA5 import download_ERA5
from read_ERA5     import Read_ERA
from messages      import header, message, error

# Import MicroHH specific tools
from DALES_tools import *
from soil import *


if (__name__ == '__main__'):

    # Start, end, and size of individual runs (dt)
    start = datetime.datetime(year=2016, month=8, day=4, hour=0)
    end   = datetime.datetime(year=2016, month=8, day=19, hour=0)
    dt    = datetime.timedelta(hours=24)

    # Working directory; individual cases are placed in yyyymmdd subdirectory
    workdir = '.'

    # Path to DALES binary
    dales_bin = '/Users/bart/meteo/models/dales/build/src/dales4'

    # Create cases for each individual `dt`
    date = start
    while date < end:

        # Dictionary with settings
        settings = {
            'central_lat' : 51.971,
            'central_lon' : 4.927,
            'area_size'   : 1,
            'case_name'   : 'cabauw',
            #'base_path'   : '/nobackup/users/stratum/ERA5/LS2D/',  # KNMI
            'base_path'   : '/Users/bart/meteo/data/LS2D/',   # Macbook
            #'base_path'   : '/home/scratch1/meteo_data/LS2D/',      # Arch
            #'base_path'   : '/home/bstratum/data/LS2D/',
            'start_date'  : date,
            'end_date'    : date+dt,
            'write_log'   : False,
            'ntasks'      : 1,
            'expnr'       : 1
            }

        header('Creating LES input')

        # Download the ERA5 data (or check whether it is available local)
        download_ERA5(settings)

        # Read ERA5 data, and calculate LES forcings, using +/-n_av grid point averages in ERA5
        e5 = Read_ERA(settings)
        e5.calculate_forcings(n_av=1)

        # Create stretched grid
        grid = Grid_stretched(kmax=160, dz0=20, nloc1=80, nbuf1=20, dz1=150)

        # Create nudge factor, controlling where nudging is aplied, and time scale
        nudge_fac = np.ones((e5.time.size, grid.z.size))

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
        #thlls = interp_time(grid.z, e5.z_mean, e5.dtthl_advec)
        #qtls  = interp_time(grid.z, e5.z_mean, e5.dtqt_advec )
        #uls   = interp_time(grid.z, e5.z_mean, e5.dtu_advec  )
        #vls   = interp_time(grid.z, e5.z_mean, e5.dtv_advec  )
        ug    = interp_time(grid.z, e5.z_mean, e5.ug         )
        vg    = interp_time(grid.z, e5.z_mean, e5.vg         )

        zero1  = np.zeros(grid.kmax)  # Dummy field (height)
        zero2  = np.zeros_like(thl)   # Dummy field (time, height)

        # Add radiative tendency
        #thlls += interp_time(grid.z, e5.z_mean, e5.dtthl_sw_mean)
        #thlls += interp_time(grid.z, e5.z_mean, e5.dtthl_lw_mean)

        #
        # Write DALES input
        #
        message('Writing forcings as LES input')

        # Docstring for DALES input files
        docstring = '{0} to {1}'.format(date, date+dt)

        # Write intial profiles
        tke = np.ones_like(grid.z)*0.1

        output = odict([('z (m)',        grid.z),
                        ('thl (K)',      thl[0,:]),
                        ('qt (kg kg-1)', qt[0,:]),
                        ('u (m s-1)',    u[0,:]),
                        ('v (m s-1)',    v[0,:]),
                        ('tke (m2 s-2)', tke)])

        write_profiles('prof.inp.{0:03d}'.format(settings['expnr']), output, grid.kmax, docstring)

        # Dummy initial conditions for the scalars
        output = odict([('z (m)',        grid.z),
                        ('qr (kg kg-1)', zero1),
                        ('nr (kg kg-1)', zero1)])

        write_profiles('scalar.inp.{0:03d}'.format(settings['expnr']), output, grid.kmax, docstring)

        # Write the large-scale forcings (ls_flux.inp.expnr)
        output_sfc = odict([('time', e5.time_sec),
                            ('p_s',  e5.ps_mean ),
                            ('T_s',  e5.Ts_mean)])

        output_ls  = odict([('time', e5.time_sec),
                            ('z',    grid.z     ),
                            ('ug',   ug         ),
                            ('vg',   vg         ),
                            ('wls',  w          )])

        write_forcings('ls_flux.inp.{0:03d}'.format(settings['expnr']), output_sfc, output_ls, docstring)

        # Also create non-time dependent file (lscale.inp), required by DALES (why?)
        output_ls2  = odict([('height', grid.z), ('ug', zero1), ('vg', zero1), ('wfls', zero1), \
                             ('dqtdxls', zero1), ('dqtdyls', zero1), ('dqtdtls', zero1), ('dthldt', zero1)])
        write_profiles('lscale.inp.{0:03d}'.format(settings['expnr']), output_ls2, grid.kmax, docstring)


        # Write nudging profiles (nudge.inp.expnr)
        output = odict([('z (m)',        grid.z   ),
                        ('factor (-)',   nudge_fac),
                        ('u (m s-1)',    u        ),
                        ('v (m s-1)',    v        ),
                        ('w (m s-1)',    zero2    ),
                        ('thl (K)',      thl      ),
                        ('qt (kg kg-1)', qt       )])

        write_time_profiles('nudge.inp.{0:03d}'.format(settings['expnr']), e5.time_sec, output, grid.kmax, docstring)

        # Write radiation background profiles
        nc_file = nc.Dataset('backrad.inp.{0:03d}.nc'.format(settings['expnr']), 'w')
        dims = nc_file.createDimension('lev', e5.nfull)

        p = nc_file.createVariable('lev', 'f4', ('lev'))
        T = nc_file.createVariable('T',   'f4', ('lev'))
        q = nc_file.createVariable('q',   'f4', ('lev'))

        p[:] = e5.p_mean [0,:]
        T[:] = e5.T_mean [0,:]
        q[:] = e5.qt_mean[0,:]

        nc_file.close()

        # Soil properties
        phi_soil = e5.phisoil_mean[0,:]
        T_soil   = e5.Tsoil_mean[0,:]

        if True:
            # Option to re-scale soil moisture content
            soil_in     = soil_med_fine      # ERA5 grid point soil type
            soil_out    = soil_fine          # ~Cabauw soil type
            old_phisoil = phi_soil
            phi_soil    = soil_in.rescale(old_phisoil, soil_out)

        # Update namelist
        namelist = 'namoptions.{0:03d}'.format(settings['expnr'])
        replace_namelist_value(namelist, 'iexpnr',   '{0:03d}'.format(settings['expnr']))
        replace_namelist_value(namelist, 'runtime',  dt.total_seconds())
        replace_namelist_value(namelist, 'trestart', dt.total_seconds())
        replace_namelist_value(namelist, 'xlat',     settings['central_lat'])
        replace_namelist_value(namelist, 'xlon',     settings['central_lon'])
        replace_namelist_value(namelist, 'xday',     date.timetuple().tm_yday)
        replace_namelist_value(namelist, 'xtime',    date.hour)
        replace_namelist_value(namelist, 'kmax',     grid.kmax)

        replace_namelist_value(namelist, 'tsoilav',  array_to_string(T_soil))
        replace_namelist_value(namelist, 'phiwav',   array_to_string(phi_soil))
        replace_namelist_value(namelist, 'tsoildeepav', T_soil[-1])  #????

        print('Setting soil properties for {} (input={})'.format(soil_out.name, soil_in.name))
        replace_namelist_value(namelist, 'gammasat', soil_out.gammasat)
        replace_namelist_value(namelist, 'nvg',      soil_out.nvg)
        replace_namelist_value(namelist, 'Lvg',      soil_out.lvg)
        replace_namelist_value(namelist, 'alphavg',  soil_out.alphavg)
        replace_namelist_value(namelist, 'phir',     soil_out.phir)
        replace_namelist_value(namelist, 'phi',      soil_out.phi_sat)
        replace_namelist_value(namelist, 'phiwp',    soil_out.phi_wp)
        replace_namelist_value(namelist, 'phifc',    soil_out.phi_fc)


        # Copy files to working directory
        dir = '{0}/{1:04d}{2:02d}{3:02d}_t{4:02d}'.format(workdir, date.year, date.month, date.day, date.hour)
        if os.path.exists(dir):
            error('Work directory {} already exists!!'.format(dir))
        else:
            os.makedirs(dir)

        strexp = '{0:03d}'.format(settings['expnr'])
        to_copy = ['namoptions.{}'.format(strexp), dales_bin, 'rrtmg_lw.nc', 'rrtmg_sw.nc']
        to_move = ['backrad.inp.{}.nc'.format(strexp), 'ls_flux.inp.{}'.format(strexp),
                   'lscale.inp.{}'.format(strexp), 'nudge.inp.{}'.format(strexp),
                   'prof.inp.{}'.format(strexp), 'scalar.inp.{}'.format(strexp)]

        for f in to_copy:
            shutil.copy(f, dir)
        for f in to_move:
            shutil.move(f, dir)

        # Advance time
        date += dt
