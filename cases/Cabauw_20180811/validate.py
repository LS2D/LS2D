import matplotlib.pyplot as pl
import xarray as xr
from datetime import datetime

pl.close('all'); pl.ion()

start = datetime(2018, 8, 11, 5)
end   = datetime(2018, 8, 11, 19)

#
# Read Cabauw observations
#
cb_path = '/home/scratch1/meteo_data/observations/cabauw/'
to_drop = ['valid_dates', 'time_bnds']

cb_rad = xr.open_dataset(
        '{0}/surface_radiation_lc1/cesar_surface_radiation_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(
            cb_path, start.year, start.month), drop_variables=to_drop)
cb_rad = cb_rad.to_dataframe()
cb_rad = cb_rad.loc[start:end]

cb_flux = xr.open_dataset(
        '{0}/surface_flux_lc1/cesar_surface_flux_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(
            cb_path, start.year, start.month), drop_variables=to_drop)
cb_flux = cb_flux.to_dataframe()
cb_flux = cb_flux.loc[start:end]


#
# Read MicroHH experiment
#
#nc = 'cabauw.default.0000000.nc'
nc = 'tmp.nc'

f  = xr.open_dataset(nc)
fd = xr.open_dataset(nc, group='default')
fr = xr.open_dataset(nc, group='radiation')
ft = xr.open_dataset(nc, group='thermo')

#
# Radiation
#
pl.figure(figsize=(8,5))
pl.subplot(221)
pl.plot(f.time, fr.sw_flux_dn[:,0], color='tab:red', label=r'$\mu$HH')
pl.scatter(cb_rad.index, cb_rad.SWD, s=10, marker='o', facecolor='none', alpha=0.5, color='k', label='Cabauw')
pl.ylabel(r'$SW_\mathrm{in}$ (W m$^{-2}$)')
pl.ylim(0,1200)
pl.legend(ncol=2)

pl.subplot(222)
pl.plot(f.time, fr.sw_flux_up[:,0], color='tab:red')
pl.scatter(cb_rad.index, cb_rad.SWU, s=10, marker='o', facecolor='none', alpha=0.5, color='k')
pl.ylabel(r'$SW_\mathrm{out}$ (W m$^{-2}$)')

pl.subplot(223)
pl.plot(f.time, fr.lw_flux_dn[:,0], color='tab:red')
pl.scatter(cb_rad.index, cb_rad.LWD, s=10, marker='o', facecolor='none', alpha=0.5, color='k')
pl.ylabel(r'$LW_\mathrm{in}$ (W m$^{-2}$)')

pl.subplot(224)
pl.plot(f.time, fr.lw_flux_up[:,0], color='tab:red')
pl.scatter(cb_rad.index, cb_rad.LWU, s=10, marker='o', facecolor='none', alpha=0.5, color='k')
pl.ylabel(r'$LW_\mathrm{out}$ (W m$^{-2}$)')

pl.tight_layout()

#
# Surface fluxes
#
pl.figure(figsize=(8,3))
pl.subplot(121)
pl.plot(f.time, ft.thl_flux[:,0]*ft.rhoh[:,0]*1004, color='tab:red', label=r'$\mu$HH')
pl.scatter(cb_flux.index, cb_flux.H, s=10, marker='o', facecolor='none', alpha=0.5, color='k', label='Cabauw')
pl.ylabel(r'$H$ (W m$^{-2}$)')
pl.legend(ncol=2)

pl.subplot(122)
pl.plot(f.time, ft.qt_flux[:,0]*ft.rhoh[:,0]*2.45e6, color='tab:red')
pl.scatter(cb_flux.index, cb_flux.LE, s=10, marker='o', facecolor='none', alpha=0.5, color='k')
pl.ylabel(r'$LE$ (W m$^{-2}$)')






if False:
    nc = 'cabauw.default.0000000.nc'
    #nc = 'tmp.nc'
    
    f  = xr.open_dataset(nc)
    fd = xr.open_dataset(nc, group='default')
    fr = xr.open_dataset(nc, group='radiation')
    ft = xr.open_dataset(nc, group='thermo')
    
    pl.figure(figsize=(8,6))
    pl.subplot(321)
    pl.plot(f.time, ft.thl[:,0])
    pl.ylabel(r'$\theta_1$ (K)')
    pl.subplot(322)
    pl.plot(f.time, ft.qt[:,0])
    pl.ylabel(r'$q_1$ (K)')
    pl.subplot(323)
    pl.plot(f.time, ft.thl_flux[:,0]*ft.rhoh[:,0]*1004)
    pl.ylabel(r'$H$ (W m$^{-2}$)')
    pl.subplot(324)
    pl.plot(f.time, ft.qt_flux[:,0]*ft.rhoh[:,0]*2.45e6)
    pl.ylabel(r'$LE$ (W m$^{-2}$)')
    pl.subplot(325)
    pl.plot(f.time, fr.sw_flux_dn[:,0], label=r'$SW_\mathrm{in}$')
    pl.plot(f.time, fr.sw_flux_up[:,0], label=r'$SW_\mathrm{out}$')
    pl.plot(f.time, fr.lw_flux_dn[:,0], label=r'$LW_\mathrm{in}$')
    pl.plot(f.time, fr.lw_flux_up[:,0], label=r'$LW_\mathrm{out}$')
    pl.legend()
    pl.ylabel(r'Flux (W m$^{-2}$)')
    pl.subplot(326)
    pl.plot(f.time, ft.ql_cover)
    pl.ylabel(r'cc (-)')
    pl.tight_layout()
