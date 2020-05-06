import matplotlib.pyplot as pl
import xarray as xr
from datetime import datetime

pl.close('all'); pl.ion()

#start = datetime(2018, 8, 11, 5)
#end   = datetime(2018, 8, 11, 19)

start = datetime(2016, 8, 17, 5)
end   = datetime(2016, 8, 17, 19)

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
fs = xr.open_dataset(nc, group='land_surface')

def calc_mean(ds, var):
    fs[var] = ds['c_low_veg']   * ds['{}_low_veg'.format(var)] +\
              ds['c_bare_soil'] * ds['{}_bare_soil'.format(var)] +\
              ds['c_wet_skin']  * ds['{}_wet_skin'.format(var)]

calc_mean(fs, 'H')
calc_mean(fs, 'LE')
calc_mean(fs, 'G')

if True:
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

if True:
    #
    # Surface fluxes
    #
    pl.figure(figsize=(8,3))
    pl.subplot(131)
    pl.plot(f.time, fs.H, color='tab:red', label=r'$\mu$HH')
    pl.scatter(cb_flux.index, cb_flux.H, s=10, marker='o', facecolor='none', alpha=0.5, color='k', label='Cabauw')
    pl.ylabel(r'$H$ (W m$^{-2}$)')
    pl.legend(ncol=2)

    pl.subplot(132)
    pl.plot(f.time, fs.LE, color='tab:red')
    pl.scatter(cb_flux.index, cb_flux.LE, s=10, marker='o', facecolor='none', alpha=0.5, color='k')
    pl.ylabel(r'$LE$ (W m$^{-2}$)')

    pl.subplot(133)
    pl.plot(f.time, fs.G, color='tab:red')
    pl.scatter(cb_flux.index, -cb_flux.G0, s=10, marker='o', facecolor='none', alpha=0.5, color='k')
    pl.ylabel(r'$G$ (W m$^{-2}$)')


if True:
    pl.figure(figsize=(8,3))
    pl.subplot(121)
    pl.plot(f.time, fs.t[:,-1], label='1 (top)')
    pl.plot(f.time, fs.t[:,-2], label='2')
    pl.plot(f.time, fs.t[:,-3], label='3')
    pl.plot(f.time, fs.t[:,-4], label='4 (bottom)')
    pl.ylabel(r'$T_\mathrm{soil}$ (K)')
    pl.legend()

    pl.subplot(122)
    pl.plot(f.time, fs.theta[:,-1])
    pl.plot(f.time, fs.theta[:,-2])
    pl.plot(f.time, fs.theta[:,-3])
    pl.plot(f.time, fs.theta[:,-4])
    pl.ylabel(r'$\theta_\mathrm{soil}$ (-)')

    pl.tight_layout()

