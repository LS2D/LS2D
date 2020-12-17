import matplotlib.pyplot as pl
import matplotlib.dates as mdates

import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime

pl.close('all'); pl.ion()

def format_ax(major, minor, fmt='%H:%M', ax=None):
    if ax is None:
        ax = pl.gca()

    minor_loc = mdates.HourLocator(interval=minor)
    major_loc = mdates.HourLocator(interval=major)

    ax.xaxis.set_minor_locator(minor_loc)
    ax.xaxis.set_major_locator(major_loc)
    ax.xaxis.set_major_formatter(mdates.DateFormatter(fmt))

    #pl.grid()

def exner(p):
    Rd = 287.04; cp = 1005.; p0 = 1e5
    return (p/p0)**(Rd/cp)

def k_cb(z):
    return {200: 0, 140: 1, 80: 2, 40: 3, 20: 4, 10: 5, 2: 6}[z]

def interp_z(array, z, z_goal):
    k0 = np.abs(z-z_goal).argmin()
    if z[k0] > z_goal and k0>0:
        k0 -= 1
    f1 = (z_goal-z[k0]) / (z[k0+1]-z[k0])
    return (1-f1)*array[:,k0] + f1*array[:,k0+1]


if __name__ == '__main__':

    start = datetime(2018, 8, 11, 5)
    end   = datetime(2018, 8, 11, 18)

    date_str = '{0:04d}-{1:02d}-{2:02d} (UTC)'.format(start.year, start.month, start.day)

    #
    # Read Cabauw observations
    #
    if 'cb_rad' not in locals():
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

        cb_tower = xr.open_dataset(
                '{0}/tower_meteo_lc1/cesar_tower_meteo_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(
                    cb_path, start.year, start.month), drop_variables=to_drop)
        cb_tower = cb_tower.loc[{'time': slice(start,end)}]

        #
        # Read MicroHH experiment
        #
        mhh_path = '/home/scratch1/meteo_data/LS2D/MicroHH_results/202005_144x144_24h/20180811/'
        nc = f'{mhh_path}/cabauw.default.0000000.nc'
        groups = ['', 'default', 'radiation', 'thermo', 'land_surface', 'tend']
        dss = []
        for group in groups:
            dss.append( xr.open_dataset(nc, group=group) )
        ds = xr.merge(dss)

        # Read column statistics
        col_ij = [52,92]
        cols1m  = []
        cols10m = []
        for i in col_ij:
            for j in col_ij:
                cols1m.append(xr.open_dataset('{0}/cabauw_column_{1:05d}_{2:05d}_0000000.nc'.format(mhh_path,i,j)))
                cols10m.append(cols1m[-1].resample(time='10min').mean())


    def plot_columns(cols, variable, surface=False):
        cc = pl.cm.RdBu(np.linspace(0,1,len(cols)))
        for i,c in enumerate(cols):
            data = c[variable][:] if surface else c[variable][:,0]
            pl.plot(c['time'], data, color=cc[i], alpha=0.4, linewidth=0.5)

    c_mhh = '#19334c'
    cols = cols1m


    if True:
        #
        # Radiation
        #
        pl.figure(figsize=(8,5))
        pl.subplot(221)
        plot_columns(cols, 'sw_flux_dn')
        pl.plot(ds.time, ds.sw_flux_dn[:,0], color=c_mhh, label=r'$\mu$HH')
        pl.scatter(cb_rad.index, cb_rad.SWD, s=10, marker='o', facecolor='none', alpha=0.5, color='k', label='Cabauw')
        pl.ylabel(r'$SW_\mathrm{in}$ (W m$^{-2}$)')
        pl.ylim(0,1000)
        pl.legend(ncol=2)
        format_ax(major=3, minor=1)
        pl.legend()

        pl.subplot(222)
        plot_columns(cols, 'sw_flux_up')
        pl.plot(ds.time, ds.sw_flux_up[:,0], color=c_mhh)
        pl.scatter(cb_rad.index, cb_rad.SWU, s=10, marker='o', facecolor='none', alpha=0.5, color='k')
        pl.ylabel(r'$SW_\mathrm{out}$ (W m$^{-2}$)')
        format_ax(major=3, minor=1)

        pl.subplot(223)
        plot_columns(cols, 'lw_flux_dn')
        pl.plot(ds.time, ds.lw_flux_dn[:,0], color=c_mhh)
        pl.scatter(cb_rad.index, cb_rad.LWD, s=10, marker='o', facecolor='none', alpha=0.5, color='k')
        pl.ylabel(r'$LW_\mathrm{in}$ (W m$^{-2}$)')
        pl.xlabel(date_str)
        format_ax(major=3, minor=1)

        pl.subplot(224)
        plot_columns(cols, 'lw_flux_up')
        pl.plot(ds.time, ds.lw_flux_up[:,0], color=c_mhh)
        pl.scatter(cb_rad.index, cb_rad.LWU, s=10, marker='o', facecolor='none', alpha=0.5, color='k')
        pl.ylabel(r'$LW_\mathrm{out}$ (W m$^{-2}$)')
        pl.xlabel(date_str)
        format_ax(major=3, minor=1)

        pl.tight_layout()
        pl.savefig('radiation.pdf')


    if True:
        #
        # Surface fluxes
        #
        pl.figure(figsize=(8,4))
        pl.subplot(131)
        plot_columns(cols, 'H', True)
        pl.plot(ds.time, ds.H, color=c_mhh, label=r'$\mu$HH')
        pl.scatter(cb_flux.index, cb_flux.H, s=10, marker='o', facecolor='none', alpha=0.5, color='k', label='Cabauw')
        pl.ylabel(r'$H$ (W m$^{-2}$)')
        pl.xlabel(date_str)
        pl.legend(ncol=1)
        format_ax(major=6, minor=1)

        pl.subplot(132)
        plot_columns(cols, 'LE', True)
        pl.plot(ds.time, ds.LE, color=c_mhh)
        pl.scatter(cb_flux.index, cb_flux.LE, s=10, marker='o', facecolor='none', alpha=0.5, color='k')
        pl.ylabel(r'$LE$ (W m$^{-2}$)')
        pl.xlabel(date_str)
        format_ax(major=6, minor=1)

        pl.subplot(133)
        plot_columns(cols, 'G', True)
        pl.plot(ds.time, ds.G, color=c_mhh)
        pl.scatter(cb_flux.index, -cb_flux.G0, s=10, marker='o', facecolor='none', alpha=0.5, color='k')
        pl.ylabel(r'$G$ (W m$^{-2}$)')
        pl.xlabel(date_str)
        format_ax(major=6, minor=1)

        pl.tight_layout()
        pl.savefig('surface_flux.pdf')


    if True:
        #
        # Atmospheric variables
        #
        T  = np.zeros((ds.time.size, cb_tower.z.size))
        qt = np.zeros((ds.time.size, cb_tower.z.size))
        U  = np.zeros((ds.time.size, cb_tower.z.size))

        z_cb = np.array([10,20,40,80,140,200])

        for k in range(z_cb.size):
            p       = interp_z(ds.phydro.values, ds.z.values, z_cb[k])
            T [:,k] = interp_z(ds.thl.values, ds.z.values, z_cb[k]) * exner(p)
            qt[:,k] = interp_z(ds.qt.values, ds.z.values, z_cb[k])
            u       = interp_z(ds.u.values, ds.z.values, z_cb[k])
            v       = interp_z(ds.v.values, ds.z.values, z_cb[k])
            U[:,k]  = np.sqrt(u**2 + v**2)

        #cc = pl.cm.PiYG(np.array([0,0.15,0.30,0.70,0.85,1]))
        cc = pl.cm.RdBu(np.array([0,0.15,0.30,0.70,0.85,1]))

        pl.figure(figsize=(8,4))
        pl.subplot(131)
        for k in range(z_cb.size):
            pl.plot(ds.time, T[:,k], color=cc[k], label=r'$z$={}m'.format(z_cb[k]))
            pl.scatter(cb_tower.time, cb_tower.TA[:,k_cb(z_cb[k])], s=10, marker='o', facecolor='none', color=cc[k])
        pl.ylabel(r'$T$ (K)')
        pl.legend()
        format_ax(major=6, minor=1)
        pl.xlabel(date_str)

        pl.subplot(132)
        for k in range(z_cb.size):
            pl.plot(ds.time, qt[:,k]*1000, color=cc[k])
            pl.scatter(cb_tower.time, cb_tower.Q[:,k_cb(z_cb[k])], s=10, marker='o', facecolor='none', color=cc[k])
        pl.ylabel(r'$q$ (g kg$^{-1}$)')
        format_ax(major=6, minor=1)
        pl.xlabel(date_str)

        pl.subplot(133)
        for k in range(z_cb.size):
            pl.plot(ds.time, U[:,k], color=cc[k])
            pl.scatter(cb_tower.time, cb_tower.F[:,k_cb(z_cb[k])], s=10, marker='o', facecolor='none', color=cc[k])
        pl.ylabel(r'$|u|$ (m s$^{-1}$)')
        format_ax(major=6, minor=1)
        pl.xlabel(date_str)

        pl.tight_layout()
        pl.savefig('time_series_atmos.pdf')




    if False:
        #
        # ASDF
        #

        pl.figure(figsize=(7,8))
        ax=pl.subplot(511)
        pl.plot(ds.time, ds.rr*1000*3600)
        pl.ylabel(r'Rain rate (mm h$^{-1}$)')
        ax.set_xticklabels([])

        ax=pl.subplot(512)
        pl.plot(ds.time, ds.c_veg, label='vegetation')
        pl.plot(ds.time, ds.c_soil, label='bare soil')
        pl.plot(ds.time, ds.c_wet, label='wet skin')
        pl.ylabel('Fraction (-)')
        ax.set_xticklabels([])
        pl.legend()
        pl.ylim(0,1)

        ax=pl.subplot(513)
        pl.plot(ds.time, ds.wl*1000)
        pl.ylabel(r'$w_1$ (mm)')
        ax.set_xticklabels([])
        pl.legend()

        ax=pl.subplot(514)
        pl.plot(ds.time, ds.LE_veg*ds.c_veg, label='vegetation')
        pl.plot(ds.time, ds.LE_soil*ds.c_soil, label='soil')
        pl.plot(ds.time, ds.LE_wet*ds.c_wet, label='wet_skin')
        pl.ylabel(r'LE (W m$^{-2}$)')
        ax.set_xticklabels([])
        pl.legend()

        pl.subplot(515)
        pl.plot(ds.time, ds.theta[:,-1])
        pl.ylabel(r'$\theta_\mathrm{soil}$ (-)')

        pl.tight_layout()









    if False:
        #
        # Soil temperature and moisture
        #
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


    if False:
        #
        # Tendencies
        #
        k = 0

        terms = ['advec', 'diff', 'ls', 'micro', 'nudge', 'rad', 'subs', 'cor', 'pres', 'total']

        pl.figure()
        pl.subplot(221)
        for t in terms:
            var = 'thlt_{}'.format(t)
            if hasattr(ds, var):
                pl.plot(ds.time, getattr(ds, var)[:,k]*3600, label=t)
        pl.legend()
        pl.ylabel(r'$\partial_t \theta_l$ (K h$^{-1}$)')

        pl.subplot(222)
        for t in terms:
            var = 'qtt_{}'.format(t)
            if hasattr(ds, var):
                pl.plot(ds.time, getattr(ds, var)[:,k]*3600000, label=t)
        pl.legend()
        pl.ylabel(r'$\partial_t q_t$ (kg kg$^{-1}$ h$^{-1}$)')

        pl.subplot(223)
        for t in terms:
            var = 'ut_{}'.format(t)
            if hasattr(ds, var):
                pl.plot(ds.time, getattr(ds, var)[:,k]*3600, label=t)
        pl.legend()
        pl.ylabel(r'$\partial_t u$ (m s$^{-1}$ h$^{-1}$)')

        pl.subplot(224)
        for t in terms:
            var = 'vt_{}'.format(t)
            if hasattr(ds, var):
                pl.plot(ds.time, getattr(ds, var)[:,k]*3600, label=t)
        pl.legend()
        pl.ylabel(r'$\partial_t v$ (m s$^{-1}$ h$^{-1}$)')
