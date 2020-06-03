import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec

import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime, timedelta

pl.close('all'); pl.ion()

def format_ax(fmt='%d/%m', ax=None):
    if ax is None:
        ax = pl.gca()

    minor_loc = mdates.HourLocator(interval=24)
    major_loc = mdates.HourLocator(interval=48)

    ax.xaxis.set_minor_locator(minor_loc)
    ax.xaxis.set_major_locator(major_loc)
    ax.xaxis.set_major_formatter(mdates.DateFormatter(fmt))

    pl.grid(which='both', axis='x')

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

def colloc(time_model, var_model, time_obs, var_obs):
    ds1 = pd.DataFrame({'model': var_model}, index=time_model)
    ds1.index = ds1.index.round('10s')

    ds2 = pd.DataFrame({'obs': var_obs}, index=time_obs)
    ds2.index = ds2.index.round('10s')

    ds = pd.concat([ds1, ds2], axis=1)
    ds.dropna(inplace=True)

    return ds

def colloc_multi(times, variables, names):
    dss = []
    for time, variable, name in zip(times, variables, names):
        ds = pd.DataFrame({name: variable}, index=time)
        ds.index = ds.index.round('10s')
        dss.append(ds)

    ds = pd.concat(dss, axis=1)
    ds.dropna(inplace=True)

    return ds

def calc_stats(model, obs):
    mae  = np.abs(model-obs).mean()
    rmse = np.sqrt(np.mean((model-obs)**2))
    bias = (model-obs).mean()
    ia   = 1 -(np.sum((obs-model)**2))/(np.sum((np.abs(model-np.mean(obs))+np.abs(obs-np.mean(obs)))**2))
    return mae, rmse, bias, ia


def scatter(obs, model, xlabel=False):
    vmin = min(obs.min(), model.min())
    vmax = max(obs.max(), model.max())
    pl.scatter(obs, model, s=1)
    pl.plot([vmin, vmax], [vmin, vmax], 'k:', linewidth=1)
    pl.plot([vmin, vmax], [0,0], 'k:', linewidth=1)
    pl.plot([0,0], [vmin, vmax], 'k:', linewidth=1)
    pl.xlim(vmin, vmax)
    pl.ylim(vmin, vmax)
    pl.ylabel('Model')
    if xlabel:
        pl.xlabel('Observations')

def plot_ts_scatter(ds, row, gs, ylabel_ts, xlabel_sc=False, legend=False):

    # Calculate statistics
    mae, rmse, bias, ia = calc_stats(ds.model, ds.obs)
    stat_str = 'mae={0:.1f}, rmse={1:.1f}, bias={2:.1f}, ia={3:.2f}'.format(mae, rmse, bias, ia)

    # Y limits with some margin for the stats string
    ymin = min(ds.model.min(), ds.obs.min())
    ymax = max(ds.model.max(), ds.obs.max())
    dy = ymax-ymin
    ymin -= 0.05*dy
    ymax += 0.10*dy

    # Plot time series
    ax=pl.subplot(gs[row,0])
    pl.plot(ds.index, ds.model, color=c_mhh, label=r'$\mu$HH')
    pl.scatter(ds.index, ds.obs, s=6, marker='o', facecolor='none', alpha=0.5, color=c_obs, label='Cabauw')
    pl.ylabel(ylabel_ts)
    format_ax()
    pl.xlim(start, end)
    pl.ylim(ymin, ymax)
    ax.set_xticklabels([])
    if legend:
        pl.legend(loc=9, ncol=3, scatterpoints=3, bbox_to_anchor=[0.5, 1.4])
    pl.text(ds.index[30], ymax, stat_str, fontsize=8, ha='left', va='top')

    # Plot scatter
    ax=pl.subplot(gs[row,1])
    vmin = min(ds.obs.min(), ds.model.min())
    vmax = max(ds.obs.max(), ds.model.max())
    ax.scatter(ds.obs, ds.model, color='tab:blue', s=1)
    ax.plot([vmin, vmax], [vmin, vmax], 'k:', linewidth=1)
    ax.plot([vmin, vmax], [0,0], 'k:', linewidth=1)
    ax.plot([0,0], [vmin, vmax], 'k:', linewidth=1)
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)
    ax.set_ylabel('Model')
    if xlabel_sc:
        ax.set_xlabel('Observations')



if __name__ == '__main__':

    start = datetime(2016, 8, 4,  0)
    end   = datetime(2016, 8, 18, 0)

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
        # Read ERA5
        #
        era5_path = '/home/scratch1/meteo_data/ERA5/LASSIE/monthly/'
        ds_era5_sfc_fc = xr.open_dataset('{}/sfc_fc/2016/ERA5_sfc_fc_201608.nc'.format(era5_path))
        ds_era5_sfc_fc = ds_era5_sfc_fc.sel(longitude=5, latitude=52., method='nearest')
        era5_sfc_fc_time = pd.to_datetime(ds_era5_sfc_fc.time.values) - timedelta(minutes=30)

        # Derived quantities
        ds_era5_sfc_fc['ssru'] = -(ds_era5_sfc_fc['ssr'] - ds_era5_sfc_fc['ssrd'])
        ds_era5_sfc_fc['stru'] = -(ds_era5_sfc_fc['str'] - ds_era5_sfc_fc['strd'])

        #
        # Read MicroHH experiments
        #
        mhh_path = '/home/scratch1/meteo_data/LS2D/MicroHH_results/202005_144x144_24h/'
        files = []
        dates = pd.date_range(start, end)[:-1]
        for date in dates:
            files.append('{0}/{1:04d}{2:02d}{3:02d}/cabauw.default.0000000.nc'.format(
                mhh_path, date.year, date.month, date.day))

        def read_group(files, group):
            to_drop = ['lw_flux_up_ref', 'lw_flux_dn_ref']

            def preproc(ds):
                return ds.isel(time=np.arange(1,ds.time.size))

            return xr.open_mfdataset(
                files, group=group, concat_dim='time',
                drop_variables=to_drop, preprocess=preproc)

        ds_base  = read_group(files, '')
        ds_def   = read_group(files, 'default')
        ds_rad   = read_group(files, 'radiation')
        ds_lsm   = read_group(files, 'land_surface')
        ds_therm = read_group(files, 'thermo')

        # Merge in the base group, to get the dimensions
        ds_def   = xr.merge([ds_def,   ds_base])
        ds_rad   = xr.merge([ds_rad,   ds_base])
        ds_lsm   = xr.merge([ds_lsm,   ds_base])
        ds_therm = xr.merge([ds_therm, ds_base])

    c_mhh = '#19334c'
    c_era = 'tab:red'
    c_obs = 'tab:blue'
    mev = 1


    if True:
        #
        # Clouds et al.
        #
        pl.figure()
        pl.subplot(311)
        pl.pcolormesh(ds_therm.time, ds_therm.z, np.log(ds_therm.ql.T))
        pl.ylim(0,10000)

        pl.subplot(312)
        pl.pcolormesh(ds_def.time, ds_def.zh, ds_def.w_2.T, vmin=0, vmax=1)
        pl.ylim(0,10000)


    if False:
        #
        # Radiation
        #
        ds_swd = colloc(ds_base.time.values, ds_rad.sw_flux_dn[:,0].values, cb_rad.index.values, cb_rad.SWD.values)
        ds_swu = colloc(ds_base.time.values, ds_rad.sw_flux_up[:,0].values, cb_rad.index.values, cb_rad.SWU.values)
        ds_lwd = colloc(ds_base.time.values, ds_rad.lw_flux_dn[:,0].values, cb_rad.index.values, cb_rad.LWD.values)
        ds_lwu = colloc(ds_base.time.values, ds_rad.lw_flux_up[:,0].values, cb_rad.index.values, cb_rad.LWU.values)

        pl.figure(figsize=(8,6.5))
        gs = gridspec.GridSpec(4, 2, width_ratios=[1,0.3])

        plot_ts_scatter(ds_swd, 0, gs, r'$SW_\mathrm{in}$ (W m$^{-2}$)', legend=True)
        plot_ts_scatter(ds_swu, 1, gs, r'$SW_\mathrm{out}$ (W m$^{-2}$)')
        plot_ts_scatter(ds_lwd, 2, gs, r'$LW_\mathrm{in}$ (W m$^{-2}$)')
        plot_ts_scatter(ds_lwu, 3, gs, r'$LW_\mathrm{out}$ (W m$^{-2}$)', xlabel_sc=True)

        pl.tight_layout()
        pl.savefig('ts_surface_radiation.pdf')


    if False:
        #
        # Surface fluxes
        #
        ds_H  = colloc(ds_base.time.values, ds_lsm.H.values,  cb_flux.index.values,  cb_flux.H.values)
        ds_LE = colloc(ds_base.time.values, ds_lsm.LE.values, cb_flux.index.values,  cb_flux.LE.values)
        ds_G  = colloc(ds_base.time.values, ds_lsm.G.values,  cb_flux.index.values, -cb_flux.G0.values)

        pl.figure(figsize=(8,5.5))
        gs = gridspec.GridSpec(3, 2, width_ratios=[1,0.3])

        plot_ts_scatter(ds_H,  0, gs, r'$H$ (W m$^{-2}$)', legend=True)
        plot_ts_scatter(ds_LE, 1, gs, r'$LE$ (W m$^{-2}$)')
        plot_ts_scatter(ds_G,  2, gs, r'$G$ (W m$^{-2}$)', xlabel_sc=True)

        pl.tight_layout()
        pl.savefig('ts_surface_flux.pdf')


    if False:
        #
        # Atmospheric variables
        #
        p010 = interp_z(ds_therm.phydro.values, ds_therm.z.values, 10)
        T010 = interp_z(ds_therm.thl.values, ds_therm.z.values, 10) * exner(p010)
        q010 = interp_z(ds_therm.qt.values, ds_therm.z.values, 10)
        u010 = interp_z(ds_def.u.values, ds_base.z.values, 10)
        v010 = interp_z(ds_def.v.values, ds_base.z.values, 10)
        U010 = np.sqrt(u010**2 + v010**2)

        p200 = interp_z(ds_therm.phydro.values, ds_therm.z.values, 200)
        T200 = interp_z(ds_therm.thl.values, ds_therm.z.values, 200) * exner(p200)
        q200 = interp_z(ds_therm.qt.values, ds_therm.z.values, 200)
        u200 = interp_z(ds_def.u.values, ds_base.z.values, 200)
        v200 = interp_z(ds_def.v.values, ds_base.z.values, 200)
        U200 = np.sqrt(u200**2 + v200**2)

        ds_T010  = colloc(ds_base.time.values, T010, cb_tower.time.values, cb_tower.TA[:,k_cb(10)])
        ds_q010  = colloc(ds_base.time.values, q010*1000, cb_tower.time.values, cb_tower.Q [:,k_cb(10)])
        ds_U010  = colloc(ds_base.time.values, U010, cb_tower.time.values, cb_tower.F [:,k_cb(10)])

        ds_T200  = colloc(ds_base.time.values, T200, cb_tower.time.values, cb_tower.TA[:,k_cb(200)])
        ds_q200  = colloc(ds_base.time.values, q200*1000, cb_tower.time.values, cb_tower.Q [:,k_cb(200)])
        ds_U200  = colloc(ds_base.time.values, U200, cb_tower.time.values, cb_tower.F [:,k_cb(200)])

        for z, ds_T, ds_q, ds_U in zip([10,200], [ds_T010, ds_T200], [ds_q010, ds_q200], [ds_U010, ds_U200]):

            pl.figure(figsize=(8,5.5))
            gs = gridspec.GridSpec(3, 2, width_ratios=[1,0.3])

            plot_ts_scatter(ds_T, 0, gs, r'$T_\mathrm{{{}m}}$ (K)'.format(z), legend=True)
            plot_ts_scatter(ds_q, 1, gs, r'$q_\mathrm{{{}m}}$ (g kg$^{{-1}}$)'.format(z))
            plot_ts_scatter(ds_U, 2, gs, r'$U_\mathrm{{{}m}}$ (m s$^{{-1}}$)'.format(z), xlabel_sc=True)

            pl.tight_layout()
            pl.savefig('ts_atmosphere_{0:03d}m.pdf'.format(z))




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


