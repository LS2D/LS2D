import matplotlib.pyplot as pl
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap

import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import os

from statistics import colocate, calc_basic_stats

from matplotlib import rc
rc('text', usetex=True)
rc('font', size=12)
rc('legend', fontsize=11)
rc('text.latex', preamble=r'\usepackage{sansmathfonts}')

pl.close('all'); pl.ion()

def format_ax(fmt='%d/%m', ax=None):
    if ax is None:
        ax = pl.gca()

    minor_loc = mdates.DayLocator(bymonthday=np.arange(1,32,1))
    major_loc = mdates.DayLocator(bymonthday=np.arange(1,32,4))

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


def plot_ts_scatter(ds, row, gs, ylabel_ts):

    def letter(i):
        return chr(ord('`')+i+1)

    nrows, ncols = gs.get_geometry()

    # Calculate statistics
    stats = calc_basic_stats(ds.model, ds.obs)
    stat_str = 'mae={0:.1f}, rmse={1:.1f}, bias={2:.1f}, ia={3:.2f}'.format(
            stats['mae'], stats['rmse'], stats['bias'], stats['ia'])

    # Y limits with some margin for the stats string
    ymin = min(ds.model.min(), ds.obs.min())
    ymax = max(ds.model.max(), ds.obs.max())
    dy = ymax-ymin
    ymin -= 0.05*dy
    ymax += 0.10*dy

    # Plot time series
    ax=pl.subplot(gs[row,0])
    pl.title('{}) {}'.format(letter(row*2), stat_str), loc='left', fontsize=8)
    pl.plot(ds.index, ds.model, color=c_mhh, label=r'$\mu$HH')
    pl.scatter(ds.index, ds.obs, s=6, marker='o', facecolor='none',
            alpha=0.5, color=c_obs, label='Cabauw', rasterized=True)
    pl.ylabel(ylabel_ts)
    format_ax()
    pl.xlim(start, end)
    pl.ylim(ymin, ymax)
    if row < nrows-1:
        ax.set_xticklabels([])
    if row == 0:
        pl.legend(loc=9, ncol=3, scatterpoints=3, bbox_to_anchor=[0.5, 1.5])
    #pl.text(ds.index[30], ymax, stat_str, fontsize=8, ha='left', va='top')

    # Plot scatter
    ax=pl.subplot(gs[row,1])
    pl.title('{})'.format(letter(row*2+1)), loc='left', fontsize=8)
    vmin = min(ds.obs.min(), ds.model.min())
    vmax = max(ds.obs.max(), ds.model.max())
    ax.scatter(ds.obs, ds.model, color='tab:blue', s=1, rasterized=True)
    ax.plot([vmin, vmax], [vmin, vmax], 'k:', linewidth=1)
    ax.plot([vmin, vmax], [0,0], 'k:', linewidth=1)
    ax.plot([0,0], [vmin, vmax], 'k:', linewidth=1)
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)
    ax.set_ylabel('Model')
    if row == nrows-1:
        ax.set_xlabel('Observations')



if __name__ == '__main__':

    start = datetime(2016, 8, 1, 0)
    end   = datetime(2016, 9, 1, 0)

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

        cb_sfc = xr.open_dataset(
                '{0}/surface_meteo_lc1/cesar_surface_meteo_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(
                    cb_path, start.year, start.month), drop_variables=to_drop)
        cb_sfc = cb_sfc.to_dataframe()
        cb_sfc = cb_sfc.loc[start:end]

        cb_flux = xr.open_dataset(
                '{0}/surface_flux_lc1/cesar_surface_flux_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(
                    cb_path, start.year, start.month), drop_variables=to_drop)
        cb_flux = cb_flux.to_dataframe()
        cb_flux = cb_flux.loc[start:end]

        cb_tower = xr.open_dataset(
                '{0}/tower_meteo_lc1/cesar_tower_meteo_lc1_t10_v1.0_{1:04d}{2:02d}.nc'.format(
                    cb_path, start.year, start.month), drop_variables=to_drop)
        cb_tower = cb_tower.loc[{'time': slice(start,end)}]

        cb_nubi = xr.open_dataset(
                '{0}/cesar_nubiscope_cloudcover_la1/cesar_nubiscope_cloudcover_la1_t10_v1.0_{1:04d}{2:02d}.nc'.format(
                    cb_path, start.year, start.month), drop_variables=to_drop)
        cb_nubi = cb_nubi.loc[{'time': slice(start,end)}]

        files = []
        date = start
        while date < end:
            f = '{0}/cloudnet/{1:04d}{2:02d}{3:02d}_cabauw_classification.nc'.format(
                cb_path, date.year, date.month, date.day)
            if os.path.exists(f):
                files.append(f)
            date += timedelta(days=1)
        cb_cloud = xr.open_mfdataset(files, combine='by_coords')

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
        #mhh_path = '/home/scratch1/meteo_data/LS2D/MicroHH_results/202005_144x144_24h_newrad_nonudge/'
        #mhh_path = '/home/scratch1/meteo_data/LS2D/MicroHH_results/202005_144x144_24h_newrad/'
        mhh_path = '/home/scratch1/meteo_data/LS2D/MicroHH_results/202006_144x144_24h_full_month/'
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
                files, group=group, concat_dim='time', combine='nested',
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




    if False:
        #
        # Cloudnet classification
        #
        ql_frac = ds_therm.ql_frac.values
        cloud_base = np.zeros(ds_therm.time.size)
        for t in range(ds_therm.time.size):
            kk = np.where(ql_frac[t,:] > 0)[0]
            if len(kk) > 0:
                cloud_base[t] = ds_therm.z[kk[0]]
            else:
                cloud_base[t] = np.nan
        cloud_base = np.ma.masked_invalid(cloud_base)

        # Calculate 10min average cloudnet data
        df_cloudbase = cb_cloud.cloud_base_height.to_dataframe()
        df_cloudbase = df_cloudbase.droplevel(level=1)
        df_cloudbase = df_cloudbase.resample('10min').mean()


        pl.figure()
        pl.subplot(211)
        pl.scatter(df_cloudbase.index, df_cloudbase.cloud_base_height, s=1, color=c_obs)
        pl.plot(ds_therm.time, cloud_base, color=c_mhh)
        pl.ylim(0,3000)


        #pl.plot(ds_therm.time, ds_therm.
        #pl.pcolormesh(cb_cloud.time, cb_cloud.height, cb_cloud.target_classification.T)

    if False:
        #
        # Cloud/rain classification
        #
        def cmap_alpha(cmap):
            my_cmap = cmap(np.arange(cmap.N))
            my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
            return ListedColormap(my_cmap)

        cm_ql = cmap_alpha(pl.cm.Reds)
        cm_qi = cmap_alpha(pl.cm.Blues)
        cm_qr = cmap_alpha(pl.cm.Purples)


        pl.figure()
        pl.subplot(111)
        pl.pcolormesh(ds_therm.time, ds_therm.z/1000, np.log(ds_therm.ql.T*1000), vmin=-6, vmax=0, cmap=cm_ql)
        pl.pcolormesh(ds_therm.time, ds_therm.z/1000, np.log(ds_therm.qi.T*1000), vmin=-6, vmax=0, cmap=cm_qi)
        pl.pcolormesh(ds_therm.time, ds_therm.z/1000, np.log(ds_therm.qr.T*1000), vmin=-12, vmax=-4, cmap=cm_qr)
        pl.colorbar()
        pl.ylim(0,10)

        #pl.subplot(312)
        #pl.pcolormesh(ds_therm.time, ds_therm.z, np.log(ds_therm.qi.T))

        #pl.subplot(313)
        #pl.pcolormesh(ds_therm.time, ds_therm.z, np.log(ds_therm.qr.T))


    if True:
        #
        # Cloud cover
        # NOTE: Nubiscope low < 2100m, high > 5400, mid = between.
        #
        low  = ds_therm.z.values < 2100
        high = ds_therm.z.values > 5400
        mid  = ~(low+high)

        cc_total = ds_therm.ql_frac[:,:   ].max(axis=1)
        cc_low   = ds_therm.ql_frac[:,low ].max(axis=1)
        cc_mid   = ds_therm.ql_frac[:,mid ].max(axis=1)
        cc_high  = ds_therm.ql_frac[:,high].max(axis=1)

        ds_cover = colocate(ds_therm.time.values, ds_therm.ql_cover.values*100, cb_nubi.time.values, cb_nubi.cldcover_total.values)
        ds_total = colocate(ds_therm.time.values, cc_total*100, cb_nubi.time.values, cb_nubi.cldcover_total.values)
        ds_low   = colocate(ds_therm.time.values, cc_low*100, cb_nubi.time.values, cb_nubi.cldcover_low.values)
        ds_mid   = colocate(ds_therm.time.values, cc_mid*100, cb_nubi.time.values, cb_nubi.cldcover_middle.values)
        ds_high  = colocate(ds_therm.time.values, cc_high*100, cb_nubi.time.values, cb_nubi.cldcover_high.values)

        pl.figure(figsize=(8,8))
        gs = gridspec.GridSpec(5, 2, width_ratios=[1,0.3])

        plot_ts_scatter(ds_cover, 0, gs, r'$cc_\mathrm{proj}$ (%)')
        plot_ts_scatter(ds_total, 1, gs, r'$cc_\mathrm{total}$ (%)')
        plot_ts_scatter(ds_low,   2, gs, r'$cc_\mathrm{low}$ (%)')
        plot_ts_scatter(ds_mid,   3, gs, r'$cc_\mathrm{middle}$ (%)')
        plot_ts_scatter(ds_high,  4, gs, r'$cc_\mathrm{high}$ (%)')

        pl.tight_layout()
        pl.savefig('figs/ts_cloud_cover.pdf')



    if False:
        #
        # Radiation
        #
        ds_swd = colocate(ds_base.time.values, ds_rad.sw_flux_dn[:,0].values, cb_rad.index.values, cb_rad.SWD.values)
        ds_swu = colocate(ds_base.time.values, ds_rad.sw_flux_up[:,0].values, cb_rad.index.values, cb_rad.SWU.values)
        ds_lwd = colocate(ds_base.time.values, ds_rad.lw_flux_dn[:,0].values, cb_rad.index.values, cb_rad.LWD.values)
        ds_lwu = colocate(ds_base.time.values, ds_rad.lw_flux_up[:,0].values, cb_rad.index.values, cb_rad.LWU.values)

        pl.figure(figsize=(8,6.5))
        gs = gridspec.GridSpec(4, 2, width_ratios=[1,0.3])

        plot_ts_scatter(ds_swd, 0, gs, r'$SW_\mathrm{in}$ (W m$^{-2}$)')
        plot_ts_scatter(ds_swu, 1, gs, r'$SW_\mathrm{out}$ (W m$^{-2}$)')
        plot_ts_scatter(ds_lwd, 2, gs, r'$LW_\mathrm{in}$ (W m$^{-2}$)')
        plot_ts_scatter(ds_lwu, 3, gs, r'$LW_\mathrm{out}$ (W m$^{-2}$)')

        pl.tight_layout()
        pl.savefig('figs/ts_surface_radiation.pdf')


    if False:
        #
        # Surface fluxes
        #
        ds_H  = colocate(ds_base.time.values, ds_lsm.H.values,  cb_flux.index.values,  cb_flux.H.values)
        ds_LE = colocate(ds_base.time.values, ds_lsm.LE.values, cb_flux.index.values,  cb_flux.LE.values)
        ds_G  = colocate(ds_base.time.values, ds_lsm.G.values,  cb_flux.index.values, -cb_flux.G0.values)

        pl.figure(figsize=(8,5.5))
        gs = gridspec.GridSpec(3, 2, width_ratios=[1,0.3])

        plot_ts_scatter(ds_H,  0, gs, r'$H$ (W m$^{-2}$)')
        plot_ts_scatter(ds_LE, 1, gs, r'$LE$ (W m$^{-2}$)')
        plot_ts_scatter(ds_G,  2, gs, r'$G$ (W m$^{-2}$)')

        pl.tight_layout()
        pl.savefig('figs/ts_surface_flux.pdf')


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

        ds_T010  = colocate(ds_base.time.values, T010, cb_tower.time.values, cb_tower.TA[:,k_cb(10)])
        ds_q010  = colocate(ds_base.time.values, q010*1000, cb_tower.time.values, cb_tower.Q [:,k_cb(10)])
        ds_U010  = colocate(ds_base.time.values, U010, cb_tower.time.values, cb_tower.F [:,k_cb(10)])

        ds_T200  = colocate(ds_base.time.values, T200, cb_tower.time.values, cb_tower.TA[:,k_cb(200)])
        ds_q200  = colocate(ds_base.time.values, q200*1000, cb_tower.time.values, cb_tower.Q [:,k_cb(200)])
        ds_U200  = colocate(ds_base.time.values, U200, cb_tower.time.values, cb_tower.F [:,k_cb(200)])

        for z, ds_T, ds_q, ds_U in zip([10,200], [ds_T010, ds_T200], [ds_q010, ds_q200], [ds_U010, ds_U200]):

            pl.figure(figsize=(8,5.5))
            gs = gridspec.GridSpec(3, 2, width_ratios=[1,0.3])

            plot_ts_scatter(ds_T, 0, gs, r'$T_\mathrm{{{}m}}$ (K)'.format(z))
            plot_ts_scatter(ds_q, 1, gs, r'$q_\mathrm{{{}m}}$ (g kg$^{{-1}}$)'.format(z))
            plot_ts_scatter(ds_U, 2, gs, r'$U_\mathrm{{{}m}}$ (m s$^{{-1}}$)'.format(z))

            pl.tight_layout()
            pl.savefig('figs/ts_atmosphere_{0:03d}m.pdf'.format(z))


    if False:
        #
        # Rain / clouds/ etc.
        #
        #if True:
        #    t0 = 65
        #    t1 = 70
        #    nt = t1-t0

        #    cc = pl.cm.jet(np.linspace(0,1,nt))

        #    pl.figure(figsize=(9,8))
        #    pl.subplot(221)
        #    for t in range(nt):
        #        tt = t+t0
        #        #pl.plot(ds_therm.thl[tt,:]*exner(ds_therm.phydro[tt,:]), ds_therm.z, color=cc[t], label=ds_therm.time[tt].values)
        #        pl.plot(ds_therm.rh[tt,:], ds_therm.z, color=cc[t], label=ds_therm.time[tt].values)
        #    pl.legend(fontsize=6)
        #    pl.ylim(0,5500)
        #    #pl.xlim(290,310)
        #    pl.grid()
        #    pl.xlabel('thl')

        #    pl.subplot(222)
        #    for t in range(nt):
        #        tt = t+t0
        #        pl.plot(ds_therm.qt[tt,:]*1000, ds_therm.z, color=cc[t])
        #    pl.ylim(0,5500)
        #    pl.grid()
        #    pl.xlabel('qt (g/kg)')

        #    pl.subplot(223)
        #    for t in range(nt):
        #        tt = t+t0
        #        #pl.plot((ds_therm.ql[tt,:]+ds_therm.qi[tt,:])*1000, ds_therm.z, color=cc[t])
        #        pl.plot(ds_therm.ql[tt,:]*1000, ds_therm.z, color=cc[t])
        #        pl.plot(ds_therm.qi[tt,:]*1000, ds_therm.z, color=cc[t], dashes=[4,2])
        #    pl.ylim(0,5500)
        #    pl.grid()
        #    pl.xlabel('ql (solid), qi (dashed) (g/kg)')

        #    pl.subplot(224)
        #    for t in range(nt):
        #        tt = t+t0
        #        pl.plot(ds_therm.qr[tt,:]*1000, ds_therm.z, color=cc[t])
        #    pl.ylim(0,5500)
        #    pl.xlabel('qr (g/kg)')

        #    pl.tight_layout()
        #    pl.grid()



        #if False:
        #    t0 = datetime(2016,8,1)
        #    t1 = datetime(2016,8,2)

        #    pl.figure()
        #    ax=pl.subplot(511)
        #    pl.pcolormesh(ds_therm.time, ds_therm.z, ds_therm.thl.T, cmap=pl.cm.magma_r)
        #    pl.colorbar()
        #    pl.xlim(t0, t1)
        #    pl.ylim(0,5000)

        #    ax=pl.subplot(512)
        #    pl.pcolormesh(ds_therm.time, ds_therm.z, ds_therm.qt.T*1e3, cmap=pl.cm.magma_r)
        #    pl.colorbar()
        #    pl.xlim(t0, t1)
        #    pl.ylim(0,5000)

        #    ax=pl.subplot(513, sharex=ax)
        #    pl.pcolormesh(ds_therm.time, ds_therm.z, np.log(ds_therm.ql.T*1e3), cmap=pl.cm.magma_r, vmin=-6, vmax=-1)
        #    pl.colorbar()
        #    pl.xlim(t0, t1)
        #    pl.ylim(0,5000)

        #    ax=pl.subplot(514, sharex=ax)
        #    pl.pcolormesh(ds_therm.time, ds_therm.z, np.log(ds_therm.qi.T*1e3), cmap=pl.cm.magma_r, vmin=-6, vmax=-1)
        #    pl.colorbar()
        #    pl.xlim(t0, t1)
        #    pl.ylim(0,5000)

        #    pl.subplot(515, sharex=ax)
        #    pl.pcolormesh(ds_therm.time, ds_therm.z, np.log(ds_therm.qr.T*1e3), cmap=pl.cm.magma_r, vmin=-9, vmax=-1)
        #    pl.colorbar()
        #    pl.xlim(t0, t1)
        #    pl.ylim(0,5000)





        pl.figure()
        pl.semilogy(ds_therm.time, ds_therm.rr.values*3600)
        pl.semilogy(cb_sfc.index, cb_sfc.RAIN/600*3600, '.')


        ds_cover = colocate(ds_therm.time.values, ds_therm.ql_cover.values*100, cb_nubi.time.values, cb_nubi.cldcover_total.values)
        ds_rain  = colocate(ds_therm.time.values, ds_therm.rr.values*3600, cb_sfc.index, cb_sfc.RAIN.values/10*60)

        pl.figure(figsize=(8, 5.5))
        gs = gridspec.GridSpec(3, 2, width_ratios=[1, 0.3])

        plot_ts_scatter(ds_cover,  0, gs, r'$cc_\mathrm{total}$ (\%)')
        plot_ts_scatter(ds_rain, 1, gs, r'$rr$ (mm h$^{-1}$)')
        #plot_ts_scatter(ds_G,  2, gs, r'$G$ (W m$^{-2}$)')

        pl.tight_layout()
        pl.savefig('figs/ts_cloud_rain.png')



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


