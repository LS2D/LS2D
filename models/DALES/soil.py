import xarray as xr
import numpy as np
import datetime

# ------------------
# Open/read soil from ERA5
# ------------------
def open_soil(date, lon, lat, path):
    file = 'soil_{0:04d}{1:02d}.nc'.format(date.year, date.month)
    ds   = xr.open_dataset('{}/{}'.format(path, file))
    return ds.sel(longitude=lon, latitude=lat, time=date, method='nearest')

def get_Tsoil_ERA5(date, lon, lat, path):
    ds = open_soil(date, lon, lat, path)
    return np.array([ds.stl1.values, ds.stl2.values, ds.stl3.values, ds.stl4.values])

def get_phisoil_ERA5(date, lon, lat, path):
    ds = open_soil(date, lon, lat, path)
    return np.array([ds.swvl1.values, ds.swvl2.values, ds.swvl3.values, ds.swvl4.values])

def array_to_string(arr):
    return str(arr).replace('[', '').replace(']', '')

# ------------------
# Soil types
# ------------------
class Soil_type:
    def __init__(self, phi_sat, phi_fc, phi_wp, 
                 gammasat, nvg, lvg, alphavg, phir, name):

        # Save soil properties
        self.phi_sat  = phi_sat
        self.phi_fc   = phi_fc
        self.phi_wp   = phi_wp
        self.gammasat = gammasat
        self.nvg      = nvg
        self.lvg      = lvg
        self.alphavg  = alphavg
        self.phir     = phir
        self.name     = name

    def calc_f2(self, phi_in):

        f2 = (phi_in - self.phi_wp) / (self.phi_fc - self.phi_wp)

        if np.any(f2 > 1):
            print('Warning: phi > phi_fc')
        if np.any(f2 < 0):
            print('Warning: phi < phi_wp')

        f2 = np.minimum(1, np.maximum(0, f2))

        return f2

    def rescale(self, phi_in, new_type):

        f2 = self.calc_f2(phi_in)
        phi_out = new_type.phi_wp + f2 * (new_type.phi_fc - new_type.phi_wp)

        return phi_out


# ECMWF soil types (a few of them at least...) and Wosten's soil
soil_fine     = Soil_type(phi_sat=0.520, phi_fc=0.448, phi_wp=0.279, 
                          gammasat=2.87e-6, nvg=1.10, lvg=-1.977, 
                          alphavg=3.67, phir=0.01, name='ECMWF fine')

soil_med_fine = Soil_type(phi_sat=0.430, phi_fc=0.383, phi_wp=0.133, 
                          gammasat=0.26e-6, nvg=1.25, lvg=-0.588,
                          alphavg=0.83, phir=0.01, name='ECMWF medium fine')

soil_wosten   = Soil_type(phi_sat=0.590, phi_fc=0.528, phi_wp=0.320,
                          gammasat=0.52e-6, nvg=1.11, lvg=-5.90,
                          alphavg=1.95, phir=0.01, name='Wosten')






if __name__ == '__main__':
    import matplotlib.pyplot as pl
    pl.close('all')

    #start = datetime.datetime(year=2017, month=2, day=1, hour=12)
    #path  = '/nobackup/users/stratum/ERA5/soil/'

    #Tsoil = get_Tsoil_ERA5  (start, 4.9, 51.97, path)
    #qsoil = get_phisoil_ERA5(start, 4.9, 51.97, path)


    phi  = np.array([0.3106892,  0.2866719,  0.29377648, 0.35078132])
    phi2 = soil_med_fine.rescale(phi, soil_fine)
    phi3 = soil_med_fine.rescale(phi, soil_wosten)
    
    f21 = soil_med_fine.calc_f2(phi)
    f22 = soil_wosten.  calc_f2(phi2)
    
    cc = pl.cm.bwr(np.linspace(0,1,4))
    
    def scatter_phi(phi, x, ha='left'):
        pl.scatter(np.ones(3)*x, phi, marker='x', color='k')
        labels = ['wp','fc','sat']
        for i in range(3):
            pl.text(x+0.05, phi[i], labels[i], ha=ha, va='center')


    pl.figure()
    ax=pl.subplot(111)
    scatter_phi([soil_med_fine.phi_wp, soil_med_fine.phi_fc, soil_med_fine.phi_sat], 1)
    scatter_phi([soil_fine.    phi_wp, soil_fine.    phi_fc, soil_fine.    phi_sat], 2)
    scatter_phi([soil_wosten.  phi_wp, soil_wosten.  phi_fc, soil_wosten.  phi_sat], 3)
    for i in range(4):
        pl.plot([1,3], [phi[i], phi[i]], color=cc[i], label='L{} original'.format(i+1))
    for i in range(4):
        pl.plot([1,2,3], [phi[i], phi2[i], phi3[i]], '--', color=cc[i], label='L{} scaled'.format(i+1))
    ax.set_xticks(np.arange(1,4))
    ax.set_xticklabels(['ECMWF medium_fine', 'ECMWF fine', 'Wosten'])
    pl.ylabel('phi (m3/m3)')
    pl.grid()
    pl.legend(ncol=2)
