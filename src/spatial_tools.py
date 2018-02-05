import numpy as np

r_earth = 6367.47 * 1000    # Radius earth (m)

def haversine(lon1, lat1, lon2, lat2):
    """ Distance between two coordinates """
    lon1,lat1,lon2,lat2 = map(np.deg2rad,[lon1,lat1,lon2,lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    return r_earth * 2 * np.arcsin(np.sqrt(a))

def distance(lat1, lon1, lat2, lon2):
    dist_we = dlon(lon1, lon2, lat1)
    dist_ns = dlat(lat1, lat2)
    return (dist_we**2 + dist_ns**2)**0.5    

def dlon(lonW, lonE, lat):
    """ Distance between longitudes in spherical coordinates """
    return r_earth * np.cos(np.deg2rad(lat)) * np.deg2rad(lonE - lonW)

def ddlon(valW, valE, lonW, lonE, lat):
    """ Gradient between longitudes in spherical coordinates """
    return (valE - valW) / dlon(lonW, lonE, lat) 

def dlat(latS, latN):
    """ Distance between latitudes in spherical coordinates """
    return r_earth * np.deg2rad(latN-latS)

def ddlat(valS, valN, latS, latN):
    """ Gradient between latitudes in spherical coordinates """
    return (valN-valS) / dlat(latS, latN)


if __name__ == '__main__':
    """ Test / example, only executed if script is called directly """

    # Test of correct gradient signs
    latC = 50.
    latN = 51.
    latS = 49.

    lonC = 5.
    lonW = 4.
    lonE = 6.

    # Positive gradient towards east
    valE = 2.
    valW = 1.

    # Positive gradient towards north
    valS = 1.
    valN = 2.

    # Positive distances from S->N and W->E
    print('dlat=', dlat(latS, latN))
    print('dlon=', dlon(lonW, lonE, latC))

    # Positive gradients
    print('ddlat=', ddlat(valS, valN, latS, latN))
    print('ddlon=', ddlon(valW, valE, lonW, lonE, latC))
