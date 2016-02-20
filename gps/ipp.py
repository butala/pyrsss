import numpy as np

def cnv_azel2latlon(az, el, site, ht=450):
    '''
    Function to convert from an azimuth/elevation grid to a
    latitude/longitude grid given an unwarping height and a site location.
    For elevations below the horizon, returns NaN.
    INPUTS:
        az - M x N array of azimuths to be converted [degrees]
        el - M x N array of elevation to be converted [degrees]
        site - 1 x 2 array containing [latitude, longitude] of the site
           [degrees]
        ht - scalar height to be used in the conversion [km] (default is
             450 [km])
    OUTPUTS:
        lat - M x N array of latitudes [degrees]
        lon - M x N array of longitudes [degrees]

    HISTORY:
        17-Oct-2006: Converted from IDL by Jonathan J. Makela
        (jmakela@uiuc.edu)
        03-Oct-2013: Converted from MATLAB to Python by Brian Harding
        (bhardin2@illinois.edu)
    '''
    Re = 6371.2 # radius of Earth in km

    # Convert inputs from degrees to radians
    el_r = np.radians(el)
    az_r = np.radians(az)
    lat_r = np.radians(site[0])
    lon_r = np.radians(site[1])

    # Calculate the differential angle, alpha
    temp = np.cos(el_r)/(1.+(ht/Re))
    alpha = np.arccos(temp) - el_r

    # Calculate the pierce point latitude
    temp = np.sin(lat_r) * np.cos(alpha) + np.cos(lat_r)*np.cos(az_r)*np.sin(alpha)
    lat_r = np.arcsin(temp)

    # Calculate the pierce point longitude
    temp = np.sin(alpha) * np.sin(az_r) / np.cos(lat_r)
    lon_r = np.arcsin(temp) + lon_r

    # Convert radian measurements to degrees
    lat = np.degrees(lat_r)
    lon = np.degrees(lon_r)

    return lat, lon
