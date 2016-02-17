from __future__ import division

import math
from collections import namedtuple

import numpy as NP
import pyproj


class Geoid(namedtuple('Geoid', 'a b inv_f')):
    pass


WGS84 = Geoid(6378137,
              6356752.314245,
              298.257223563)


def compute_ipp():
    """
    ???
    """


def shell_mapping(el_deg,
                  h=450,
                  R_mean=6371):
    """
    ???
    """
    el_rad = math.radians(el_deg)
    return 1 / math.sqrt(1 - (R_mean * math.cos(el_rad) / (h + R_mean))**2)


if __name__ == '__main__':
    az_deg = 67.2951969
    el_deg = 61.6965387
    # JPLM location
    stn_xyz = [-2493.304303, -4655.215676, 3565.497548]
    # EFlat = 34.8054245
    # EFLon = -115.8767082

    h_km = 450
    # h_km = 350
    # h_km = 300

    Re_km = WGS84.a / 1e3

    R = NP.linalg.norm(stn_xyz)
    stn_lon_rad = math.atan2(stn_xyz[1], stn_xyz[0])
    stn_lat_rad = math.asin(stn_xyz[2] / R)

    # print(math.degrees(stn_lat_rad))
    # print(math.degrees(stn_lon_rad))

    az_rad = math.radians(az_deg)
    el_rad = math.radians(el_deg)
    phi = math.pi / 2 - el_rad - math.sin(Re_km / (Re_km + h_km) * math.cos(el_rad))

    ipp_lat_rad = math.asin(math.sin(stn_lat_rad) * math.cos(phi) + \
                            math.cos(stn_lat_rad) * math.sin(phi) * math.cos(az_rad))

    ipp_lon_rad = stn_lon_rad + math.asin(math.sin(phi) * math.sin(az_rad) / math.cos(ipp_lat_rad))

    h_m = h_km * 1e3
    Re_m = Re_km * 1e3

    x = (Re_m + h_m) * math.cos(ipp_lat_rad) * math.cos(ipp_lon_rad)
    y = (Re_m + h_m) * math.cos(ipp_lat_rad) * math.sin(ipp_lon_rad)
    z = (Re_m + h_m) * math.sin(ipp_lat_rad)

    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    print(math.degrees(ipp_lat_rad))
    print(math.degrees(ipp_lon_rad))

    lat, lon, alt = pyproj.transform(ecef, lla, x, y, z, radians=False)


    print(lat, lon, alt)

    print(shell_mapping(el_deg))
    print(el_deg)

    print(shell_mapping(10))
    print(shell_mapping(20))
