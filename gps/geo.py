import math
from collections import namedtuple, Iterable

import numpy as NP
import pyproj


class Geoid(namedtuple('Geoid', 'a b inv_f')):
    pass


WGS84 = Geoid(6378137,
              6356752.314245,
              298.257223563)
"""
WGS84 ellipsoid parameters.
"""

ECEF = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
"""
ECEF WGS84 coordinate conversion object.
"""

LLA = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
"""
Geodetic WGS84 coordinate conversion object.
"""

def xyz2geodetic(x, y, z):
    """
    Convert ECEF *x*, *y*, and *z* (all in [m]) to geodetic
    coordinates (using the WGS84 ellipsoid). Return lat [deg], lon
    [deg], and alt [m]. Multiple coordinates are acceptable, i.e.,
    *x*, *y*, and *z* may be equal length containers.
    """
    lon, lat, alt = pyproj.transform(ECEF, LLA, x, y, z, radians=False)
    if isinstance(lon, Iterable):
        lon = [x - 360 if x > 180 else x for x in lon]
    else:
        if lon > 180:
            lon -= 360
    return lat, lon, alt
