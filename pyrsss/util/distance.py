from __future__ import division

import numpy as NP
import pyproj

from angle import convert_lon


WGS84 = pyproj.Geod(ellps='WGS84')
"""WGS84 ellipsoid."""


def distance(lat1, lon1, lat2, lon2, geod=WGS84, units='km'):
    """
    Compute and return the geodetic distance from (*lat1*, *lon1*) to
    (*lat2*, *lon2*). Note that the coordinates may be scalars or
    equal length containers. Use ellipsoid *geod*. Return output in
    *units*.
    """
    def listify(x):
        try:
            iter(x)
            return x
        except TypeError:
            return [x]
    lat1 = listify(lat1)
    lon1 = map(convert_lon, listify(lon1))
    lat2 = listify(lat2)
    lon2 = map(convert_lon, listify(lon2))
    assert len(lat1) == len(lon1) == len(lat2) == len(lon2)
    _, _, d = geod.inv(lon1, lat1, lon2, lat2)
    d = NP.array(d)
    if units == 'km':
        d /= 1e3
    elif units == 'm':
        pass
    else:
        raise ValueError('unknown units {}'.format(units))
    if len(d) == 1:
        d = d[0]
    return d
