from collections import OrderedDict

import numpy as np
from apexpy import Apex


def mag_parallels(date, parallels=range(-75, 76, 15), height=350, N=1000):
    """
    Return a mapping between magnetic latitudes specified by
    *parallels* to the tuple of mapped geographic latitudes and
    longitudes. The mapping is made across *N* uniformly spaced
    geographic longitudes, on :class:`datetime` *date*, and at
    *height* (in [km]) in apex geomagnetic coordinates. If *date* is
    None, use the current date and time in the coordinate
    transformation.
    """
    apex = Apex(date=date)
    parallel_map = OrderedDict()
    lons = np.linspace(-180, 180, N)
    for parallel in parallels:
        glat, glon = apex.convert(parallel,
                                  lons,
                                  source='apex',
                                  dest='geo')
        # sort by geographic longitude
        glat, glon = zip(*sorted(zip(glat, glon), key=lambda x: x[1]))
        parallel_map[parallel] = glat, glon
    return parallel_map
