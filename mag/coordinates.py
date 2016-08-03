from collections import OrderedDict

import numpy as NP
from apexpy import Apex


def mag_parallels(date=None, parallels=range(-75, 76, 15), height=350, N=1000):
    """ ??? """
    apex = Apex(date=date)
    parallel_map = OrderedDict()
    lons = NP.linspace(-180, 180, N)
    for parallel in parallels:
        # parallel_map[parallel] = apex2geo([parallel] * N,
        #                                   lons,
        #                                   [height] * N)
        # glat, glon, _ = apex.apex2geo([parallel] * N,
        #                               lons,
        #                               [height] * N)
        glat, glon = apex.convert(parallel,
                                  lons,
                                  source='apex',
                                  dest='geo')
        # sort by geographic longitude
        glat, glon = zip(*sorted(zip(glat, glon), key=lambda x: x[1]))
        parallel_map[parallel] = glat, glon
    return parallel_map
