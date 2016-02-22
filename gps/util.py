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


def shell_mapping(el_deg,
                  h=450,
                  R_mean=6371):
    """
    ???
    """
    el_rad = math.radians(el_deg)
    return 1 / math.sqrt(1 - (R_mean * math.cos(el_rad) / (h + R_mean))**2)
