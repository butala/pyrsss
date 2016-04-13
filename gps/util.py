from __future__ import division

import math


def shell_mapping(el_deg,
                  h=450,
                  R_mean=6371):
    """
    ???
    """
    el_rad = math.radians(el_deg)
    return 1 / math.sqrt(1 - (R_mean * math.cos(el_rad) / (h + R_mean))**2)
