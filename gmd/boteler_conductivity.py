from __future__ import division

from collections import OrderedDict

from intervals import FloatInterval

"""
Conductivity models given in D. H. Boteler and R. J. Pirjola, "The
complex-image method for calculating the magnetic and electric fields
produced at the surface of the Earth by the auroral electrojet,"
Geophysics Journal International, Vol. 132, pp. 31--40, 1998.
"""

"""
Note there are more 1-D models available in R. Pirjola, D. Boteler,
and L. Trichtchenko, "Ground effects of space weather investigated by
the surface impedance," Earth, Planets and Space, Vol. 61,
pp. 249--261, 2009.
"""


INF = float('inf')


"""
From Figure 2. Units are [m] and [Ohm/m].
"""
QUEBEC_RESISTIVITY = [(15e3,  20e3),
                      (10e3,   2e2),
                      (125e3,  1e3),
                      (200e3,  1e2),
                      (INF,      3)]


"""
From Figure 2. Units are [m] and [Ohm/m].
"""
BC_RESISTIVITY = [(4e3,   500),
                  (6e3,   150),
                  (5e3,    20),
                  (20e3,  100),
                  (65e3,  300),
                  (300e3, 100),
                  (200e3,  10),
                  (INF,     1)]


def parse_resistivity(model):
    """
    Given a Boteler resitivity specification *model* (list of depth
    [m] / resistivity [Ohm/m] tuples), return a conductivity map
    suitable for func:`conductivity.surface_impedance_1D`.
    """
    conductivity_model = OrderedDict()
    last_depth = 0
    for depth_i, r in model:
        bound = FloatInterval.closed_open(last_depth,
                                          last_depth + depth_i)
        conductivity_model[bound] = 1 / r
        last_depth += depth_i
    return conductivity_model
