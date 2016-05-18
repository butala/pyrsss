from __future__ import division

import numpy as NP


def centers2edges(x):
    """
    Given the array of (assumed) uniformly spaced coordinate centers
    *x*, return the array of cell edge coordinates.
    """
    delta = x[1] - x[0]
    return NP.linspace(x[0] - delta / 2,
                       x[-1] + delta / 2,
                       len(x) + 1)
