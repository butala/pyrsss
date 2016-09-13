from __future__ import division

import numpy as NP
import pylab as PL
from mpl_toolkits.axes_grid1 import make_axes_locatable


def centers2edges(x):
    """
    Given the array of (assumed) uniformly spaced coordinate centers
    *x*, return the array of cell edge coordinates.
    """
    delta = x[1] - x[0]
    return NP.linspace(x[0] - delta / 2,
                       x[-1] + delta / 2,
                       len(x) + 1)


def add_colorbar(ax, im, side='right', size='5%', pad=0.1):
    """
    Add colorbar to the axes *ax* with colors corresponding to the
    color mappable object *im*. Place the colorbar at the *side* of
    *ax* (options are `'right'`, `'left'`, `'top'`, or
    `'bottom'`). The width (or height) of the colorbar is specified by
    *size* and is relative to *ax*. Add space *pad* between *ax* and
    the colorbar. Return the colorbar instance.

    Reference: http://matplotlib.org/mpl_toolkits/axes_grid/users/overview.html
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(side, size=size, pad=pad)
    return PL.colorbar(im, cax=cax)
