from __future__ import division

import numpy as NP
import pylab as PL
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy
import cartopy.crs as ccrs


def centers2edges(x):
    """
    Given the array of (assumed) uniformly spaced coordinate centers
    *x*, return the array of cell edge coordinates.
    """
    delta = x[1] - x[0]
    return NP.linspace(x[0] - delta / 2,
                       x[-1] + delta / 2,
                       len(x) + 1)


def add_colorbar(ax, im, side='right', size='5%', pad=0.1, **kwds):
    """
    Add colorbar to the axes *ax* with colors corresponding to the
    color mappable object *im*. Place the colorbar at the *side* of
    *ax* (options are `'right'`, `'left'`, `'top'`, or
    `'bottom'`). The width (or height) of the colorbar is specified by
    *size* and is relative to *ax*. Add space *pad* between *ax* and
    the colorbar. The remaining keyword arguments *kwds* are passed to
    the call to :func:`colorbar`. Return the colorbar instance.

    Reference: http://matplotlib.org/mpl_toolkits/axes_grid/users/overview.html
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(side, size=size, pad=pad)
    cb = PL.colorbar(im, cax=cax, **kwds)
    PL.axes(ax)
    return cb


def plot_map(ax=None, alpha=0.3, zorder=0):
    """
    Add map features (coastlines, national boundaries, etc.) to *a*
    with transparency level *alpha* and *zorder*. Return *ax*.
    """
    if ax is None:
        ax = PL.axes(projection=ccrs.PlateCarree())
    # national boundaries
    boundaries_50m = cartopy.feature.NaturalEarthFeature(category='cultural',
                                                         name='admin_0_boundary_lines_land',
                                                         scale='50m',
                                                         edgecolor='k',
                                                         facecolor='none')
    ax.add_feature(boundaries_50m,
                   alpha=alpha,
                   zorder=zorder)
    # states
    states_50m = cartopy.feature.NaturalEarthFeature(category='cultural',
                                                     name='admin_1_states_provinces_lines',
                                                     scale='50m',
                                                     edgecolor='k',
                                                     facecolor='none')
    ax.add_feature(states_50m,
                   alpha=alpha,
                   zorder=zorder)
    # coastlines
    coastline_50m = cartopy.feature.NaturalEarthFeature('physical',
                                                        'coastline',
                                                        '50m',
                                                        edgecolor='k',
                                                        facecolor='none')
    ax.add_feature(coastline_50m,
                   alpha=alpha,
                   zorder=zorder)
    # lakes
    lakes_110m = cartopy.feature.NaturalEarthFeature('physical',
                                                     'lakes',
                                                     '110m',
                                                     edgecolor='k',
                                                     facecolor='none')
    # add all shape objects
    ax.add_feature(lakes_110m,
                   alpha=alpha,
                   zorder=zorder)

    return ax
