import numpy as NP
import scipy as SP
import scipy.signal
import pylab as PL


def pzplot(b,
           a,
           ax=None,
           c='C0',
           guide_opts={'color':  '0.8',
                       'ls':     '--',
                       'zorder': -10}):
    """
    Create pole-zero plot of the rational transfer function defined by
    the polynomial coefficients *b* and *a*. Plot in axis *ax* (or
    create a new axis if `None`).  Plot poles and zeros in color
    *c*. Use the plot options *guide_opts* to style the unit circle
    and real and imaginary axis lines. Return the tuple containing the
    tuple of zeros, poles, and gain (see :func:`scipy.signal.tf2zpk`)
    and the axis *ax*.
    """
    if ax is None:
        ax = PL.subplot(111)
    z, p, k = SP.signal.tf2zpk(b, a)
    ax.add_patch(PL.Circle((0, 0), 1, fill=False, **guide_opts))
    ax.axhline(0, **guide_opts)
    ax.axvline(0, **guide_opts)
    ax.scatter(NP.real(z),
               NP.imag(z),
               marker=PL.matplotlib.markers.MarkerStyle('o', 'none'),
               facecolors='none',
               edgecolors=c)
    ax.scatter(NP.real(p),
               NP.imag(p),
               marker='x',
               color=c)
    xlim = PL.xlim()
    PL.xlim(xmin=max(xlim[0], -1.1))
    PL.xlim(xmax=max(xlim[1],  1.1))
    ylim = PL.ylim()
    PL.ylim(ymin=max(ylim[0], -1.1))
    PL.ylim(ymax=max(ylim[1],  1.1))
    ax.set_aspect('equal')
    PL.xlabel('Real axis')
    PL.ylabel('Imaginary axis')
    return (z, p, k), ax
