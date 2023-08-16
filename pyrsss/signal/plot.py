import numpy as np
import scipy as sp
from matplotlib.patches import Circle
from matplotlib.markers import MarkerStyle


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
    z, p, k = sp.signal.tf2zpk(b, a)
    ax.add_patch(Circle((0, 0), 1, fill=False, **guide_opts))
    ax.axhline(0, **guide_opts)
    ax.axvline(0, **guide_opts)
    ax.scatter(np.real(z),
               np.imag(z),
               marker=MarkerStyle('o', 'none'),
               facecolors=c)
    ax.scatter(np.real(p),
               np.imag(p),
               marker='x',
               color=c)
    xlim = ax.get_xlim()
    ax.set_xlim(xmin=max(xlim[0], -1.1))
    ax.set_xlim(xmax=max(xlim[1],  1.1))
    ylim = ax.get_ylim()
    ax.set_ylim(ymin=max(ylim[0], -1.1))
    ax.set_ylim(ymax=max(ylim[1],  1.1))
    ax.set_aspect('equal')
    ax.set_xlabel('Real axis')
    ax.set_ylabel('Imaginary axis')
    return (z, p, k), ax


if __name__ == '__main__':
    import matplotlib.pylab as plt

    # Compare with Figure 3.5 in the 3rd edition of Discrete Time
    # Signal Processing by Oppenheim and Schafer

    b = [2, -1/6, 0]
    a = [1, -1/6, -1/6]

    fig, ax = plt.subplots(1, 1)
    pzplot(b, a, ax=ax)
    plt.show()
