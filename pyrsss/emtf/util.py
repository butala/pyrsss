import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import mu_0

from .calc_e_3d import parse_xml, Zw_interpolator


def plot_3D_reponse(xml_fname,
                    fmin=1e-3,
                    fmax=1/60/2,
                    N=100,
                    figsize=(10,7),
                    interpolator_opts={},
                    **kwds):
    """
    Generate two figures, each a 2x2 grid plot (Zxx, Zxy, Zyx, and
    Zyy) showing the magnitude (figure 1) and phase (figure 2) for the
    3-D EMTF given in *xml_fname*. Plot the frequency response for $N$
    uniformly spaced frequencies over the interval from *fmin* to
    *fmax* (in Hz). Use *figsize* sized plots and pass *kwds* as
    keyword arguments to all calls to :func:`plt.plot` and
    :func:`plt.scatter`. Return a tuple with the handles to the two
    generated figures.
    """
    f = np.linspace(fmin, fmax, N)
    f_mHz = f * 1e3
    fmin_mHz = fmin * 1e3
    fmax_mHz = fmax * 1e3
    omega = 2 * np.pi * f
    # setup surface impedance function
    Z_map = parse_xml(xml_fname)
    interp = Zw_interpolator(Z_map, **interpolator_opts)
    # mu_0 * 1e3 converts from [mv / km] / [nT] to [Ohm]
    Zxx_function = lambda omega: interp(omega, 'xx') * mu_0 * 1e3
    Zxy_function = lambda omega: interp(omega, 'xy') * mu_0 * 1e3
    Zyx_function = lambda omega: interp(omega, 'yx') * mu_0 * 1e3
    Zyy_function = lambda omega: interp(omega, 'yy') * mu_0 * 1e3
    # compute transer function response at specified frequencies
    Zxx = Zxx_function(omega)
    Zxy = Zxy_function(omega)
    Zyx = Zyx_function(omega)
    Zyy = Zyy_function(omega)
    # retrieve EMTF values (those provided in the XML file before interpolation)
    Zxx_xml = np.array([x[0,0] for x in Z_map.values()]) * mu_0 * 1e3
    Zxy_xml = np.array([x[0,1] for x in Z_map.values()]) * mu_0 * 1e3
    Zyx_xml = np.array([x[1,0] for x in Z_map.values()]) * mu_0 * 1e3
    Zyy_xml = np.array([x[1,1] for x in Z_map.values()]) * mu_0 * 1e3
    T_xml = np.array(list(Z_map.keys()))
    f_xml = 1 / T_xml
    f_xml_mHz = f_xml * 1e3
    # plot magnitude response in Figure 1
    fig1 = plt.figure(figsize=figsize)
    plt.subplot(221)
    plt.plot(f_mHz,
            np.abs(Zxx) * 1e3,
            label='$|Z_{xx}|$',
            **kwds)
    plt.scatter(f_xml_mHz,
               np.abs(Zxx_xml) * 1e3,
               **kwds)
    plt.xlabel('Frequency (mHz)')
    plt.ylabel('Impedance (m$\Omega$)')
    plt.xlim(fmin_mHz, fmax_mHz)
    plt.legend(loc='lower right')
    plt.subplot(222)
    plt.plot(f_mHz,
            np.abs(Zxy) * 1e3,
            label='$|Z_{xy}|$',
            **kwds)
    plt.scatter(f_xml_mHz,
               np.abs(Zxy_xml) * 1e3,
               **kwds)
    plt.xlabel('Frequency (mHz)')
    plt.ylabel('Impedance (m$\Omega$)')
    plt.xlim(fmin_mHz, fmax_mHz)
    plt.legend(loc='lower right')
    plt.subplot(223)
    plt.plot(f_mHz,
            np.abs(Zyx) * 1e3,
            label='$|Z_{yx}|$',
            **kwds)
    plt.scatter(f_xml_mHz,
               np.abs(Zyx_xml) * 1e3,
               **kwds)
    plt.xlabel('Frequency (mHz)')
    plt.ylabel('Impedance (m$\Omega$)')
    plt.xlim(fmin_mHz, fmax_mHz)
    plt.legend(loc='lower right')
    plt.subplot(224)
    plt.plot(f_mHz,
            np.abs(Zyy) * 1e3,
            label='$|Z_{yy}|$',
            **kwds)
    plt.scatter(f_xml_mHz,
               np.abs(Zyy_xml) * 1e3,
               **kwds)
    plt.xlabel('Frequency (mHz)')
    plt.ylabel('Impedance (m$\Omega$)')
    plt.xlim(fmin_mHz, fmax_mHz)
    plt.legend(loc='lower right')
    plt.suptitle('Magnitude Reponse')
    # plot phase response in Figure 2
    fig2 = plt.figure(figsize=figsize)
    plt.subplot(221)
    plt.plot(f_mHz,
            np.unwrap(np.angle(Zxx, deg=True), discont=180),
            label='$\\angle Z_{xx}$',
            **kwds)
    plt.scatter(f_xml_mHz,
               np.angle(Zxx_xml, deg=True),
               **kwds)
    plt.xlabel('Frequency (mHz)')
    plt.ylabel('Phase (degrees)')
    plt.xlim(fmin_mHz, fmax_mHz)
    plt.legend(loc='upper right')
    plt.subplot(222)
    plt.plot(f_mHz,
            np.unwrap(np.angle(Zxy, deg=True), discont=180),
            label='$\\angle Z_{xy}$',
            **kwds)
    plt.scatter(f_xml_mHz,
               np.angle(Zxy_xml, deg=True),
               **kwds)
    plt.xlabel('Frequency (mHz)')
    plt.ylabel('Phase (degrees)')
    plt.xlim(fmin_mHz, fmax_mHz)
    plt.legend(loc='upper right')
    plt.suptitle('Phase Response')
    plt.subplot(223)
    plt.plot(f_mHz,
            np.unwrap(np.angle(Zyx, deg=True), discont=180),
            label='$\\angle Z_{yx}$',
            **kwds)
    plt.scatter(f_xml_mHz,
               np.angle(Zyx_xml, deg=True),
               **kwds)
    plt.xlabel('Frequency (mHz)')
    plt.ylabel('Phase (degrees)')
    plt.xlim(fmin_mHz, fmax_mHz)
    plt.legend(loc='upper right')
    plt.suptitle('Phase Response')
    plt.subplot(224)
    plt.plot(f_mHz,
            np.angle(Zyy, deg=True),
            label='$\\angle Z_{yy}$',
            **kwds)
    plt.scatter(f_xml_mHz,
               np.unwrap(np.angle(Zyy_xml, deg=True), discont=180),
               **kwds)
    plt.xlabel('Frequency (mHz)')
    plt.ylabel('Phase (degrees)')
    plt.xlim(fmin_mHz, fmax_mHz)
    plt.legend(loc='upper right')
    plt.suptitle('Phase Response')
    plt.show()
    return fig1, fig2
