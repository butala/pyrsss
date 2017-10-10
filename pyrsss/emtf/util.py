import numpy as NP
import pylab as PL
from scipy.constants import mu_0

from pyrsss.emtf.calc_e_3d import parse_xml, Zw_interpolator


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
    keyword arguments to all calls to :func:`PL.plot` and
    :func:`PL.scatter`. Return a tuple with the handles to the two
    generated figures.
    """
    f = NP.linspace(fmin, fmax, N)
    f_mHz = f * 1e3
    fmin_mHz = fmin * 1e3
    fmax_mHz = fmax * 1e3
    omega = 2 * NP.pi * f
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
    Zxx_xml = NP.array([x[0,0] for x in Z_map.values()]) * mu_0 * 1e3
    Zxy_xml = NP.array([x[0,1] for x in Z_map.values()]) * mu_0 * 1e3
    Zyx_xml = NP.array([x[1,0] for x in Z_map.values()]) * mu_0 * 1e3
    Zyy_xml = NP.array([x[1,1] for x in Z_map.values()]) * mu_0 * 1e3
    T_xml = NP.array(Z_map.keys())
    f_xml = 1 / T_xml
    f_xml_mHz = f_xml * 1e3
    # plot magnitude response in Figure 1
    fig1 = PL.figure(figsize=figsize)
    PL.subplot(221)
    PL.plot(f_mHz,
            NP.abs(Zxx) * 1e3,
            label='$|Z_{xx}|$',
            **kwds)
    PL.scatter(f_xml_mHz,
               NP.abs(Zxx_xml) * 1e3,
               **kwds)
    PL.xlabel('Frequency (mHz)')
    PL.ylabel('Impedance (m$\Omega$)')
    PL.xlim(fmin_mHz, fmax_mHz)
    PL.legend(loc='lower right')
    PL.subplot(222)
    PL.plot(f_mHz,
            NP.abs(Zxy) * 1e3,
            label='$|Z_{xy}|$',
            **kwds)
    PL.scatter(f_xml_mHz,
               NP.abs(Zxy_xml) * 1e3,
               **kwds)
    PL.xlabel('Frequency (mHz)')
    PL.ylabel('Impedance (m$\Omega$)')
    PL.xlim(fmin_mHz, fmax_mHz)
    PL.legend(loc='lower right')
    PL.subplot(223)
    PL.plot(f_mHz,
            NP.abs(Zyx) * 1e3,
            label='$|Z_{yx}|$',
            **kwds)
    PL.scatter(f_xml_mHz,
               NP.abs(Zyx_xml) * 1e3,
               **kwds)
    PL.xlabel('Frequency (mHz)')
    PL.ylabel('Impedance (m$\Omega$)')
    PL.xlim(fmin_mHz, fmax_mHz)
    PL.legend(loc='lower right')
    PL.subplot(224)
    PL.plot(f_mHz,
            NP.abs(Zyy) * 1e3,
            label='$|Z_{yy}|$',
            **kwds)
    PL.scatter(f_xml_mHz,
               NP.abs(Zyy_xml) * 1e3,
               **kwds)
    PL.xlabel('Frequency (mHz)')
    PL.ylabel('Impedance (m$\Omega$)')
    PL.xlim(fmin_mHz, fmax_mHz)
    PL.legend(loc='lower right')
    PL.suptitle('Magnitude Reponse')
    # plot phase response in Figure 2
    fig2 = PL.figure(figsize=figsize)
    PL.subplot(221)
    PL.plot(f_mHz,
            NP.unwrap(NP.angle(Zxx, deg=True), discont=180),
            label='$\\angle Z_{xx}$',
            **kwds)
    PL.scatter(f_xml_mHz,
               NP.angle(Zxx_xml, deg=True),
               **kwds)
    PL.xlabel('Frequency (mHz)')
    PL.ylabel('Phase (degrees)')
    PL.xlim(fmin_mHz, fmax_mHz)
    PL.legend(loc='upper right')
    PL.subplot(222)
    PL.plot(f_mHz,
            NP.unwrap(NP.angle(Zxy, deg=True), discont=180),
            label='$\\angle Z_{xy}$',
            **kwds)
    PL.scatter(f_xml_mHz,
               NP.angle(Zxy_xml, deg=True),
               **kwds)
    PL.xlabel('Frequency (mHz)')
    PL.ylabel('Phase (degrees)')
    PL.xlim(fmin_mHz, fmax_mHz)
    PL.legend(loc='upper right')
    PL.suptitle('Phase Response')
    PL.subplot(223)
    PL.plot(f_mHz,
            NP.unwrap(NP.angle(Zyx, deg=True), discont=180),
            label='$\\angle Z_{yx}$',
            **kwds)
    PL.scatter(f_xml_mHz,
               NP.angle(Zyx_xml, deg=True),
               **kwds)
    PL.xlabel('Frequency (mHz)')
    PL.ylabel('Phase (degrees)')
    PL.xlim(fmin_mHz, fmax_mHz)
    PL.legend(loc='upper right')
    PL.suptitle('Phase Response')
    PL.subplot(224)
    PL.plot(f_mHz,
            NP.angle(Zyy, deg=True),
            label='$\\angle Z_{yy}$',
            **kwds)
    PL.scatter(f_xml_mHz,
               NP.unwrap(NP.angle(Zyy_xml, deg=True), discont=180),
               **kwds)
    PL.xlabel('Frequency (mHz)')
    PL.ylabel('Phase (degrees)')
    PL.xlim(fmin_mHz, fmax_mHz)
    PL.legend(loc='upper right')
    PL.suptitle('Phase Response')
    PL.show()
    return fig1, fig2
