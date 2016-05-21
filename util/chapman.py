from __future__ import division

from datetime import datetime

import numpy as NP
import scipy.optimize
import sympy as SYM
import pylab as PL


def chapman(z, Nm, Hm, H_O, exp=NP.exp):
    """
    Return Chapman function electron density at height *z*, maximum
    electron density at the F-peak *Nm*, height at the maximum *Hm*,
    and the scale height of atomic oxygen *H_O*. The exponential
    function can be overridden with *exp*.
    """
    return Nm * exp((1 - ((z - Hm) / H_O) - exp(-(z - Hm) / H_O)) / 2)


def chapman_sym(z, Nm, Hm, H_O):
    """
    Return symbolic (i.e., :module:`sympy`) Chapman electron density
    profile function.
    """
    return chapman(z, Nm, Hm, H_O, exp=SYM.exp)


def chapman_fit(alt,
                ne,
                x0=[1e6, 300, 50],
                bounds=[(0, None),
                        (150, 500),
                        (30, 80)],
                verbose=False,
                **kwds):
    """
    """
    # optimization setup
    z, Nm, Hm, H_O = SYM.symbols('z Nm Hm H_O')
    chapman = chapman_sym(z, Nm, Hm, H_O)
    dNm = SYM.diff(chapman, Nm)
    dHm = SYM.diff(chapman, Hm)
    dH_O = SYM.diff(chapman, H_O)
    chapman_f = SYM.lambdify((z, Nm, Hm, H_O),
                             chapman,
                             modules='numexpr')
    dNm_f = SYM.lambdify((z, Nm, Hm, H_O),
                         dNm,
                         modules='numexpr')
    dHm_f = SYM.lambdify((z, Nm, Hm, H_O),
                         dHm,
                         modules='numexpr')
    dH_O_f = SYM.lambdify((z, Nm, Hm, H_O),
                          dH_O,
                          modules='numexpr')
    # define cost function
    y = NP.asarray(ne)
    def J(x):
        Nm, Hm, H_O = x
        if verbose:
            print('-' * 80)
            print(x)
        y_hat = NP.array([chapman_f(z, Nm, Hm, H_O) for z in alt])
        diff = y - y_hat
        J1 = NP.array([dNm_f(z, Nm, Hm, H_O) for z in alt])
        J2 = NP.array([dHm_f(z, Nm, Hm, H_O) for z in alt])
        J3 = NP.array([dH_O_f(z, Nm, Hm, H_O) for z in alt])
        return (NP.dot(diff, diff),
                NP.array([-2 * NP.sum(diff * J1),
                          -2 * NP.sum(diff * J2),
                          -2 * NP.sum(diff * J3)]))
    # minimize cost function
    x_star, f, d = scipy.optimize.fmin_l_bfgs_b(J,
                                                x0,
                                                bounds=bounds,
                                                **kwds)
    return x_star


if __name__ == '__main__':
    from pyglow.pyglow import Point

    N = 200
    alt = NP.linspace(100, 1500, N)

    dt = datetime(2000, 1, 1)
    lat = 0
    lon = 0

    iri_ne = []
    for alt_i in alt:
        point = Point(dt, lat, lon, alt_i)
        point.run_iri()
        iri_ne.append(point.ne)

    Nm_star, Hm_star, H_O_star = chapman_fit(alt, iri_ne, verbose=True)

    chapman_ne = [chapman(z, Nm_star, Hm_star, H_O_star) for z in alt]

    fig = PL.figure(figsize=(6,10))
    PL.plot(iri_ne,
            alt,
            color='b',
            label='IRI')
    PL.plot(chapman_ne,
            alt,
            color='g',
            label='Chapman fit')
    PL.legend()
    PL.xlabel('Electron density [cm$^{-3}$]')
    PL.ylabel('Height [km]')
    PL.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    PL.axis('tight')

    PL.show()
