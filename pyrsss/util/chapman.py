from datetime import datetime

import numpy as np
import scipy.optimize
import sympy as sym
import matplotlib.pyplot as plt


def chapman(z, Nm, Hm, H_O, exp=np.exp):
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
    return chapman(z, Nm, Hm, H_O, exp=sym.exp)


def chapman_vec(z_vec, Nm_vec, Hm_vec, H_O_vec):
    """
    Vectorized implementation of the Chapman function evaluation
    routine :func:`chapman`. The input arguments must be sequences
    with the same length and the output is an :class:`np.ndarray` with
    that length.
    """
    try:
        chapman_vec._chapman_sym_f
    except AttributeError:
        sym_vars = sym.symbols('z Nm Hm H_O')
        chapman_vec._chapman_sym_f = sym.lambdify(sym_vars,
                                                  chapman_sym(*sym_vars),
                                                  modules='numexpr')
    return chapman_vec._chapman_sym_f(z_vec,
                                      Nm_vec,
                                      Hm_vec,
                                      H_O_vec)


def chapman_fit(alt,
                ne,
                x0=[1e6, 300, 50],
                bounds=[(1, None),
                        (150, 500),
                        (30, 80)],
                verbose=False,
                **kwds):
    """
    """
    # optimization setup
    z, Nm, Hm, H_O = sym.symbols('z Nm Hm H_O')
    chapman = chapman_sym(z, Nm, Hm, H_O)
    dNm = sym.diff(chapman, Nm)
    dHm = sym.diff(chapman, Hm)
    dH_O = sym.diff(chapman, H_O)
    chapman_f = sym.lambdify((z, Nm, Hm, H_O),
                             chapman,
                             modules='numexpr')
    dNm_f = sym.lambdify((z, Nm, Hm, H_O),
                         dNm,
                         modules='numexpr')
    dHm_f = sym.lambdify((z, Nm, Hm, H_O),
                         dHm,
                         modules='numexpr')
    dH_O_f = sym.lambdify((z, Nm, Hm, H_O),
                          dH_O,
                          modules='numexpr')
    # define cost function
    y = np.asarray(ne)
    def J(x):
        Nm, Hm, H_O = x
        if verbose:
            print('-' * 80)
            print(x)
        y_hat = np.array([chapman_f(z, Nm, Hm, H_O) for z in alt])
        diff = y - y_hat
        J1 = np.array([dNm_f(z, Nm, Hm, H_O) for z in alt])
        J2 = np.array([dHm_f(z, Nm, Hm, H_O) for z in alt])
        J3 = np.array([dH_O_f(z, Nm, Hm, H_O) for z in alt])
        return (np.dot(diff, diff),
                np.array([-2 * np.sum(diff * J1),
                          -2 * np.sum(diff * J2),
                          -2 * np.sum(diff * J3)]))
    # minimize cost function
    x_star, f, d = scipy.optimize.fmin_l_bfgs_b(J,
                                                x0,
                                                bounds=bounds,
                                                **kwds)
    assert d['warnflag'] == 0
    return x_star


if __name__ == '__main__':
    from pyglow.pyglow import Point

    N = 200
    alt = np.linspace(100, 1500, N)

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

    fig = plt.figure(figsize=(6,10))
    plt.plot(iri_ne,
            alt,
            color='b',
            label='IRI')
    plt.plot(chapman_ne,
            alt,
            color='g',
            label='Chapman fit')
    plt.legend()
    plt.xlabel('Electron density [cm$^{-3}$]')
    plt.ylabel('Height [km]')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.axis('tight')

    plt.show()
