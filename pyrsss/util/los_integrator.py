from collections.abc import Iterable

import sympy as sym
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

from ..util.chapman import chapman_sym
from ..gnsstk import PyPosition


class SlantIntegrator(object):
    def __init__(self, fun, stn_pos, height1=50, height2=2000):
        """
        Return an object for computing line-of-site integrals of the
        function *f* starting at the point *stn_pos* (a
        :class:`PyPosition`). Furthermore, integration is restricted
        to *height1* [km] <= h <= *height2* [km] where h is the
        geodetic height along the line-of-site, i.e., *f* is assumed
        to be 0 outside this bound. See __call__ for more details
        regarding *fun*.
        """
        self.fun = fun
        self.stn_pos = stn_pos
        self.height1 = height1
        self.height2 = height2

    def __call__(self, sat_pos, args=(), **kwds):
        """
        Return the definite integral of *self.fun*(pos, *args*[0],
        ..., *args*[-1]) for the line of site from *stn_pos* to
        *sat_pos* (a :class:`PyPosition`) where pos is a
        :class:`PyPosition` on the line of site (and note the
        integration bounds on h defined in __init__). The remaining
        *kwds* are passed to the quadrature routine (:py:func:`quad`).
        """
        diff = np.array(sat_pos.xyz) - np.array(self.stn_pos.xyz)
        S_los = np.linalg.norm(diff) / 1e3
        def pos(s):
            """
            Return the ECEF vector a distance *s* along the line-of-site (in
            [km]).
            """
            return PyPosition(*(np.array(self.stn_pos.xyz) + (s / S_los) * diff))
        # determine integration bounds
        # distance along of line of site at which the geodetic height
        # is self.height1
        s1 = minimize_scalar(lambda l: (pos(l).height / 1e3 - self.height1)**2,
                             bounds=[0, S_los],
                             method='Bounded').x
        # distance along of line of site at which the geodetic height
        # is self.height2
        s2 = minimize_scalar(lambda l: (pos(l).height / 1e3 - self.height2)**2,
                             bounds=[0, S_los],
                             method='Bounded').x
        def wrapper(s, *args):
            return self.fun(pos(s), *args)
        return quad(wrapper, s1, s2, args=args, **kwds)[0]


def chapman_sym_scaled(z, Nm, Hm, H_O):
    """
    Return symbolic (i.e., :module:`sympy`) Chapman electron density
    profile function. The returned value is scaled so that integrated
    quantities are in [TECU].
    """
    return chapman_sym(z, Nm, Hm, H_O) / 1e7


class ChapmanSI(SlantIntegrator):
    def __init__(self, stn_pos, **kwds):
        """
        Setup Chapman electron density profile slant integrator.
        """
        z, Nm, Hm, H_O = sym.symbols('z Nm Hm H_O')
        f_sym = chapman_sym_scaled(z, Nm, Hm, H_O)
        f = sym.lambdify((z, Nm, Hm, H_O),
                         f_sym,
                         modules='numexpr')
        wrapper = lambda pos, *args: f(pos.height / 1e3, *args)
        super(ChapmanSI, self).__init__(wrapper, stn_pos, **kwds)

    def __call__(self, sat_pos, Nm, Hm, H_O):
        """
        Return the slant integral of the Chapman electron density profile
        with parameters Nm [cm^-3], Hm [km], and H_O [km].
        """
        return super(ChapmanSI, self).__call__(sat_pos,
                                               args=(Nm, Hm, H_O))


class DNmChapmanSI(SlantIntegrator):
    def __init__(self, stn_pos, **kwds):
        """
        Setup Chapman electron density profile F2 peak density derivative
        slant integrator.
        """
        z, Nm, Hm, H_O = sym.symbols('z Nm Hm H_O')
        f_sym = chapman_sym_scaled(z, Nm, Hm, H_O)
        DNm_sym = sym.diff(f_sym, Nm)
        f = sym.lambdify((z, Hm, H_O),
                         DNm_sym,
                         modules='numexpr')
        wrapper = lambda pos, *args: f(pos.height / 1e3, *args)
        super(DNmChapmanSI, self).__init__(wrapper, stn_pos, **kwds)

    def __call__(self, sat_pos, Nm, Hm, H_O):
        """
        Return the slant integral of the peak density derivative of the
        Chapman electron density profile with parameters Nm [cm^-3],
        Hm [km], and H_O [km].
        """
        return super(DNmChapmanSI, self).__call__(sat_pos,
                                                  args=(Hm, H_O))


class DHmChapmanSI(SlantIntegrator):
    def __init__(self, stn_pos, **kwds):
        """
        Setup Chapman electron density profile F2 peak height derivative
        slant integrator.
        """
        z, Nm, Hm, H_O = sym.symbols('z Nm Hm H_O')
        f_sym = chapman_sym_scaled(z, Nm, Hm, H_O)
        DHm_sym = sym.diff(f_sym, Hm)
        f = sym.lambdify((z, Nm, Hm, H_O),
                         DHm_sym,
                         modules='numexpr')
        wrapper = lambda pos, *args: f(pos.height / 1e3, *args)
        super(DHmChapmanSI, self).__init__(wrapper, stn_pos, **kwds)

    def __call__(self, sat_pos, Nm, Hm, H_O):
        """
        Return the slant integral of the peak height derivative of the
        Chapman electron density profile with parameters Nm [cm^-3],
        Hm [km], and H_O [km].
        """
        return super(DHmChapmanSI, self).__call__(sat_pos,
                                                  args=(Nm, Hm, H_O))


class DH_OChapmanSI(SlantIntegrator):
    def __init__(self, stn_pos, **kwds):
        """
        Setup Chapman electron density profile plasma scale height
        derivative slant integrator.
        """
        z, Nm, Hm, H_O = sym.symbols('z Nm Hm H_O')
        f_sym = chapman_sym_scaled(z, Nm, Hm, H_O)
        DH_O_sym = sym.diff(f_sym, H_O)
        f = sym.lambdify((z, Nm, Hm, H_O),
                         DH_O_sym,
                         modules='numexpr')
        wrapper = lambda pos, *args: f(pos.height / 1e3, *args)
        super(DH_OChapmanSI, self).__init__(wrapper, stn_pos, **kwds)


    def __call__(self, sat_pos, Nm, Hm, H_O):
        """
        Return the slant integral of the plasma scale height derivative of
        the Chapman electron density profile with parameters Nm
        [cm^-3], Hm [km], and H_O [km].
        """
        return super(DH_OChapmanSI, self).__call__(sat_pos,
                                                   args=(Nm, Hm, H_O))
