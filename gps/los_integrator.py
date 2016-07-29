from collections import Iterable

import sympy as SYM
import numpy as NP
from scipy.integrate import fixed_quad
from scipy.optimize import minimize_scalar

from ..util.chapman import chapman_sym
from ..gpstk import PyPosition


class SlantIntegrator(object):
    def __init__(self, f, stn_xyz, height1=50, height2=2000):
        """
        Return an object for computing line-of-site integrals of the
        function *f* starting at the point *stn_xyz* (3 element XYZ
        ECEF coordinates in [m]). Furthermore, integration is
        restricted to *height1* [km] <= h <= *height2* [km] where h is
        the geodetic height along the line-of-site, i.e., *f* is
        assumed to be 0 outside this bound. See __call__ for more
        details regarding *f*.
        """
        self.f = f
        self.stn_xyz = NP.array(stn_xyz)
        self.height1 = height1
        self.height2 = height2

    def __call__(self, sat_xyz, args=(), n=50, **kwds):
        """
        Return the definite integral of *self.f*(h, *args*[0], ...,
        *args*[-1]) for the line of site from *stn_xyz* to *sat_xyz*
        (3 element XYZ ECEF coordinates in [m]) where h is the
        geodetic height [km] (and note the integration bounds on h
        defined in __init__). The order of the fixed quadrature
        evaluation is *n* (and is the trade-off between accuracy and
        computation). The remaining *kwds* are passed to the
        quadrature routine (:py:func:`fixed_quad`).
        """
        sat_xyz = NP.array(sat_xyz)
        S_los = NP.linalg.norm(sat_xyz - self.stn_xyz) / 1e3
        def los(s):
            """
            Return the ECEF vector a distance *s* along the line-of-site (in
            [km]).
            """
            return self.stn_xyz + (s / S_los) * (sat_xyz - self.stn_xyz)
        # determine integration bounds
        def los_height(s):
            """
            Return the geodetic height for the vector at a distance *s* along
            the line-of-site (in [km]). If the parameter *s* is
            iterable, return a list of heights.
            """
            if isinstance(s, Iterable):
                return [PyPosition(*los(s_i)).height / 1e3 for s_i in s]
            else:
                return PyPosition(*los(s)).height / 1e3
        # distance along of line of site at which the geodetic height
        # is self.height1
        s1 = minimize_scalar(lambda l: (los_height(l) - self.height1)**2,
                             bounds=[0, S_los],
                             method='Bounded').x
        # distance along of line of site at which the geodetic height
        # is self.height2
        s2 = minimize_scalar(lambda l: (los_height(l) - self.height2)**2,
                             bounds=[0, S_los],
                             method='Bounded').x
        def wrapper(s):
            wrapper_args = [los_height(s)] + list(args)
            return self.f(*wrapper_args)
        # scipy.integrate.quad was often fragile, reporting
        # integration warning for reasons that were not fully
        # diagnosed
        return fixed_quad(wrapper, s1, s2, n=n, **kwds)[0]


def chapman_sym_scaled(z, Nm, Hm, H_O):
    """
    Return symbolic (i.e., :module:`sympy`) Chapman electron density
    profile function. The returned value is scaled so that integrated
    quantities are in [TECU].
    """
    return chapman_sym(z, Nm, Hm, H_O) / 1e7


class ChapmanSI(SlantIntegrator):
    def __init__(self, stn_xyz, **kwds):
        """
        Setup Chapman electron density profile slant integrator.
        """
        z, Nm, Hm, H_O = SYM.symbols('z Nm Hm H_O')
        f_sym = chapman_sym_scaled(z, Nm, Hm, H_O)
        f = SYM.lambdify((z, Nm, Hm, H_O),
                         f_sym,
                         modules='numexpr')
        super(ChapmanSI, self).__init__(f, stn_xyz, **kwds)

    def __call__(self, sat_xyz, Nm, Hm, H_O):
        """
        Return the slant integral of the Chapman electron density profile
        with parameters Nm [cm^-3], Hm [km], and H_O [km].
        """
        return super(ChapmanSI, self).__call__(sat_xyz,
                                               (Nm, Hm, H_O))


class DNmChapmanSI(SlantIntegrator):
    def __init__(self, stn_xyz, **kwds):
        """
        Setup Chapman electron density profile F2 peak density derivative
        slant integrator.
        """
        z, Nm, Hm, H_O = SYM.symbols('z Nm Hm H_O')
        f_sym = chapman_sym_scaled(z, Nm, Hm, H_O)
        DNm_sym = SYM.diff(f_sym, Nm)
        f = SYM.lambdify((z, Hm, H_O),
                         DNm_sym,
                         modules='numexpr')
        super(DNmChapmanSI, self).__init__(f, stn_xyz, **kwds)

    def __call__(self, sat_xyz, Nm, Hm, H_O):
        """
        Return the slant integral of the peak density derivative of the
        Chapman electron density profile with parameters Nm [cm^-3],
        Hm [km], and H_O [km].
        """
        return super(DNmChapmanSI, self).__call__(sat_xyz,
                                                  (Hm, H_O))

class DHmChapmanSI(SlantIntegrator):
    def __init__(self, stn_xyz, **kwds):
        """
        Setup Chapman electron density profile F2 peak height derivative
        slant integrator.
        """
        z, Nm, Hm, H_O = SYM.symbols('z Nm Hm H_O')
        f_sym = chapman_sym_scaled(z, Nm, Hm, H_O)
        DHm_sym = SYM.diff(f_sym, Hm)
        f = SYM.lambdify((z, Nm, Hm, H_O),
                         DHm_sym,
                         modules='numexpr')
        super(DHmChapmanSI, self).__init__(f, stn_xyz, **kwds)

    def __call__(self, sat_xyz, Nm, Hm, H_O):
        """
        Return the slant integral of the peak height derivative of the
        Chapman electron density profile with parameters Nm [cm^-3],
        Hm [km], and H_O [km].
        """
        return super(DHmChapmanSI, self).__call__(sat_xyz,
                                                  (Nm, Hm, H_O))

class DH_OChapmanSI(SlantIntegrator):
    def __init__(self, stn_xyz, **kwds):
        """
        Setup Chapman electron density profile plasma scale height
        derivative slant integrator.
        """
        z, Nm, Hm, H_O = SYM.symbols('z Nm Hm H_O')
        f_sym = chapman_sym_scaled(z, Nm, Hm, H_O)
        DH_O_sym = SYM.diff(f_sym, H_O)
        f = SYM.lambdify((z, Nm, Hm, H_O),
                         DH_O_sym,
                         modules='numexpr')
        super(DH_OChapmanSI, self).__init__(f, stn_xyz, **kwds)


    def __call__(self, sat_xyz, Nm, Hm, H_O):
        """
        Return the slant integral of the plasma scale height derivative of
        the Chapman electron density profile with parameters Nm
        [cm^-3], Hm [km], and H_O [km].
        """
        return super(DH_OChapmanSI, self).__call__(sat_xyz,
                                                   (Nm, Hm, H_O))
