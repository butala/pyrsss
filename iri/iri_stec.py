from __future__ import division

from collections import Iterable

import numpy as NP
import scipy.optimize
import scipy.integrate

from pyglow.pyglow import Point

from ..gpstk import PyPosition


ORDER = 80
""" ??? """

def iri_stec(dt, p1, p2, alt1=50, alt2=2000, order=ORDER):
    """
    Integrate IRI electron density for :class:`datetime` *dt* starting
    at :class:`PyPosition` *p1* and ending at *p2*. Further, limit the
    altitude range from *alt1* to *alt2* (given in [km]). Use a fixed
    order quadrature routine with *order*. Return the integrated TEC
    in [TECU].
    """
    xyz1 = NP.array([p1.x, p1.y, p1.z])
    xyz2 = NP.array([p2.x, p2.y, p2.z])
    def los(l):
        xyz = (1 - l) * xyz1 + l * xyz2
        return PyPosition(*xyz)
    l1 = scipy.optimize.minimize_scalar(lambda l: (los(l).height / 1e3 - alt1)**2,
                                        bounds=[0, 1],
                                        method='Bounded').x
    l2 = scipy.optimize.minimize_scalar(lambda l: (los(l).height / 1e3 - alt2)**2,
                                        bounds=[0, 1],
                                        method='Bounded').x
    def iri_tec(l):
        p = los(l)
        lat = p.geodeticLatitude
        lon = p.longitude
        alt = p.height / 1e3
        point = Point(dt, lat, lon, alt)
        point.run_iri()
        ne = point.ne
        if ne < 0:
            return 0
        else:
            return ne
    p1 = los(l1)
    p2 = los(l2)
    delta = p1.distance(p2) / (l2 - l1) * 1e2
    def integrand(l):
        iterable = isinstance(l, Iterable)
        if not iterable:
            l = [l]
        ne = [iri_tec(l_i) * delta for l_i in l]
        if iterable:
            return ne
        else:
            return ne[0]
    stec_cm2, _ = scipy.integrate.fixed_quad(integrand, l1, l2, n=order)
    return stec_cm2 / 1e12


if __name__ == '__main__':
    from datetime import datetime

    dt = datetime(2010, 1, 1)
    stn_xyz = NP.array([4696.986004,    723.992717,   4239.681595]) * 1e3
    sat_xyz = NP.array([10741.320824,  12456.414622,  21019.082339]) * 1e3

    stec = iri_stec(dt,
                    PyPosition(*stn_xyz),
                    PyPosition(*sat_xyz))
    print(stec)
