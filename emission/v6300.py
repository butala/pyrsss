from __future__ import division

import math
from datetime import datetime

from pyglow.pyglow import Point
# using back-port for python2 (enum34 package) --- this is in the
# python3 standard library (>=3.4)
from enum import Enum


OplusType = Enum('O2type', 'ne charge_neutrality')


"""
The following (alpha1 through BETA_1D) are from Table 2.2 from
Makela's dissertation. They are in turn from:

Link and Cogger, A reexamination of the O I 6300 \AA nightglow,
Journal of Geophysical Research, vol. 93, pp. 9883--9892, 1988.
"""

def alpha1(Te):
    """[cm^3 / s]"""
    return 1.95e-7 * (Te / 300.)**(-0.7)


def alpha2(Te):
    """[cm^3 / s]"""
    return 4.00e-7 * (Te / 300.)**(-0.9)


def k1(Ti):
    """[cm^3 / s]"""
    return 3.23e-12 * math.exp(3.72/(Ti/300) - 1.87/(Ti/300)**2)


def k2(Ti):
    """[cm^3 / s]"""
    return 2.78e-13 * math.exp(2.07/(Ti/300) - 0.61/(Ti/300)**2)


def k3(Tn):
    """[cm^3 / s]"""
    return 2.0e-11 * math.exp(111.8/Tn)


def k4(Tn):
    """[cm^3 / s]"""
    return 2.9e-11 * math.exp(67.5/Tn)


def k5(Tn):
    """[cm^3 / s]"""
    return 1.6e-12 * Tn**(0.91)


A_1D = 7.45e-3
"""[1/s]"""


A_6300 = 5.63e-3
"""[1/s]"""

BETA_1D = 1.1




def Oplus_simple(ne):
    """
    """
    return ne


def Oplus(ne,
          Te,
          Ti,
          O2,
          N2):
    """
    """
    return ne / (1 \
                 + k1(Ti) * O2 / (alpha1(Te) * ne) \
                 + k2(Ti) * N2 / (alpha2(Te) * ne))

def emission_v6300(ne,
                   Te,
                   Ti,
                   Tn,
                   O2,
                   N2,
                   oplus_type=OplusType.charge_neutrality):
    """
    """
    if oplus_type == OplusType.ne:
        oplus = Oplus_simple(ne)
    elif oplus_type == OplusType.charge_neutrality:
        oplus = Oplus(ne, Te, Ti, O2, N2)
    else:
        raise NotImplemented('oplus_type = ' + str(oplus_type))
    N = (A_1D / A_6300) * BETA_1D * k1(Ti) * O2 * oplus
    D = 1 + (k3(Tn) * N2 + k4(Tn) * O2 + k5(Tn) * ne) / A_1D
    return N / D
