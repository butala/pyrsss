from __future__ import division

from datetime import datetime
# from math import exp

import sympy as SP
from sympy import exp
import numpy as NP
import pylab as PL
from pyglow.pyglow import Point
# using back-port for python2 (enum34 package) --- this is in the
# python3 standard library (>=3.4)
from enum import Enum


OplusType = Enum('O2type', 'ne charge_neutrality')



class V6300(object):
    @staticmethod
    def alpha1(Te):
        return 1.95e-7 * (Te / 300)**(-0.7)

    @staticmethod
    def alpha2(Te):
        return 4.00e-7 * (Te / 300)**(-0.9)

    @staticmethod
    def k1(Ti):
        return 3.23e-12 * exp(3.72/(Ti/300) - 1.87/(Ti/300)**2)

    @staticmethod
    def k2(Ti):
        return 2.78e-13 * exp(2.07/(Ti/300) - 0.61/(Ti/300)**2)

    @staticmethod
    def k3(Tn):
        return 2.0e-11 * exp(111.8/Tn)

    @staticmethod
    def k4(Tn):
        return 2.9e-11 * exp(67.5/Tn)

    @staticmethod
    def k5(Tn):
        return 1.6e-12 * Tn**(0.91)

    A_1D = 6.81e-3
    beta_1D = 1.1

    @staticmethod
    def Oplus_simple(ne):
        return ne

    @staticmethod
    def Oplus(ne,
              Te,
              Ti,
              O2,
              N2):
        return ne / (1 \
                     + V6300.k1(Ti) * O2 / (V6300.alpha1(Te) * ne) \
                     + V6300.k2(Ti) * N2 / (V6300.alpha2(Te) * ne))

    @staticmethod
    def emission(ne,
                 Te,
                 Ti,
                 Tn,
                 O2,
                 N2,
                 oplus_type=OplusType.charge_neutrality):
        if oplus_type == OplusType.ne:
            oplus = V6300.Oplus_simple(ne)
        elif oplus_type == OplusType.charge_neutrality:
            oplus = V6300.Oplus(ne, Te, Ti, O2, N2)
        else:
            raise NotImplemented('oplus_type = ' + str(oplus_type))
        N = 0.76 * V6300.beta_1D * V6300.k1(Ti) * O2 * oplus
        D = 1 + (V6300.k3(Tn) * N2 + V6300.k4(Tn) * O2 + V6300.k5(Tn) * ne) / V6300.A_1D
        return N / D


if __name__ == '__main__':
    ne, Te, Ti, Tn, O2, N2 = SP.symbols('ne Te Ti Tn O2 N2')
    v6300 = V6300.emission(ne, Te, Ti, Tn, O2, N2)
    v6300_f = SP.lambdify((ne, Te, Ti, Tn, O2, N2),
                          v6300)

    d_v6300_ne = SP.simplify(SP.diff(v6300, ne))
    d_v6300_ne_f = SP.lambdify((ne, Te, Ti, Tn, O2, N2),
                               d_v6300_ne)

    dn = datetime(1999, 3, 13)
    lat = 25
    lon = 0
    alts = NP.linspace(100, 1000, 256)

    v6300_ne = []
    v6300_neutrality = []
    d_v6300_neutrality = []
    for alt in alts:
        pt = Point(dn, lat, lon, alt)
        pt.run_msis()
        pt.run_iri()
        v6300_ne.append(V6300.emission(pt.ne,
                                       pt.Te,
                                       pt.Ti,
                                       pt.Tn_msis,
                                       pt.nn['O2'],
                                       pt.nn['N2'],
                                       oplus_type=OplusType.ne))
        # v6300_neutrality.append(V6300.emission(pt.ne,
        #                                        pt.Te,
        #                                        pt.Ti,
        #                                        pt.Tn_msis,
        #                                        pt.nn['O2'],
        #                                        pt.nn['N2'],
        #                                        oplus_type=OplusType.charge_neutrality))
        v6300_neutrality.append(v6300_f(pt.ne,
                                        pt.Te,
                                        pt.Ti,
                                        pt.Tn_msis,
                                        pt.nn['O2'],
                                        pt.nn['N2']))
        d_v6300_neutrality.append(d_v6300_ne_f(pt.ne,
                                               pt.Te,
                                               pt.Ti,
                                               pt.Tn_msis,
                                               pt.nn['O2'],
                                               pt.nn['N2']))


    PL.close('all')

    fig = PL.figure(figsize=(8.5, 11))
    PL.plot(v6300_ne,
            alts,
            color='b',
            label='[O$^+$] = [e]')
    PL.plot(v6300_neutrality,
            alts,
            color='r',
            label='[e] = [O$^+$] + [$\mathrm{O}_2^+$] + [NO$^+$]')
    PL.xlabel('V$_{6300}$ [ph/cm$^3$/s]')
    PL.ylabel('Altitude [km]')
    PL.legend()
    PL.title('630.0 [nm] Emission vs. Altitude')

    fig = PL.figure(figsize=(8.5, 11))
    PL.plot(d_v6300_neutrality,
            alts,
            color='b')

    PL.show()
