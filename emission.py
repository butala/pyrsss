from datetime import datetime
from math import exp

import numpy as NP
import pylab as PL
from pyglow.pyglow import Point
# using back-port for python2 (enum34 package) --- this is in the
# python3 standard library (>=3.4)
from enum import Enum


OplusType = Enum('O2type', 'ne charge_neutrality')


def Oplus(pt, oplus=OplusType.charge_neutrality):
    if oplus == OplusType.ne:
        return pt.ne
    elif oplus == OplusType.charge_neutrality:
        alpha1 = 1.95e-7 * (pt.Te / 300)**(-0.7)
        alpha2 = 4.00e-7 * (pt.Te / 300)**(-0.9)
        k1 = 3.23e-12 * exp(3.72/(pt.Ti/300) - 1.87/(pt.Ti/300)**2)
        k2 = 2.78e-13 * exp(2.07/(pt.Ti/300) - 0.61/(pt.Ti/300)**2)
        return pt.ne / (1 \
                        + k1*pt.nn['O2']/(alpha1*pt.ne) \
                        + k2*pt.nn['N2']/(alpha2*pt.ne))
    else:
        raise NotImplemented('Unknown OplusType')


def v6300(pt, oplus=OplusType.charge_neutrality):
    '''Using approximation [O+] = [e]'''
    pt.run_msis()
    pt.run_iri()
    alpha1 = 1.95e-7 * (pt.Te / 300)**(-0.7)
    alpha2 = 4.00e-7 * (pt.Te / 300)**(-0.9)
    beta1D = 1.1
    k1 = 3.23e-12 * exp(3.72/(pt.Ti/300) - 1.87/(pt.Ti/300)**2)
    k2 = 2.78e-13 * exp(2.07/(pt.Ti/300) - 0.61/(pt.Ti/300)**2)
    k3 = 2.0e-11 * exp(111.8/pt.Tn_msis)
    k4 = 2.9e-11 * exp(67.5/pt.Tn_msis)
    k5 = 1.6e-12 * pt.Tn_msis**(0.91)
    A1D = 6.81e-3
    N = 0.76 * beta1D * k1 * pt.nn['O2'] * Oplus(pt, oplus=oplus)
    D = 1 + (k3 * pt.nn['N2'] + k4 * pt.nn['O2'] + k5 * pt.ne) / A1D
    return N / D


if __name__ == '__main__':
    dn = datetime(1999, 3, 13)
    lat = 25
    lon = 0
    alts = NP.linspace(100, 1000, 256)

    v6300_ne = []
    v6300_neutrality = []
    for alt in alts:
        pt = Point(dn, lat, lon, alt)
        v6300_ne.append(v6300(pt, oplus=OplusType.ne))
        v6300_neutrality.append(v6300(pt, oplus=OplusType.charge_neutrality))

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

    PL.show()
