from __future__ import division

from collections import OrderedDict

import numpy as NP
import scipy.constants
from intervals import FloatInterval


"""
USGS ground conductivity models are available at
ftp://hazards.cr.usgs.gov/Rigler/Conductivity_Latest/
"""


"""
Interior plains - Michigan region.
"""
IP_3 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*IP3 earth conductivity model
*/    150.0,   1950.0,  12900.0,   5000.0,  23000.0,  57000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    100. ,    50.0 ,   8000. ,    200. ,    625. ,   1000. ,    209. ,    50.1 ,    20.0 ,    5.62 ,    1.58 ,    .891 ,/ !Resistivities in Ohm-m
11                             Number of layers from surface

0.0100000                      Conductivity in S/m (layer 1)
1.500e+02                      Layer thickness in m (layer 1)

0.0200000                      Conductivity in S/m (layer 2)
1.950e+03                      Layer thickness in m (layer 2)

0.0001250                      Conductivity in S/m (layer 3)
1.290e+04                      Layer thickness in m (layer 3)

0.0050000                      Conductivity in S/m (layer 4)
5.000e+03                      Layer thickness in m (layer 4)

0.0016000                      Conductivity in S/m (layer 5)
2.300e+04                      Layer thickness in m (layer 5)

0.0010000                      Conductivity in S/m (layer 6)
5.700e+04                      Layer thickness in m (layer 6)

0.0047860                      Conductivity in S/m (layer 7)
1.500e+05                      Layer thickness in m (layer 7)

0.0199520                      Conductivity in S/m (layer 8)
1.600e+05                      Layer thickness in m (layer 8)

0.0501180                      Conductivity in S/m (layer 9)
1.100e+05                      Layer thickness in m (layer 9)

0.1778270                      Conductivity in S/m (layer10)
1.500e+05                      Layer thickness in m (layer10)

0.6309570                      Conductivity in S/m (layer11)
2.300e+05                      Layer thickness in m (layer11)

1.1220100                      Semi-infinite earth conductivity
"""


"""
Piedmont (SE Appalachians)
"""
PT_1 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*PT1 earth conductivity model
*/   6000.0,   7000.0,   1500.0,   3500.0,  21000.0,  29000.0,  32000.0,  45000.0, 105000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/   1000. ,    625. ,    200. ,   1000. ,    800. ,   4000. ,    800. ,   8000. ,    400. ,    50.1 ,    20.0 ,    5.62 ,    1.58 ,    .891 ,/ !Resistivities in Ohm-m
13                             Number of layers from surface

0.0010000                      Conductivity in S/m (layer 1)
6.000e+03                      Layer thickness in m (layer 1)

0.0016000                      Conductivity in S/m (layer 2)
7.000e+03                      Layer thickness in m (layer 2)

0.0050000                      Conductivity in S/m (layer 3)
1.500e+03                      Layer thickness in m (layer 3)

0.0010000                      Conductivity in S/m (layer 4)
3.500e+03                      Layer thickness in m (layer 4)

0.0012500                      Conductivity in S/m (layer 5)
2.100e+04                      Layer thickness in m (layer 5)

0.0002500                      Conductivity in S/m (layer 6)
2.900e+04                      Layer thickness in m (layer 6)

0.0012500                      Conductivity in S/m (layer 7)
3.200e+04                      Layer thickness in m (layer 7)

0.0001250                      Conductivity in S/m (layer 8)
4.500e+04                      Layer thickness in m (layer 8)

0.0025000                      Conductivity in S/m (layer 9)
1.050e+05                      Layer thickness in m (layer 9)

0.0199520                      Conductivity in S/m (layer10)
1.600e+05                      Layer thickness in m (layer10)

0.0501180                      Conductivity in S/m (layer11)
1.100e+05                      Layer thickness in m (layer11)

0.1778270                      Conductivity in S/m (layer12)
1.500e+05                      Layer thickness in m (layer12)

0.6309570                      Conductivity in S/m (layer13)
2.300e+05                      Layer thickness in m (layer13)

1.1220100                      Semi-infinite earth conductivity
"""


"""
Interior Plains (Great Plains)
"""
IP_4 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*IP4 earth conductivity model
*/    100.0,    500.0,   1000.0,    900.0,  16500.0,   9000.0,  19000.0,  53000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    10.0 ,    15.0 ,    1.70 ,    20.0 ,   2801. ,    40.0 ,    265. ,    500. ,    209. ,    50.1 ,    20.0 ,    5.62 ,    1.58 ,    .891 ,/ !Resistivities in Ohm-m
13                             Number of layers from surface

0.1000000                      Conductivity in S/m (layer 1)
1.000e+02                      Layer thickness in m (layer 1)

0.0666000                      Conductivity in S/m (layer 2)
5.000e+02                      Layer thickness in m (layer 2)

0.5882000                      Conductivity in S/m (layer 3)
1.000e+03                      Layer thickness in m (layer 3)

0.0500000                      Conductivity in S/m (layer 4)
9.000e+02                      Layer thickness in m (layer 4)

0.0003570                      Conductivity in S/m (layer 5)
1.650e+04                      Layer thickness in m (layer 5)

0.0250000                      Conductivity in S/m (layer 6)
9.000e+03                      Layer thickness in m (layer 6)

0.0037700                      Conductivity in S/m (layer 7)
1.900e+04                      Layer thickness in m (layer 7)

0.0020000                      Conductivity in S/m (layer 8)
5.300e+04                      Layer thickness in m (layer 8)

0.0047860                      Conductivity in S/m (layer 9)
1.500e+05                      Layer thickness in m (layer 9)

0.0199520                      Conductivity in S/m (layer10)
1.600e+05                      Layer thickness in m (layer10)

0.0501180                      Conductivity in S/m (layer11)
1.100e+05                      Layer thickness in m (layer11)

0.1778270                      Conductivity in S/m (layer12)
1.500e+05                      Layer thickness in m (layer12)

0.6309570                      Conductivity in S/m (layer13)
2.300e+05                      Layer thickness in m (layer13)

1.1220100                      Semi-infinite earth conductivity
"""


"""
Coastal Plain (Georgia)
"""
CP_2 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*CP2 earth conductivity model
*/     15.0,   1485.0,  13500.0,  17000.0,  36000.0,  77000.0, 105000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    50.0 ,    100. ,   2000. ,   5000. ,   2000. ,    800. ,    5.00 ,    50.1 ,    20.0 ,    5.62 ,    1.58 ,    .891 ,/ !Resistivities in Ohm-m
11                             Number of layers from surface

0.0200000                      Conductivity in S/m (layer 1)
1.500e+01                      Layer thickness in m (layer 1)

0.0100000                      Conductivity in S/m (layer 2)
1.485e+03                      Layer thickness in m (layer 2)

0.0005000                      Conductivity in S/m (layer 3)
1.350e+04                      Layer thickness in m (layer 3)

0.0002000                      Conductivity in S/m (layer 4)
1.700e+04                      Layer thickness in m (layer 4)

0.0005000                      Conductivity in S/m (layer 5)
3.600e+04                      Layer thickness in m (layer 5)

0.0012500                      Conductivity in S/m (layer 6)
7.700e+04                      Layer thickness in m (layer 6)

0.2000000                      Conductivity in S/m (layer 7)
1.050e+05                      Layer thickness in m (layer 7)

0.0199520                      Conductivity in S/m (layer 8)
1.600e+05                      Layer thickness in m (layer 8)

0.0501180                      Conductivity in S/m (layer 9)
1.100e+05                      Layer thickness in m (layer 9)

0.1778270                      Conductivity in S/m (layer10)
1.500e+05                      Layer thickness in m (layer10)

0.6309570                      Conductivity in S/m (layer11)
2.300e+05                      Layer thickness in m (layer11)

1.1220100                      Semi-infinite earth conductivity
"""


"""
Mapping between region name and conductivity model.
"""
NAME_MAP = {'CP_2': CP_2,
            'IP_3': IP_3,
            'IP_4': IP_4,
            'PT_1': PT_1}


def parse_conductivity(fid):
    """
    Parse a USGS model conductivity file-like object *fid*. Return an
    ordered mapping between layer depth interval, [upper depth,
    lower_depth) in [m], and the conductivity [(Ohm-m)^{-1}]
    """
    for line in fid:
        if line.startswith('*/'):
            if 'thicknesses' in line.split('!')[1].lower():
                thicknesses = map(float, line[1:].split('/')[1].split(',')[:-1])
            elif 'resistivities' in line.split('!')[1].lower():
                resistivites = map(float, line[1:].split('/')[1].split(',')[:-1])
            else:
                break
    last_depth = 0
    bounds = []
    for depth in thicknesses:
        bounds.append(FloatInterval.closed_open(last_depth, last_depth + depth))
        last_depth += depth
    return OrderedDict(zip(bounds,
                           1 / NP.array(resistivites)))


def surface_impedance_1D(conductivity_map, omega):
    """
    Calculate the surface impedance given the 1-D conductivity model
    *conductivity_map* at angular frequencies *omega*.
    """
    # start at bottom layer
    sigma = conductivity_map.values()[-1]
    # (5)
    k = NP.sqrt(1j * omega * scipy.constants.mu_0 * sigma)
    # (6)
    Z = 1j * omega * scipy.constants.mu_0 / k
    # iterate in reversed order (i.e., interior to exterior) and skip the bottom layer
    # print(Z)
    for interval_i, sigma_i in conductivity_map.items()[-1::-1]:
        k = NP.sqrt(1j * omega * scipy.constants.mu_0 * sigma_i)
        # (7)
        A = k * Z / (1j * omega * scipy.constants.mu_0)
        r = (1 - A) / (1 + A)
        # (8)
        d = interval_i.length
        B = r * NP.exp(-2 * k * d)
        Z = 1j * omega * scipy.constants.mu_0 * (1 - B) / (k * (1 + B))
    return Z


if __name__ == '__main__':
    from StringIO import StringIO

    import pylab as PL

    fid_pt1 = StringIO(PT_1)
    usgs_map_pt1 = parse_conductivity(fid_pt1)

    fid_ip4 = StringIO(IP_4)
    usgs_map_ip4 = parse_conductivity(fid_ip4)


    N = 1000
    depths = NP.logspace(2, 6, N)

    resistivites_pt1 = NP.empty(N)
    for interval, sigma in usgs_map_pt1.iteritems():
        I = [i for i, d in enumerate(depths) if d in interval]
        resistivites_pt1[I] = 1 / sigma

    resistivites_ip4 = NP.empty(N)
    for interval, sigma in usgs_map_ip4.iteritems():
        I = [i for i, d in enumerate(depths) if d in interval]
        resistivites_ip4[I] = 1 / sigma

    fig = PL.figure(figsize=(8.5, 11))
    PL.loglog(resistivites_pt1,
              NP.array(depths) / 1e3)
    PL.gca().invert_yaxis()
    PL.xlabel('Resistivity [$\Omega$ / m]')
    PL.ylabel('Depth [km]')
    PL.title('1-D Resistivity Model for Piedmont (SE Appalachians) Model PT-1')

    fig = PL.figure(figsize=(8.5, 11))
    PL.loglog(resistivites_ip4,
              NP.array(depths) / 1e3)
    PL.gca().invert_yaxis()
    PL.xlabel('Resistivity [$\Omega$ / m]')
    PL.ylabel('Depth [km]')
    PL.title('1-D Resistivity Model for Interior Plains (Great Plains) Model IP-4')

    f = NP.logspace(-5, -1, N)
    omega = 2 * NP.pi * f

    Z_pt1 = surface_impedance_1D(usgs_map_pt1, omega)
    Z_ip4 = surface_impedance_1D(usgs_map_ip4, omega)

    fig = PL.figure(figsize=(11, 8.5))
    PL.semilogx(f,
                NP.abs(Z_pt1),
                label='PT-1',
                c='b')
    PL.semilogx(f,
                NP.abs(Z_ip4),
                label='IP-4',
                c='g')
    PL.grid(which='both')
    PL.xlabel('Frequency [Hz]')
    PL.ylabel('$|Z(\omega)|$ [$\Omega$]')
    PL.legend(loc='upper left')
    PL.title('Frequency response of two layered Earth conductivity models')

    PL.show()
