from collections import OrderedDict

import numpy as np
import scipy.constants
from intervals import FloatInterval


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
                resistivites = list(map(float, line[1:].split('/')[1].split(',')[:-1]))
            else:
                break
    last_depth = 0
    bounds = []
    for depth in thicknesses:
        bounds.append(FloatInterval.closed_open(last_depth, last_depth + depth))
        last_depth += depth
    return OrderedDict(zip(bounds,
                           1 / np.array(resistivites)))


def surface_impedance_1D(conductivity_map, omega):
    """
    Calculate the surface impedance [Ohm] given the 1-D conductivity
    model *conductivity_map* at angular frequencies *omega* [rad].
    """
    # check that bottom layer is an open half space
    assert list(conductivity_map.keys())[-1].length == float('inf')
    # start at bottom layer
    sigma = list(conductivity_map.values())[-1]
    # (5) in NERC, Application Guide: Computing
    # geomagnetically-induced current in the bulk power-system, 2013.
    #
    # http://www.nerc.com/pa/Stand/Project201303GeomagneticDisturbanceMitigation/GIC%20Application%20Guide%202013_approved.pdf
    k = np.sqrt(1j * omega * scipy.constants.mu_0 * sigma)
    # (6)
    Z = 1j * omega * scipy.constants.mu_0 / k
    # iterate in reversed order (i.e., interior to exterior) and skip the bottom layer
    for interval_i, sigma_i in list(conductivity_map.items())[-2::-1]:
        k = np.sqrt(1j * omega * scipy.constants.mu_0 * sigma_i)
        # (7)
        A = k * Z / (1j * omega * scipy.constants.mu_0)
        r = (1 - A) / (1 + A)
        # (8)
        d = interval_i.length
        B = r * np.exp(-2 * k * d)
        Z = 1j * omega * scipy.constants.mu_0 * (1 - B) / (k * (1 + B))
    return Z


if __name__ == '__main__':
    from .usgs_conductivity import USGS_MODEL_MAP

    import matplotlib.pyplot as plt

    usgs_map_pt1 = USGS_MODEL_MAP['PT_1']
    usgs_map_ip4 = USGS_MODEL_MAP['IP_4']

    N = 1000
    depths = np.logspace(2, 6, N)

    resistivites_pt1 = np.empty(N)
    for interval, sigma in usgs_map_pt1.items():
        I = [i for i, d in enumerate(depths) if d in interval]
        resistivites_pt1[I] = 1 / sigma

    resistivites_ip4 = np.empty(N)
    for interval, sigma in usgs_map_ip4.items():
        I = [i for i, d in enumerate(depths) if d in interval]
        resistivites_ip4[I] = 1 / sigma

    fig = plt.figure(figsize=(8.5, 11))
    plt.loglog(resistivites_pt1,
              np.array(depths) / 1e3)
    plt.gca().invert_yaxis()
    plt.xlabel('Resistivity [$\Omega$ / m]')
    plt.ylabel('Depth [km]')
    plt.title('1-D Resistivity Model for Piedmont (SE Appalachians) Model PT-1')

    fig = plt.figure(figsize=(8.5, 11))
    plt.loglog(resistivites_ip4,
              np.array(depths) / 1e3)
    plt.gca().invert_yaxis()
    plt.xlabel('Resistivity [$\Omega$ / m]')
    plt.ylabel('Depth [km]')
    plt.title('1-D Resistivity Model for Interior Plains (Great Plains) Model IP-4')

    f = np.logspace(-5, -1, N)
    omega = 2 * np.pi * f

    Z_pt1 = surface_impedance_1D(usgs_map_pt1, omega)
    Z_ip4 = surface_impedance_1D(usgs_map_ip4, omega)

    fig = plt.figure(figsize=(11, 8.5))
    plt.semilogx(f,
                np.abs(Z_pt1),
                label='PT-1',
                c='b')
    plt.semilogx(f,
                np.abs(Z_ip4),
                label='IP-4',
                c='g')
    plt.grid(which='both')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('$|Z(\omega)|$ [Ohm]')
    plt.legend(loc='upper left')
    plt.title('Frequency response of two layered Earth conductivity models')

    plt.show()
