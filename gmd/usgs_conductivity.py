from collections import OrderedDict

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


def parse_conductivity(fid):
    """ ??? """
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
                           resistivites))


if __name__ == '__main__':
    from StringIO import StringIO

    import numpy as NP
    import pylab as PL

    fid = StringIO(IP_3)
    usgs_map = parse_conductivity(fid)

    N = 1000
    depths = NP.logspace(2, 6, N)
    resistivites = NP.empty(N)

    for interval, r in usgs_map.iteritems():
        I = [i for i, d in enumerate(depths) if d in interval]
        resistivites[I] = r

    fig = PL.figure(figsize=(8.5, 11))
    PL.loglog(resistivites,
              NP.array(depths) / 1e3)
    PL.gca().invert_yaxis()
    PL.xlabel('Resistivity [$\Omega$ / m]')
    PL.ylabel('Depth [km]')
    PL.title('1-D Resistivity Model for Central Lowland-East Lake Section Model IP-3')

    PL.show()
