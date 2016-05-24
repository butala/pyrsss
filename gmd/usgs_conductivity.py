"""
USGS ground conductivity models are available at
ftp://hazards.cr.usgs.gov/Rigler/Conductivity_Latest/
"""


"""
Adirondack Mountains - central core.
"""
AK_1A = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*AK1A earth conductivity model
*/      3.0,   4997.0,  5000.0,  12000.0,  15000.0,  60000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    100. ,  18215. , 10000. ,   1000. ,    25.0 ,    244. ,    159. ,    28.9 ,    7.95 ,    2.40 ,    .891 ,    .479 ,/ !Resistivities in Ohm-m
11                             Number of layers from surface

0.0100000                      Conductivity in S/m (layer 1)
3.000e+00                      Layer thickness in m (layer 1)

0.0000549                      Conductivity in S/m (layer 2)
4.997e+03                      Layer thickness in m (layer 2)

0.0001000                      Conductivity in S/m (layer 3)
5.000e+03                      Layer thickness in m (layer 3)

0.0010000                      Conductivity in S/m (layer 4)
1.200e+04                      Layer thickness in m (layer 4)

0.0400000                      Conductivity in S/m (layer 5)
1.500e+04                      Layer thickness in m (layer 5)

0.0041000                      Conductivity in S/m (layer 6)
6.000e+04                      Layer thickness in m (layer 6)

0.0063000                      Conductivity in S/m (layer 7)
1.500e+05                      Layer thickness in m (layer 7)

0.0346000                      Conductivity in S/m (layer 8)
1.600e+05                      Layer thickness in m (layer 8)

0.1258000                      Conductivity in S/m (layer 9)
1.100e+05                      Layer thickness in m (layer 9)

0.4168000                      Conductivity in S/m (layer10)
1.500e+05                      Layer thickness in m (layer10)

1.1220000                      Conductivity in S/m (layer11)
2.300e+05                      Layer thickness in m (layer11)

2.0892000                      Semi-infinite earth conductivity
"""

"""
Adirondack Mountains - excludes meta-anorthosite.
"""
AK_1B = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*AK1B earth conductivity model
*/      3.0,   9997.0,  12000.0,  15000.0,  60000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    100. ,  10000. ,   1000. ,    25.0 ,    244. ,    159. ,    28.9 ,    7.95 ,    2.40 ,    .891 ,    .479 ,/ !Resistivities in Ohm-m
10                             Number of layers from surface

0.0100000                      Conductivity in S/m (layer 1)
3.000e+00                      Layer thickness in m (layer 1)

0.0001000                      Conductivity in S/m (layer 2)
9.997e+03                      Layer thickness in m (layer 2)

0.0010000                      Conductivity in S/m (layer 3)
1.200e+04                      Layer thickness in m (layer 3)

0.0400000                      Conductivity in S/m (layer 4)
1.500e+04                      Layer thickness in m (layer 4)

0.0041000                      Conductivity in S/m (layer 5)
6.000e+04                      Layer thickness in m (layer 5)

0.0063000                      Conductivity in S/m (layer 6)
1.500e+05                      Layer thickness in m (layer 6)

0.0346000                      Conductivity in S/m (layer 7)
1.600e+05                      Layer thickness in m (layer 7)

0.1258000                      Conductivity in S/m (layer 8)
1.100e+05                      Layer thickness in m (layer 8)

0.4168000                      Conductivity in S/m (layer 9)
1.500e+05                      Layer thickness in m (layer 9)

1.1220000                      Conductivity in S/m (layer10)
2.300e+05                      Layer thickness in m (layer10)

2.0892000                      Semi-infinite earth conductivity
"""

"""
Appalachian Plateaus - southern portion of the Appalachian Plateau.
"""
AP_1 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*AP1 earth conductivity model
*/   4000.0,  12000.0,  25000.0,  14000.0,  45000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    80.0 ,    80.0 ,    20.0 ,    303. ,    100. ,    10.0 ,    50.1 ,    20.0 ,    5.62 ,    1.58 ,    .891 ,/ !Resistivities in Ohm-m
10                             Number of layers from surface

0.0125000                      Conductivity in S/m (layer 1)
4.000e+03                      Layer thickness in m (layer 1)

0.0125000                      Conductivity in S/m (layer 2)
1.200e+04                      Layer thickness in m (layer 2)

0.0500000                      Conductivity in S/m (layer 3)
2.500e+04                      Layer thickness in m (layer 3)

0.0033000                      Conductivity in S/m (layer 4)
1.400e+04                      Layer thickness in m (layer 4)

0.0100000                      Conductivity in S/m (layer 5)
4.500e+04                      Layer thickness in m (layer 5)

0.1000000                      Conductivity in S/m (layer 6)
1.500e+05                      Layer thickness in m (layer 6)

0.0199520                      Conductivity in S/m (layer 7)
1.600e+05                      Layer thickness in m (layer 7)

0.0501180                      Conductivity in S/m (layer 8)
1.800e+05                      Layer thickness in m (layer 8)

0.1778270                      Conductivity in S/m (layer 9)
8.000e+04                      Layer thickness in m (layer 9)

0.6309570                      Conductivity in S/m (layer10)
2.300e+05                      Layer thickness in m (layer10)

1.1220100                      Semi-infinite earth conductivity
"""

"""
Northern Appalachian Plateaus - northern portion of the Appalachian Plateau (southern New York).
"""
AP_2 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*AP2 earth conductivity model
*/     25.0,   2725.0,  12250.0,  10000.0,   9000.0,  66000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    100. ,    303. ,  10000. ,    303. ,   1000. ,    244. ,    159. ,    28.9 ,    7.95 ,    2.40 ,    .891 ,   .479  ,/ !Resistivities in Ohm-m
11                             Number of layers from surface

0.0100000                      Conductivity in S/m (layer 1)
2.500e+01                      Layer thickness in m (layer 1)

0.0033000                      Conductivity in S/m (layer 2)
2.725e+03                      Layer thickness in m (layer 2)

0.0001000                      Conductivity in S/m (layer 3)
1.225e+04                      Layer thickness in m (layer 3)

0.0033000                      Conductivity in S/m (layer 4)
1.000e+04                      Layer thickness in m (layer 4)

0.0010000                      Conductivity in S/m (layer 5)
9.000e+03                      Layer thickness in m (layer 5)

0.0041000                      Conductivity in S/m (layer 6)
6.600e+04                      Layer thickness in m (layer 6)

0.0063000                      Conductivity in S/m (layer 7)
1.500e+05                      Layer thickness in m (layer 7)

0.0346000                      Conductivity in S/m (layer 8)
1.600e+05                      Layer thickness in m (layer 8)

0.1258000                      Conductivity in S/m (layer 9)
1.100e+05                      Layer thickness in m (layer 9)

0.4168000                      Conductivity in S/m (layer10)
1.500e+05                      Layer thickness in m (layer10)

1.1220000                      Conductivity in S/m (layer11)
2.300e+05                      Layer thickness in m (layer11)

2.0892000                      Semi-infinite earth conductivity
"""

"""
Northwest Basin and Range - northwestern margin of the Basin and
Range geophysical province.
"""
BR_1 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*BR1 earth conductivity model
*/     30.0,    970.0,   2000.0,   9000.0,   4000.0,   3000.0,   4000.0,   5000.0,   6000.0,   6000.0,  11000.0,  13000.0,  13000.0,  27000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    10.0 ,    100. ,    10.0 ,    100. ,    114. ,    89.0 ,    40.0 ,    13.3 ,    5.67 ,    4.82 ,    7.62 ,    16.0 ,    32.0 ,    50.0 ,    209. ,    50.1 ,    20.0 ,    5.62 ,    1.58 ,    .891 ,/ !Resistivities in Ohm-m
19                             Number of layers from surface

0.1000000                      Conductivity in S/m (layer 1)
3.000e+01                      Layer thickness in m (layer 1)

0.0100000                      Conductivity in S/m (layer 2)
9.700e+02                      Layer thickness in m (layer 2)

0.1000000                      Conductivity in S/m (layer 3)
2.000e+03                      Layer thickness in m (layer 3)

0.0100000                      Conductivity in S/m (layer 4)
9.000e+03                      Layer thickness in m (layer 4)

0.0087500                      Conductivity in S/m (layer 5)
4.000e+03                      Layer thickness in m (layer 5)

0.0112500                      Conductivity in S/m (layer 6)
3.000e+03                      Layer thickness in m (layer 6)

0.0250000                      Conductivity in S/m (layer 7)
4.000e+03                      Layer thickness in m (layer 7)

0.0750000                      Conductivity in S/m (layer 8)
5.000e+03                      Layer thickness in m (layer 8)

0.1762500                      Conductivity in S/m (layer 9)
6.000e+03                      Layer thickness in m (layer 9)

0.2075000                      Conductivity in S/m (layer10)
6.000e+03                      Layer thickness in m (layer10)

0.1312500                      Conductivity in S/m (layer11)
1.100e+04                      Layer thickness in m (layer11)

0.0625000                      Conductivity in S/m (layer12)
1.300e+04                      Layer thickness in m (layer12)

0.0312500                      Conductivity in S/m (layer13)
1.300e+04                      Layer thickness in m (layer13)

0.0200000                      Conductivity in S/m (layer14)
2.700e+04                      Layer thickness in m (layer14)

0.0047860                      Conductivity in S/m (layer15)
1.500e+05                      Layer thickness in m (layer15)

0.0199520                      Conductivity in S/m (layer16)
1.600e+05                      Layer thickness in m (layer16)

0.0501180                      Conductivity in S/m (layer17)
1.100e+05                      Layer thickness in m (layer17)

0.1778270                      Conductivity in S/m (layer18)
1.500e+05                      Layer thickness in m (layer18)

0.6309570                      Conductivity in S/m (layer19)
2.300e+05                      Layer thickness in m (layer19)

1.1220100                      Semi-infinite earth conductivity
"""

"""
Colorado Plateau - Approximately centered in the Four Corners
region of the southwestern United States, the Colorado Plateau
encompasses southeastern Utah, northern Arizona, northwestern New
Mexico and the western edge of Colorado.
"""
CL_1 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*CL1 earth conductivity model
*/    100.0,   1900.0,  31000.0,  12000.0, 105000.0, 100000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    10.0 ,    50.0 ,   3030. ,    80.0 ,    400. ,    69.9 ,    50.1 ,    20.0 ,    5.62 ,    1.58 ,    .891 ,/ !Resistivities in Ohm-m
10                             Number of layers from surface

0.1000000                      Conductivity in S/m (layer 1)
1.000e+02                      Layer thickness in m (layer 1)

0.0200000                      Conductivity in S/m (layer 2)
1.900e+03                      Layer thickness in m (layer 2)

0.0003300                      Conductivity in S/m (layer 3)
3.100e+04                      Layer thickness in m (layer 3)

0.0125000                      Conductivity in S/m (layer 4)
1.200e+04                      Layer thickness in m (layer 4)

0.0025000                      Conductivity in S/m (layer 5)
1.050e+05                      Layer thickness in m (layer 5)

0.0143000                      Conductivity in S/m (layer 6)
1.000e+05                      Layer thickness in m (layer 6)

0.0199520                      Conductivity in S/m (layer 7)
1.600e+05                      Layer thickness in m (layer 7)

0.0501180                      Conductivity in S/m (layer 8)
1.100e+05                      Layer thickness in m (layer 8)

0.1778270                      Conductivity in S/m (layer 9)
1.500e+05                      Layer thickness in m (layer 9)

0.6309570                      Conductivity in S/m (layer10)
2.300e+05                      Layer thickness in m (layer10)

1.1220100                      Semi-infinite earth conductivity
"""

"""
Columbia Plateau - straddles eastern Washington and Oregon, with
an arcuate extension into southwestern Idaho.
"""
CO_1 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*CO1 earth conductivity model
*/     30.0,   1470.0,   1500.0,   5000.0,   4000.0,   3000.0,   4000.0,   5000.0,   6000.0,   6000.0,  11000.0,  13000.0,  13000.0,  27000.0,  50000.0, 100000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    50.0 ,    100. ,    15.2 ,    303. ,    52.5 ,    53.3 ,    37.2 ,    21.1 ,    15.5 ,    17.0 ,    26.0 ,    42.7 ,    60.4 ,    64.0 ,    152. ,    209. ,    50.1 ,    20.0 ,    5.62  ,   1.58 ,    .891 ,/ !Resistivities in Ohm-m
20                             Number of layers from surface

0.0200000                      Conductivity in S/m (layer 1)
3.000e+01                      Layer thickness in m (layer 1)

0.0100000                      Conductivity in S/m (layer 2)
1.470e+03                      Layer thickness in m (layer 2)

0.0666000                      Conductivity in S/m (layer 3)
1.500e+03                      Layer thickness in m (layer 3)

0.0033000                      Conductivity in S/m (layer 4)
5.000e+03                      Layer thickness in m (layer 4)

0.0190600                      Conductivity in S/m (layer 5)
4.000e+03                      Layer thickness in m (layer 5)

0.0187500                      Conductivity in S/m (layer 6)
3.000e+03                      Layer thickness in m (layer 6)

0.0268700                      Conductivity in S/m (layer 7)
4.000e+03                      Layer thickness in m (layer 7)

0.0475000                      Conductivity in S/m (layer 8)
5.000e+03                      Layer thickness in m (layer 8)

0.0646800                      Conductivity in S/m (layer 9)
6.000e+04                      Layer thickness in m (layer 9)

0.0587500                      Conductivity in S/m (layer10)
6.000e+03                      Layer thickness in m (layer10)

0.0384300                      Conductivity in S/m (layer11)
1.100e+04                      Layer thickness in m (layer11)

0.0234300                      Conductivity in S/m (layer12)
1.300e+04                      Layer thickness in m (layer12)

0.0165600                      Conductivity in S/m (layer13)
1.300e+04                      Layer thickness in m (layer13)

0.0156200                      Conductivity in S/m (layer14)
2.700e+04                      Layer thickness in m (layer14)

0.0066000                      Conductivity in S/m (layer15)
5.000e+04                      Layer thickness in m (layer15)

0.0047860                      Conductivity in S/m (layer16)
1.000e+05                      Layer thickness in m (layer16)

0.0199520                      Conductivity in S/m (layer17)
1.600e+05                      Layer thickness in m (layer17)

0.0501180                      Conductivity in S/m (layer18)
1.100e+05                      Layer thickness in m (layer18)

0.1778270                      Conductivity in S/m (layer19)
1.500e+05                      Layer thickness in m (layer19)

0.6309570                      Conductivity in S/m (layer20)
2.300e+05                      Layer thickness in m (layer20)

1.1220100                      Semi-infinite earth conductivity
"""

"""
Coastal Plain - South Carolina.
"""
CP_1 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*CP1 earth conductivity model
*/      5.0,    395.0,   9600.0,  14000.0,  13000.0,  31000.0,  77000.0, 105000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    400. ,    31.3 ,    250. ,    625. ,   6250. ,   1000. ,    800. ,    5.00 ,    50.1 ,    20.0 ,    5.62 ,    1.58 ,    .891 ,/ !Resistivities in Ohm-m
12                             Number of layers from surface

0.0025000                      Conductivity in S/m (layer 1)
5.000e+00                      Layer thickness in m (layer 1)

0.0320000                      Conductivity in S/m (layer 2)
3.950e+02                      Layer thickness in m (layer 2)

0.0040000                      Conductivity in S/m (layer 3)
9.600e+03                      Layer thickness in m (layer 3)

0.0016000                      Conductivity in S/m (layer 4)
1.400e+04                      Layer thickness in m (layer 4)

0.0001600                      Conductivity in S/m (layer 5)
1.300e+04                      Layer thickness in m (layer 5)

0.0010000                      Conductivity in S/m (layer 6)
3.100e+04                      Layer thickness in m (layer 6)

0.0012500                      Conductivity in S/m (layer 7)
7.700e+04                      Layer thickness in m (layer 7)

0.2000000                      Conductivity in S/m (layer 8)
1.050e+05                      Layer thickness in m (layer 8)

0.0199520                      Conductivity in S/m (layer 9)
1.600e+05                      Layer thickness in m (layer 9)

0.0501180                      Conductivity in S/m (layer10)
1.100e+05                      Layer thickness in m (layer10)

0.1778270                      Conductivity in S/m (layer11)
1.500e+05                      Layer thickness in m (layer11)

0.6309570                      Conductivity in S/m (layer12)
2.300e+05                      Layer thickness in m (layer12)

1.1220100                      Semi-infinite earth conductivity
"""

"""
Coastal Plain - Georgia.
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
Cascade-Sierra Mountains - central Oregon and Washington states,
but extending southward into California.
"""
CS_1 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*CS1 earth conductivity model
*/     30.0,    470.0,   1500.0,  18000.0,  10000.0,  13000.0,  57000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    60.2 ,    200. ,    15.2 ,    182. ,    15.2 ,    152. ,    120. ,    209. ,    50.1 ,    20.0 ,    5.62 ,    1.58 ,    .891 ,/ !Resistivities in Ohm-m
12                             Number of layers from surface

0.0166000                      Conductivity in S/m (layer 1)
3.000e+01                      Layer thickness in m (layer 1)

0.0050000                      Conductivity in S/m (layer 2)
4.700e+02                      Layer thickness in m (layer 2)

0.0660000                      Conductivity in S/m (layer 3)
1.500e+03                      Layer thickness in m (layer 3)

0.0055000                      Conductivity in S/m (layer 4)
1.800e+04                      Layer thickness in m (layer 4)

0.0660000                      Conductivity in S/m (layer 5)
1.000e+04                      Layer thickness in m (layer 5)

0.0066000                      Conductivity in S/m (layer 6)
1.300e+04                      Layer thickness in m (layer 6)

0.0083000                      Conductivity in S/m (layer 7)
5.700e+04                      Layer thickness in m (layer 7)

0.0047860                      Conductivity in S/m (layer 8)
1.500e+05                      Layer thickness in m (layer 8)

0.0199520                      Conductivity in S/m (layer 9)
1.600e+05                      Layer thickness in m (layer 9)

0.0501180                      Conductivity in S/m (layer10)
1.100e+05                      Layer thickness in m (layer10)

0.1778270                      Conductivity in S/m (layer11)
1.500e+05                      Layer thickness in m (layer11)

0.6309570                      Conductivity in S/m (layer12)
2.300e+05                      Layer thickness in m (layer12)

1.1220100                      Semi-infinite earth conductivity
"""

"""
Florida Peninsula - Florida.
"""
FL_1 = """
* Lines starting with * are just comments.
* Text after the numbers is ignored
*FL1 earth conductivity model
*/     1000,   5000,   34000,  60000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    67. ,   224. ,   3162. ,   155. ,   138. ,   34. ,    15.5. ,    4.2 ,    1.2 ,    0.87 ,/ !Resistivities in Ohm-m
9                             Number of layers from surface

0.0149000                      Conductivity in S/m (layer 1)
1.000e+03                      Layer thickness in m (layer 1)

0.0044720		       Conductivity in S/m (layer 2)
5.0000+03		       Layer thickness in m (layer 2)

0.0003160                      Conductivity in S/m (layer 3)
3.400e+04                      Layer thickness in m (layer 3)

0.0064500                      Conductivity in S/m (layer 4)
6.000e+04                      Layer thickness in m (layer 4)

0.0072000                      Conductivity in S/m (layer 5)
1.500e+05                      Layer thickness in m (layer 5)

0.0295000                      Conductivity in S/m (layer 6)
1.600e+05                      Layer thickness in m (layer 6)

0.0644000                      Conductivity in S/m (layer 7)
1.100e+05                      Layer thickness in m (layer 7)

0.2370000                      Conductivity in S/m (layer 8)
1.500e+05                      Layer thickness in m (layer 8)

0.8540000                      Conductivity in S/m (layer 9)
2.300e+05                      Layer thickness in m (layer 9)

1.6880000                      Semi-infinite earth conductivity
"""

"""
Interior Plains - typical of the eastern portion of North Dakota
only.
"""
IP_1 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*IP1 earth conductivity model
*/    330.0,   1670.0,  13000.0,  10000.0,  20000.0,  55000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    30.0 ,    20.0 ,   4545. ,    654.,    3030. ,   5000. ,    209. ,    50.1 ,    20.0 ,    5.62 ,    1.58 ,    .891 ,/ !Resistivities in Ohm-m
11                             Number of layers from surface

0.0333000                      Conductivity in S/m (layer 1)
3.300e+02                      Layer thickness in m (layer 1)

0.0500000                      Conductivity in S/m (layer 2)
1.670e+03                      Layer thickness in m (layer 2)

0.0002200                      Conductivity in S/m (layer 3)
1.300e+04                      Layer thickness in m (layer 3)

0.0015300                      Conductivity in S/m (layer 4)
1.000e+04                      Layer thickness in m (layer 4)

0.0003300                      Conductivity in S/m (layer 5)
2.000e+04                      Layer thickness in m (layer 5)

0.0002000                      Conductivity in S/m (layer 6)
5.500e+04                      Layer thickness in m (layer 6)

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
Interior Plains (North American Conductive Anomaly) - a narrow
region of low resistivity extending across western South and North
Dakota, and north into Canada along the Saskatchewan-Manitoba
provincial boundary.
"""
IP_2 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*IP2 earth conductivity model
*/    330.0,   4570.0,  10100.0,   2000.0,   2000.0,   6000.0,  20000.0,  55000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    30.0 ,    50.0 ,    455. ,    250. ,    1.00 ,    250. ,    250. ,    800. ,    209. ,    50.1 ,    20.0 ,    5.62 ,    1.58 ,    .891 ,/ !Resistivities in Ohm-m
13                             Number of layers from surface

0.0333000                      Conductivity in S/m (layer 1)
3.300e+02                      Layer thickness in m (layer 1)

0.0200000                      Conductivity in S/m (layer 2)
4.570e+03                      Layer thickness in m (layer 2)

0.0022000                      Conductivity in S/m (layer 3)
1.010e+04                      Layer thickness in m (layer 3)

0.0040000                      Conductivity in S/m (layer 4)
2.000e+03                      Layer thickness in m (layer 4)

1.0000000                      Conductivity in S/m (layer 5)
2.000e+03                      Layer thickness in m (layer 5)

0.0040000                      Conductivity in S/m (layer 6)
6.000e+03                      Layer thickness in m (layer 6)

0.0040000                      Conductivity in S/m (layer 7)
2.000e+04                      Layer thickness in m (layer 7)

0.0012500                      Conductivity in S/m (layer 8)
5.500e+04                      Layer thickness in m (layer 8)

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
New England - southwest portion of Maine, just north of Portland.
"""
NE_1 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*NE1 earth conductivity model
*/     30.0,  9970.0,   5000.0,  10000.0,  11000.0,  64000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    60.2 ,  2000. ,   2000. ,    303. ,    303. ,    244. ,    159. ,    28.9 ,    7.95 ,    2.40 ,    .891 ,    .479 ,/ !Resistivities in Ohm-m
11                             Number of layers from surface

0.0166000                      Conductivity in S/m (layer 1)
3.000e+01                      Layer thickness in m (layer 1)

0.0005000                      Conductivity in S/m (layer 2)
9.970e+03                      Layer thickness in m (layer 2)

0.0005000                      Conductivity in S/m (layer 3)
5.000e+03                      Layer thickness in m (layer 3)

0.0033000                      Conductivity in S/m (layer 4)
1.000e+04                      Layer thickness in m (layer 4)

0.0033000                      Conductivity in S/m (layer 5)
1.100e+04                      Layer thickness in m (layer 5)

0.0041000                      Conductivity in S/m (layer 6)
6.400e+04                      Layer thickness in m (layer 6)

0.0063000                      Conductivity in S/m (layer 7)
1.500e+05                      Layer thickness in m (layer 7)

0.0346000                      Conductivity in S/m (layer 8)
1.600e+05                      Layer thickness in m (layer 8)

0.1258000                      Conductivity in S/m (layer 9)
1.100e+05                      Layer thickness in m (layer 9)

0.4168000                      Conductivity in S/m (layer10)
1.500e+05                      Layer thickness in m (layer10)

1.1220100                      Conductivity in S/m (layer11)
2.300e+05                      Layer thickness in m (layer11)

2.0892000                      Semi-infinite earth conductivity
"""

"""
Pacific Border (Willamette Valley) - the Pacific Border
physiographic province, extending southward through Oregon and into
western California.
"""
PB_1 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*PB1 earth conductivity model
*/      3.0,    397.0,   3600.0,   4000.0,  22000.0,   4000.0,  66000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    100. ,    7.50 ,    25.0 ,    85.5 ,    400. ,    30.3 ,    400. ,    209. ,    50.1 ,    20.0 ,    5.62 ,    1.58 ,    .891 ,/ !Resistivities in Ohm-m
12                             Number of layers from surface

0.0100000                      Conductivity in S/m (layer 1)
3.000e+00                      Layer thickness in m (layer 1)

0.1330000                      Conductivity in S/m (layer 2)
3.970e+02                      Layer thickness in m (layer 2)

0.0400000                      Conductivity in S/m (layer 3)
3.600e+03                      Layer thickness in m (layer 3)

0.0117000                      Conductivity in S/m (layer 4)
4.000e+03                      Layer thickness in m (layer 4)

0.0025000                      Conductivity in S/m (layer 5)
2.200e+04                      Layer thickness in m (layer 5)

0.0330000                      Conductivity in S/m (layer 6)
4.000e+03                      Layer thickness in m (layer 6)

0.0025000                      Conductivity in S/m (layer 7)
6.600e+04                      Layer thickness in m (layer 7)

0.0047860                      Conductivity in S/m (layer 8)
1.500e+05                      Layer thickness in m (layer 8)

0.0199520                      Conductivity in S/m (layer 9)
1.600e+05                      Layer thickness in m (layer 9)

0.0501180                      Conductivity in S/m (layer10)
1.100e+05                      Layer thickness in m (layer10)

0.1778270                      Conductivity in S/m (layer11)
1.500e+05                      Layer thickness in m (layer11)

0.6309570                      Conductivity in S/m (layer12)
2.300e+05                      Layer thickness in m (layer12)

1.1220100                      Semi-infinite earth conductivity
"""

"""
Pacific Border (Puget Lowlands) - the Puget Lowlands portion of
the Pacific Border physiographic province.
"""
PB_2 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*PB2 earth conductivity model
*/      3.0,   1097.0,   6900.0,  17000.0,   7000.0,  13000.0,  20000.0,  15000.0,  70000.0, 100000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    100. ,    10.0 ,    20.0 ,    100. ,    30.3 ,    250. ,    400. ,    667. ,    800. ,    209. ,    50.1 ,    20.0 ,    5.62 ,    1.58 ,    .891 ,/ !Resistivities in Ohm-m
14                             Number of layers from surface

0.0100000                      Conductivity in S/m (layer 1)
3.000e+00                      Layer thickness in m (layer 1)

0.1000000                      Conductivity in S/m (layer 2)
1.097e+03                      Layer thickness in m (layer 2)

0.0500000                      Conductivity in S/m (layer 3)
6.900e+03                      Layer thickness in m (layer 3)

0.0100000                      Conductivity in S/m (layer 4)
1.700e+04                      Layer thickness in m (layer 4)

0.0330000                      Conductivity in S/m (layer 5)
7.000e+03                      Layer thickness in m (layer 5)

0.0040000                      Conductivity in S/m (layer 6)
1.300e+04                      Layer thickness in m (layer 6)

0.0025000                      Conductivity in S/m (layer 7)
2.000e+04                      Layer thickness in m (layer 7)

0.0015000                      Conductivity in S/m (layer 8)
1.500e+04                      Layer thickness in m (layer 8)

0.0012500                      Conductivity in S/m (layer 9)
7.000e+04                      Layer thickness in m (layer 9)

0.0047860                      Conductivity in S/m (layer10)
1.000e+05                      Layer thickness in m (layer10)

0.0199520                      Conductivity in S/m (layer11)
1.600e+05                      Layer thickness in m (layer11)

0.0501180                      Conductivity in S/m (layer12)
1.100e+05                      Layer thickness in m (layer12)

0.1778270                      Conductivity in S/m (layer13)
1.500e+05                      Layer thickness in m (layer13)

0.6309570                      Conductivity in S/m (layer14)
2.300e+05                      Layer thickness in m (layer14)

1.1220100                      Semi-infinite earth conductivity
"""

"""
Piedmont (SE Appalachians) - the Piedmont physiographic province,
and is located between CP-1 and AP-1.
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
St. Lawrence Lowlands - upper New York state
"""
SL_1 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*SL1 earth conductivity model
*/      3.0,    175.0,   9825.0,  12000.0,  18000.0, 600000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    4.00 ,    500. ,  10000. ,    667. ,    25.0 ,    244. ,    159. ,    28.9 ,    7.95 ,    2.40 ,    .891 ,    .479 ,/ !Resistivities in Ohm-m
11                             Number of layers from surface

0.2500000                      Conductivity in S/m (layer 1)
3.000e+00                      Layer thickness in m (layer 1)

0.0020000                      Conductivity in S/m (layer 2)
1.720e+02                      Layer thickness in m (layer 2)

0.0001000                      Conductivity in S/m (layer 3)
9.825e+03                      Layer thickness in m (layer 3)

0.0015000                      Conductivity in S/m (layer 4)
1.200e+04                      Layer thickness in m (layer 4)

0.0400000                      Conductivity in S/m (layer 5)
1.800e+04                      Layer thickness in m (layer 5)

0.0041000                      Conductivity in S/m (layer 6)
6.000e+04                      Layer thickness in m (layer 6)

0.0063000                      Conductivity in S/m (layer 7)
1.500e+05                      Layer thickness in m (layer 7)

0.0346000                      Conductivity in S/m (layer 8)
1.600e+05                      Layer thickness in m (layer 8)

0.1258000                      Conductivity in S/m (layer 9)
1.100e+05                      Layer thickness in m (layer 9)

0.4168000                      Conductivity in S/m (layer10)
1.500e+05                      Layer thickness in m (layer10)

1.1220000                      Conductivity in S/m (layer11)
2.300e+05                      Layer thickness in m (layer11)

2.0892000                      Semi-infinite earth conductivity
"""

"""
Superior Upland - the northern portions of Minnesota, Wisconsin,
and Michigan's Upper Peninsula.
"""
SU_1 = """\
* Lines starting with * are just comments.
* Text after the numbers is ignored
*SU1 earth conductivity model
*/     30.0,  13470.0,  11500.0,  13000.0,  62000.0, 150000.0, 160000.0, 110000.0, 150000.0, 230000.0,      INF,/ ! layer thicknesses in m
*/    100. ,   5988. ,    200. ,    303. ,   1000. ,    209. ,    50.1 ,    20.0 ,    5.62 ,    1.58 ,    .891 ,/ !Resistivities in Ohm-m
10                             Number of layers from surface

0.0100000                      Conductivity in S/m (layer 1)
3.000e+01                      Layer thickness in m (layer 1)

0.0001670                      Conductivity in S/m (layer 2)
1.347e+04                      Layer thickness in m (layer 2)

0.0050000                      Conductivity in S/m (layer 3)
1.150e+04                      Layer thickness in m (layer 3)

0.0033000                      Conductivity in S/m (layer 4)
1.300e+04                      Layer thickness in m (layer 4)

0.0010000                      Conductivity in S/m (layer 5)
6.200e+04                      Layer thickness in m (layer 5)

0.0047860                      Conductivity in S/m (layer 6)
1.500e+05                      Layer thickness in m (layer 6)

0.0199520                      Conductivity in S/m (layer 7)
1.600e+05                      Layer thickness in m (layer 7)

0.0501180                      Conductivity in S/m (layer 8)
1.100e+05                      Layer thickness in m (layer 8)

0.1778270                      Conductivity in S/m (layer 9)
1.500e+05                      Layer thickness in m (layer 9)

0.6309570                      Conductivity in S/m (layer10)
2.300e+05                      Layer thickness in m (layer10)

1.1220100                      Semi-infinite earth conductivity
"""


"""
Mapping between USGS region name and conductivity model.
"""
USGS_CONDUCTIVITY_MAP = {
    'AK_1A': AK_1A,
    'AK_1B': AK_1B,
    'AP_1': AP_1,
    'AP_2': AP_2,
    'BR_1': BR_1,
    'CL_1': CL_1,
    'CO_1': CO_1,
    'CP_1': CP_1,
    'CP_2': CP_2,
    'CS_1': CS_1,
    'FL_1': FL_1,
    'IP_1': IP_1,
    'IP_2': IP_2,
    'IP_3': IP_3,
    'IP_4': IP_4,
    'NE_1': NE_1,
    'PB_1': PB_1,
    'PB_2': PB_2,
    'PT_1': PT_1,
    'SL_1': SL_1,
    'SU_1': SU_1}
