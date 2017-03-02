from collections import OrderedDict, namedtuple


"""
From http://www.carisma.ca/station-information
"""
STATION_INFO = \
"""\
ANNA	Ann Arbor	42.417	276.098	52.88	349.51	2.79
BACK	Back Lake	57.707	265.794	67.32	333.52	6.83
CONT	Contwoyto	65.754	248.750	72.82	304.82	11.64
DAWS	Dawson City	64.048	220.890	65.90	273.89	6.09
ESKI	Eskimo Point	61.106	265.950	70.52	333.15	9.13
FCHP	Fort Chipewyan	58.769	248.894	66.23	308.55	6.25
FCHU	Fort Churchill	58.763	265.920	68.32	333.54	7.44
FSIM	Fort Simpson	61.756	238.770	67.23	294.29	6.78
FSMI	Fort Smith	60.017	248.050	67.28	306.90	6.81
GILL	Gillam	56.376	265.360	66.03	333.05	6.15
GULL	Gull Lake	50.061	251.739	58.20	314.92	3.66
ISLL	Island Lake	53.856	265.340	63.62	333.36	5.15
LGRR	Little Grand Rapids	52.035	264.537	61.81	332.38	4.55
MCMU	Fort McMurray	56.657	248.790	64.17	309.20	5.35
MSTK	Ministik Lake	53.351	247.026	60.61	307.99	4.22
NORM	Norman Wells	65.257	233.311	69.53	285.74	8.31
OSAK	Osakis	45.871	264.917	55.81	333.20	3.22
OXFO	Oxford House	54.929	264.713	64.60	332.28	5.52
PINA	Pinawa	50.199	263.960	59.98	331.75	4.06
POLS	Polson	47.664	245.791	54.71	308.04	3.04
RABB	Rabbit Lake	58.222	256.320	66.85	319.11	6.57
RANK	Rankin Inlet	62.824	267.890	72.22	335.97	10.89
SACH	Sachs Harbour	71.980	234.770	76.17	280.31	N/A
TALO	Taloyoak	69.540	266.450	78.28	330.93	N/A
THRF	Thief River Falls	48.027	263.635	57.82	331.49	3.58
VULC	Vulcan	50.367	247.020	57.65	308.84	3.55
WEYB	Weyburn	49.693	256.200	58.55	320.93	3.73
WGRY	Wells Gray	51.883	239.974	57.75	299.92	3.57
"""


class Info(namedtuple('Info', 'name lat lon mlat mlon L')):
    pass


def get_station_info(station_info=STATION_INFO):
    """
    Parse CARISMA station information record and return a mapping
    between site IDs and :class:`Info`.
    """
    info_map = OrderedDict()
    for line in station_info.splitlines():
        toks = line.split('\t')
        key = toks[0]
        name = toks[1]
        lat, lon, mlat, mlon = map(float, toks[2:6])
        L = float('nan') if toks[6] == 'N/A' else float(toks[6])
        info_map[key] = Info(name, lat, lon, mlat, mlon, L)
    return info_map
