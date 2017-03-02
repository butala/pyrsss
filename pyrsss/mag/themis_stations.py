from urllib2 import urlopen
from contextlib import closing
from collections import OrderedDict, namedtuple


INFO_URL = 'http://themis.ssl.berkeley.edu/gmag/gmag_groups.php'


class Info(namedtuple('Info', 'lat lon name mlat mlon')):
    pass


PARSE_MAP = {'ccode':   ('key', str),
             'lat':     ('lat', float),
             'lng':     ('lon', float),
             'name':    ('name', str),
             'mag_lat': ('mlat', float),
             'mag_lng': ('mlon', float)}


def get_station_info(info_url=INFO_URL, parse_map=PARSE_MAP):
    """
    Parse information for magnetometer sites that report data to the
    THEMIS project. Returns a mapping between station IDs and
    :class:`Info` regarding the site.
    """
    station_info = OrderedDict()
    with closing(urlopen(info_url)) as fid:
        stn_data = {}
        for line in fid:
            if line.startswith('};'):
                key = stn_data.pop('key')
                if 'mlat' not in stn_data:
                    stn_data['mlat'] = float('nan')
                if 'mlon' not in stn_data:
                    stn_data['mlon'] = float('nan')
                station_info[key] = Info(**stn_data)
                stn_data = {}
            line = line.lstrip()
            for search_key, (key, convert) in parse_map.iteritems():
                if line.startswith(search_key):
                    stn_data[key] = convert(line.split('"')[1])
    return station_info
