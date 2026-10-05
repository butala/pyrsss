from urllib.request import urlopen
from io import TextIOWrapper
from contextlib import closing
from collections import namedtuple


INFO_URL = 'http://www.intermagnet.org/imos/imotblobs-eng.php'


class Info(namedtuple('Info', 'name country lat lon institute GIN')):
    pass


def get_value(fid):
    """
    Return the string value found on the next line of *fid*.
    """
    line = next(fid)
    return line.lstrip().split('<td>')[1].split('</td>')[0].replace('&deg;', '')



def get_stations(info_url=INFO_URL):
    """
    Parse INTERMAGET station information and return a mapping between
    IDs and :class:`Info`.
    """
    station_info = dict()
    with closing(TextIOWrapper(urlopen(info_url),
                           encoding="utf-8",
                           errors="replace")) as fid:
        record = None
        while True:
            try:
                line = next(fid)
            except StopIteration:
                break
            if 'iaga_code=' in line:
                if record:
                    pass
                toks = line.split('<span lang="">')
                key = toks[1].split('</span>')[0]
                name = next(fid).lstrip()
                name = name.split('<span lang="">')[1].split('</span>')[0]
                country = get_value(fid)
                colat = float(get_value(fid))
                lat = 90 - colat
                lon = float(get_value(fid))
                institute = next(fid).strip()[4:-5]
                GIN = get_value(fid)
                station_info[key] = Info(name,
                                         country,
                                         lat,
                                         lon,
                                         institute,
                                         GIN)
    return station_info
