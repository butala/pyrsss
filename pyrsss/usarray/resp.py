from datetime import datetime
from contextlib import closing
from urllib2 import urlopen

from intervals import DateTimeInterval


def check_same(var, value, parser):
    """
    If *var* is not `None` first check that *var* and *value* (first
    converted with *parser*) match. Return the parsed *value*.
    """
    parsed_value = parser(value)
    if var is not None:
        assert var == parsed_value
    return parsed_value


def check_same_date(var, value):
    """
    The version of :func:`check_same` tailored to resp file formatted
    date values.
    """
    parser = lambda x: datetime.strptime(x, '%Y,%j,%H:%M:%S')
    return check_same(var, value, parser)


def check_same_str(var, value):
    """
    The version of :func:`check_same` tailored to string values.
    """
    return check_same(var, value, lambda x: x)


def parse_station_resp(fid):
    """
    Gather information from a single station IRIS response file
    *fid*. Return the information as a mapping.
    """
    station = None
    network = None
    start_date = None
    end_date = None
    resp_map = {}
    for line in fid:
        if line.startswith('#'):
            continue
        toks = line.split()
        key = '_'.join(toks[1:-1])[:-1]
        value = toks[-1]
        if key == 'Channel':
            channel = value
        elif key == 'Station':
            station = check_same_str(station, value)
        elif key == 'Network':
            network = check_same_str(network, value)
        elif key == 'Start_date':
            start_date = check_same_date(start_date, value)
        elif key == 'End_date':
            end_date = check_same_date(end_date, value)
        elif key == 'Stage_sequence_number':
            stage_number = int(value)
        elif key == 'Sensitivity':
            if stage_number == 0:
                resp_map[channel] = float(value)
    resp_map['station'] = station
    resp_map['network'] = network
    resp_map['interval'] = DateTimeInterval.closed(start_date, end_date)
    return resp_map


def build_url(station, date):
    """
    Return the URL to fetch the response record for USArray MT station
    identifier *station* for *date*.
    """
    return 'http://service.iris.edu/irisws/resp/1/query?net=EM&sta={}&loc=--&cha=*&time={:%Y-%m-%dT%H:%M:%S}'.format(station, date)


def get_station_resp(station, date):
    """
    For the given USArray MT station ID *station* and *date*, return
    the instrument response function.
    """
    url = build_url(station, date)
    with closing(urlopen(url)) as fid:
        return parse_station_resp(fid)


if __name__ == '__main__':
    date = datetime(2013, 7, 2)
    station = 'KSP33'

    print(get_station_resp(station, date))
