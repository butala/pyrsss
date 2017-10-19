from datetime import datetime
from contextlib import closing
from urllib2 import urlopen
from collections import OrderedDict
from enum import Enum

from intervals import DateTimeInterval



def skip_block_header_comments(fid):
    """
    Skip past the initial three line header of the response record
    found in *fid*
    """
    line1 = fid.next()
    line2 = fid.next()
    line3 = fid.next()
    assert line1 == line3 == '#\n'
    assert line2 == '#' * 83 + '\n'
    return fid


def parse_block_header(fid):
    """
    Parse the block header found in *fid*. Return a mapping containing
    the header information.
    """
    block_header = {}
    for line in fid:
        if line.startswith('#'):
            break
        toks = line.split()
        key = '_'.join(toks[1:-1])[:-1]
        value = toks[-1]
        if key.endswith('date'):
            value = datetime.strptime(value, '%Y,%j,%H:%M:%S')
        block_header[key] = value
    return block_header


class BlockType(Enum):
    """
    Block information type identifier.
    """
    eof = 0
    end = 1
    decimation = 2
    sensitivity = 3


def parse_block_info(fid):
    """
    Parse the decimation or sensitivity information record found in
    the block in *fid*. Return the tuple containing the
    :class:`BlockType` and the mapping of block information.
    """
    try:
        line1 = fid.next()
    except StopIteration:
        return BlockType.eof, {}
    if line1 == '#' * 83 + '\n':
        line2 = fid.next()
        assert line2 == '#\n'
        return BlockType.end, {}
    assert line1 == '#                  +-----------------------------------+\n'
    line2 = fid.next()
    if line2 == '#                  |            Decimation             |\n':
        block_type = BlockType.decimation
    elif line2 == '#                  |      Channel Sensitivity/Gain     |\n':
        block_type = BlockType.sensitivity
    else:
        raise RuntimeError('error parsing data block in {}'.format(fid.name))
    for i in range(4):
        fid.next()
    block_info = {}
    for line in fid:
        if line == '#\n':
            break
        toks = line.split()
        key = '_'.join(toks[1:-1])[:-1]
        try:
            val = int(toks[-1])
        except ValueError:
            val = float(toks[-1])
        block_info[key] = val
    return block_type, block_info


def parse_block(fid):
    """
    Parse a single information block found in *fid*. Return the tuple
    containing the block header and the list of block information.
    """
    # header
    block_header = parse_block_header(fid)
    block = []
    while True:
        block_type, block_info = parse_block_info(fid)
        if block_type == BlockType.eof:
            return None, None
        elif block_type == BlockType.end:
            return block_header, block
        elif block_type == BlockType.decimation:
            continue
        elif block_type == BlockType.sensitivity:
            block.append(block_info)
        else:
            assert False
    raise RuntimeError('error parsing {}'.format(fid.name))


class RespMap(OrderedDict):
    """
    Object in which to store IRIS station response information for a
    single site.
    """
    def __call__(self, dt):
        for key, value in self.iteritems():
            if dt in key:
                return value
        raise KeyError('{:%Y-%m-%d %H:%M:%S} does not intersect with any stored datetime interval'.format(dt))

    def sensitivity(self, dt, channel):
        for stage in self(dt)[channel]:
            if stage['Stage_sequence_number'] == 0:
                return stage['Sensitivity']
        assert False


def check(header, key, value):
    """
    If *var* is not `None` first check that *var* and *value*
    match. Return the parsed *value*.
    """
    if value is None:
        return header[key]
    assert header[key] == value
    return value


def parse_station_resp(fid):
    """
    Gather information from a single station IRIS response file
    *fid*. Return the information as a :class:`RespMap`.
    """
    resp_map = RespMap()
    # sanity check initialization
    network = None
    stn = None
    location = None
    # skip initial header comment block
    skip_block_header_comments(fid)
    while True:
        block_header, block = parse_block(fid)
        if block_header is None and block is None:
            break
        # sanity check (same network, station, and location across recorded blocks)
        network = check(block_header, 'Network', network)
        stn = check(block_header, 'Station', stn)
        location = check(block_header, 'Location', location)
        # store block information
        interval = DateTimeInterval.closed_open(block_header['Start_date'],
                                                block_header['End_date'])
        resp_map.setdefault(interval, {})[block_header['Channel']] = block
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
