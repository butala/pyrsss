import os
import sys
import logging
from datetime import datetime
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import timedelta
from cStringIO import StringIO

import sh

from .constants import GPS_EPOCH
from .path import GPSTK_BUILD_PATH
from .teqc import rinex_info
from .preprocess import normalize_rinex
from .observation import Observation, ObsTimeSeries, ObsMap
from .receiver_types import ReceiverTypes
from .p1c1 import P1C1Table
from ..util.path import SmartTempDir, replace_path, tail

logger = logging.getLogger('pyrsss.gps.rinex')



GPS_RECEIVER_TYPES = ReceiverTypes()
"""
Global scope table of GPS receiver types.
"""


P1C1_TABLE = P1C1Table()
"""
Global scope table of CODE derived P1-C1 DCBs.
"""


RIN_DUMP = os.path.join(GPSTK_BUILD_PATH,
                        'core',
                        'apps',
                        'Rinextools',
                        'RinDump')



RINDUMP_OBS_MAP = {'GC1C': 'C1',
                   'GC1W': 'P1',
                   'GL1C': 'L1',
                   'GC2W': 'P2',
                   'GL2W': 'L2',
                   'RC1C': 'C1',
                   'RC1P': 'P1',
                   'RL1C': 'L1',
                   'RC2P': 'P2',
                   'RL2C': 'L2',
                   'ELE':  'el',
                   'AZI':  'az',
                   'SVX':  'satx',
                   'SVY':  'saty',
                   'SVZ':  'satz'}
"""
???

Unsure why the above do not correlate with the output of RinSum (no,
they do!).
"""


GPS_KEYS = ['GC1C', 'GC1W', 'GL1C', 'GC2W', 'GL2W', 'ELE', 'AZI', 'SVX', 'SVY', 'SVZ']


def fname2date(rinex_fname):
    """
    Return the :class:`datetime` associated with the RIENX file
    *rinex_fname* named according to the standard convention.
    """
    basename = os.path.basename(rinex_fname)
    doy = basename[4:7]
    daily_or_hour = basename[7]
    yy = basename[9:11]
    dt = datetime.strptime(doy + yy, '%j%y')
    if daily_or_hour == '0':
        return dt
    elif daily_or_hour in [chr(x) for x in range(ord('a'), ord('x') + 1)]:
        return dt + timedelta(hours=ord(daily_or_hour) - ord('a'))
    else:
        raise ValueError('could not parse date from RINEX file name '
                         '{}'.format(rinex_fname))


def get_receiver_position(rinex_fname, nav_fname):
    """
    Return the receiver position (via teqc) from information found in
    *rinex_fname* and *nav_fname*.
    """
    return rinex_info(rinex_fname, nav_fname)['xyz']


def get_receiver_type(rinex_fname):
    """
    Return the receiver type (header line REC # / TYPE / VERS) found
    in *rinex_fname*.
    """
    with open(rinex_fname) as fid:
        for line in fid:
            if line.rstrip().endswith('END OF HEADER'):
                break
            elif line.rstrip().endswith('REC # / TYPE / VERS'):
                return line[20:40].strip()
    raise ValueError('receiver type not found in header of RINEX file '
                     '{}'.format(rinex_fname))


def append_station_id(dump_fname,
                      rinex_fname):
    """
    Append the line "# Station ID: {station_id}" to *dump_fname*. The
    station identifier is determined from *rinex_fname*. Return
    *dump_fname*.
    """
    with open(dump_fname, 'a') as fid:
        stn_id = os.path.basename(rinex_fname)[:4]
        fid.write('# Station ID: {}\n'.format(stn_id))
    return dump_fname


def append_receiver_type(dump_fname,
                         rinex_fname,
                         receiver_types=GPS_RECEIVER_TYPES):
    """
    Append the lines "# Receiver type: {receiver_type}" and "#
    Receiver p1c1 type: {p1c1_type}" from *rinex_fname* to
    *dump_fname*. Return *dump_fname*.
    """
    with open(dump_fname, 'a') as fid:
        receiver_type = get_receiver_type(rinex_fname)
        fid.write('# Receiver type: {}\n'.format(receiver_type))
        fid.write('# Receiver p1c1 type: '
                  '{}\n'.format(receiver_types[receiver_type].c1p1))
    return dump_fname


def append_p1c1_date_table(dump_fname,
                           p1c1_date_table):
    """
    Append lines "# P1-C1 [m]: G{prn}: {p1c1_bias}" (in [m]) to
    *dump_fname* from the information found in
    *p1c1_date_table*. Return *dump_fname*.
    """
    with open(dump_fname, 'a') as fid:
        for prn in sorted(p1c1_date_table['prn']):
            fid.write('# P1-C1 [m]: G{:02d}: '
                      '{:8.3f}\n'.format(prn, p1c1_date_table['prn'][prn]))
    return dump_fname


def dump_rinex(dump_fname,
               rinex_fname,
               nav_fname,
               data_keys=GPS_KEYS,
               p1c1_table=P1C1_TABLE,
               receiver_position=None,
               rin_dump=RIN_DUMP):
    """
    ???

    currently only dumps GPS observables

    receiver position in [m]
    """
    rin_dump_command = sh.Command(rin_dump)
    stderr_buffer = StringIO()
    if receiver_position is None:
        receiver_position = get_receiver_position(rinex_fname,
                                                  nav_fname)
    logger.info('dumping {} to {}'.format(rinex_fname,
                                          dump_fname))
    args = ['--nav', nav_fname,
            '--ref', ','.join(map(str, receiver_position)),
            rinex_fname] + data_keys
    rin_dump_command(*args,
                     _out=dump_fname,
                     _err=stderr_buffer)
    stderr = stderr_buffer.read()
    if len(stderr) > 0:
        raise RuntimeError('error dumping the contents of '
                           '{} with {} ({})'.format(rinex_fname,
                                                    rin_dump,
                                                    stderr))
    append_station_id(dump_fname,
                      rinex_fname)
    append_receiver_type(dump_fname,
                         rinex_fname)
    append_p1c1_date_table(dump_fname,
                           p1c1_table(fname2date(rinex_fname)))
    return dump_fname


"""
MAKE CONFIG ROBUST

To incorporate RinDump we need:
- RINEX obs file
- RINEX nav file
- site location (ideally from a more robust source than RINEX header)
- desired observables (hook into GPSTk code that queries RINEX for, e.g., L1 and returns GL1C as appropriate)
"""


"""
Can use GPSTk to get robust RINEX information (e.g., interval, receiver type, data type mapping):
"""


"""
Can use GPSTk to find station position, e.g.:

./build-shaolin-v2.8/core/apps/positioning/PRSolve --obs ~/src/absolute_tec/jplm0010.14o --nav ~/src/absolute_tec/jplm0010.14n --sol GPS:12:W

Note that teqc also does this!
"""

"""
Note the GPSTk can compute ionospheric pierce points (see
core/lib/GNSSCore/Position.cpp). It would not be hard to convert this
t pure python (its just trig and coordinate transformations).
"""


"""
Ideally we would have a cython interface to the GPSTk RINEX reader
routines.
"""


class P1C1ObsTimeSeries(ObsTimeSeries):
    def __init__(self, receiver_type, p1c1_bias, replace_p1_with_c1=True):
        """ ??? """
        super(P1C1ObsTimeSeries, self).__init__()
        self.receiver_type = receiver_type
        self.p1c1_bias = p1c1_bias
        self.replace_p1_with_c1 = replace_p1_with_c1

    def __missing__(self, key):
        prn = int(key[1:])
        self[key] = ObsTimeSeries(self.receiver_p1c1_type,
                                  self.p1c1_table[prn])
        return self[key]

    def __setitem__(self, key, value):
        """ ??? """
        if self.receiver_type == 1:
            # C1 -> C1 + b
            # P2 -> P2 + b
            if value[0] != 0.0:
                value[0] += self.p1c1_bias
            if value[2] != 0.0:
                value[2] += self.p1c1_bias
        elif self.receiver_type == 2:
            # C1 -> C1 + b
            if value[0] != 0.0:
                value[0] += self.p1c1_bias
        elif self.receiver_type == 3:
            pass
        else:
            raise ValueError('unknown receiver type {}'.format(self.receiver_type))
        if value[1] == 0.0 and self.replace_p1_with_c1:
            # replace P1 with C1 (with bias correction if necessary)
            value[1] = value[0]
        # replace empty values (==0.0) with None
        value = [None if x == 0.0 else x for x in value]
        super(P1C1ObsTimeSeries, self).__setitem__(key,
                                                   Observation(*value))


class P1C1ObsMap(ObsMap):
    def __init__(self, receiver_type, receiver_p1c1_type, p1c1_table, h5_fname=None):
        """ ??? """
        super(P1C1ObsMap, self).__init__(h5_fname=h5_fname)
        self.receiver_type = receiver_type
        if receiver_p1c1_type not in [1, 2, 3]:
            raise ValueError('receiver P1-C1 type {} is unknown '
                             '(should be 1, 2, or 3 --- see the '
                             'GPS_Receiver_Type file header)')
        self.receiver_p1c1_type = receiver_p1c1_type
        self.p1c1_table = p1c1_table

    def __missing__(self, key):
        prn = int(key[1:])
        self[key] = P1C1ObsTimeSeries(self.receiver_p1c1_type,
                                      self.p1c1_table[prn])
        return self[key]


def read_rindump_footer(rindump_fname):
    """
    Return the receiver tuple and P1-C1 bias information found at the
    end of *rindump_fname*.
    """
    receiver_type = None
    receiver_p1c1_type = None
    p1c1_table = {}
    with open(rindump_fname) as fid:
        for line in tail(fid, window=100):
            if line.startswith('# Receiver type:'):
                receiver_type = line[17:].strip()
            elif line.startswith('# Receiver p1c1 type:'):
                receiver_p1c1_type = int(line[21:].strip())
            elif line.startswith('# P1-C1 [m]:'):
                prn = int(line[14:16])
                p1c1_table[prn] = float(line[17:])
            else:
                continue
    return receiver_type, receiver_p1c1_type, p1c1_table


def read_rindump(rindump_fname):
    """
    ???
    """
    obs_map = P1C1ObsMap(*read_rindump_footer(rindump_fname))
    with open(rindump_fname) as fid:
        for line in fid:
            if line.startswith('# wk'):
                # data header line
                cols = line.split()
                data_index_map = {RINDUMP_OBS_MAP[data_id]: i for i, data_id in enumerate(cols[4:])}
                column_mapping = []
                for x in Observation._fields:
                    try:
                        column_mapping.append(data_index_map[x])
                    except KeyError:
                        raise RuntimeError('could not find {} observable in {}'.format(x, rindump_fname))
                def reorder(l):
                    return [l[i] for i in column_mapping]
            elif line.startswith('# Refpos'):
                cols = line.split()
                # [m, m, m]
                obs_map.xyz = map(float, cols[3:6])
                lat = float(cols[8][:-1])
                lon = float(cols[9][:-1])
                if lon > 180:
                    lon -= 360
                alt = float(cols[10])
                # [deg, deg, m]
                obs_map.llh = [lat, lon, alt]
            elif line.startswith('#'):
                # skip other header lines
                continue
            else:
                cols = line.split()
                gps_week = int(cols[0])
                seconds = float(cols[1])
                sat = cols[2]
                dt = GPS_EPOCH + timedelta(days=7 * gps_week,
                                       seconds=seconds)
                obs_map[sat][dt] = reorder(map(float, cols[3:]))
    return obs_map


def dump_preprocessed_rinex(dump_fname,
                            obs_fname,
                            nav_fname,
                            work_path=None,
                            decimate=None):
    """
    Dump RINEX *obs_fname* and *nav_fname* to *dump_fname*. Preprocess
    the RNIEX file (i.e., normalization). Use *work_path* for
    intermediate files (use an automatically cleaned up area if not
    specified). Reduce the time interval to *decimate* [s] if
    given. Return *dump_fname*.
    """
    with SmartTempDir(work_path) as work_path:
        output_rinex_fname = replace_path(work_path, obs_fname)
        normalize_rinex(output_rinex_fname,
                        obs_fname,
                        decimate=decimate)
        dump_rinex(dump_fname,
                   output_rinex_fname,
                   nav_fname)
    return dump_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Dump RINEX observation file.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('dump_fname',
                        type=str,
                        help='output dump file file.')
    parser.add_argument('obs_fname',
                        type=str,
                        help='input RINEX obs file.')
    parser.add_argument('nav_fname',
                        type=str,
                        help='input RINEX nav file.')
    preprocess = parser.add_argument_group('RINEX preprocessing options')
    preprocess.add_argument('--decimate',
                            '-d',
                            type=int,
                            default=None,
                            help='decimate to time interval in [s]')
    args = parser.parse_args(argv[1:])

    dump_preprocessed_rinex(args.dump_fname,
                            args.obs_fname,
                            args.nav_fname,
                            decimate=args.decimate)

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)
    sys.exit(main())
