import os
import sys
import logging
from datetime import datetime
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import timedelta
from io import StringIO

import sh

from collections import OrderedDict, defaultdict

import pandas as pd

from .constants import week_sec2dt
from .path import get_gnsstk_build_path
from .teqc import rinex_info
from .preprocess import normalize_rinex
from .receiver_types import ReceiverTypes
from .p1c1 import P1C1Table
from ..util.path import SmartTempDir, replace_path, tail

logger = logging.getLogger('pyrsss.gps.rinex')



_gps_receiver_types = None


def get_receiver_types():
    """
    Return the shared :class:`ReceiverTypes` table (built on first use
    so that import does not trigger a network fetch).
    """
    global _gps_receiver_types
    if _gps_receiver_types is None:
        _gps_receiver_types = ReceiverTypes()
    return _gps_receiver_types
"""
Global scope table of GPS receiver types.
"""


_p1c1_table = None


def get_p1c1_table():
    """
    Return the shared :class:`P1C1Table` (built on first use so that
    import does not trigger a network fetch).
    """
    global _p1c1_table
    if _p1c1_table is None:
        _p1c1_table = P1C1Table()
    return _p1c1_table
"""
Global scope table of CODE derived P1-C1 DCBs.
"""


RIN_DUMP_RELPATH = os.path.join('core',
                                'apps',
                                'Rinextools',
                                'RinDump')
"""Path to the RinDump tool relative to the GNSSTk build directory."""


def get_rin_dump():
    """Return the full path to the RinDump tool."""
    return os.path.join(get_gnsstk_build_path(), RIN_DUMP_RELPATH)



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
                         receiver_types=None):
    """
    Append the lines "# Receiver type: {receiver_type}" and "#
    Receiver p1c1 type: {p1c1_type}" from *rinex_fname* to
    *dump_fname*. Return *dump_fname*.
    """
    if receiver_types is None:
        receiver_types = get_receiver_types()
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
               p1c1_table=None,
               receiver_position=None,
               rin_dump=None):
    """
    Run GNSSTk RinDump on *rinex_fname* and write the result to
    *dump_fname*. Currently only dumps GPS observables. Receiver
    position in [m] (*receiver_position* is parsed from the navigation
    file when not given).
    """
    if rin_dump is None:
        rin_dump = get_rin_dump()
    if p1c1_table is None:
        p1c1_table = get_p1c1_table()
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
- desired observables (hook into GNSSTk code that queries RINEX for, e.g., L1 and returns GL1C as appropriate)
"""


"""
Can use GNSSpTk to get robust RINEX information (e.g., interval, receiver type, data type mapping):
"""


"""
Can use GNSSTk to find station position, e.g.:

./build-shaolin-v2.8/core/apps/positioning/PRSolve --obs ~/src/absolute_tec/jplm0010.14o --nav ~/src/absolute_tec/jplm0010.14n --sol GPS:12:W

Note that teqc also does this!
"""

"""
Note the GNSSTk can compute ionospheric pierce points (see
core/lib/GNSSCore/Position.cpp). It would not be hard to convert this
t pure python (its just trig and coordinate transformations).
"""


"""
Ideally we would have a cython interface to the GNSSTk RINEX reader
routines.
"""


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


class RinexDump(pd.DataFrame):
    _metadata = ['xyz',  # in [m]
                 'llh',  # in [ddm]
                 'stn',
                 'recv_type',
                 'recv_p1c1',
                 'p1c1_table']

    @property
    def _constructor(self):
        return RinexDump

    @classmethod
    def load(cls, rindump_fname, replace_p1_with_c1=True, p1c1=True):
        """
        Parse a teqc RINEX dump file and return a :class:`RinexDump`.
        Apply P1C1 bias corrections when *p1c1* (see :func:`correct_p1c1`,
        which also honors *replace_p1_with_c1*).
        """
        with open(rindump_fname) as fid:
            columns = None
            # parse up to "# Data" line
            for line in fid:
                if line.startswith('# Data'):
                    columns = ['gps_time', 'sat'] + [RINDUMP_OBS_MAP[x] for x in line.rstrip().split(' ')[2:]]
                    break
            if columns is None:
                raise ValueError('# Data line not found in {}'.format(rindump_fname))
            data_map = defaultdict(list)
            p1c1_table = OrderedDict()
            # parse remaining lines
            for line in fid:
                if line.startswith('# Refpos'):
                    toks = line.split(' ')
                    assert toks[2] == 'XYZ(m):'
                    assert toks[7] == 'LLH(ddm):'
                    xyz = list(map(float, toks[3:6]))
                    llh = toks[8:11]
                    assert llh[0][-1] == 'N'
                    assert llh[1][-1] == 'E'
                    llh = list(map(float, [llh[0][:-1],
                                          llh[1][:-1],
                                          llh[2]]))
                elif line.startswith('# Station ID:'):
                    stn = line.split(':')[1].strip()
                elif line.startswith('# Receiver type:'):
                    recv_type = line[17:].rstrip()
                elif line.startswith('# Receiver p1c1 type:'):
                    recv_p1c1 = int(line[22:])
                elif line.startswith('# P1-C1 [m]:'):
                    toks = line.split()
                    p1c1_table[toks[3][:-1]] = float(toks[4])
                elif line.startswith('#'):
                    # skip comment lines
                    pass
                else:
                    # parse data line
                    toks = line.replace(' 0.000 ', ' nan ').split()
                    gps_week = int(toks[0])
                    seconds = float(toks[1])
                    gps_time = week_sec2dt(gps_week, seconds)
                    sat = toks[2]
                    data = list(map(float, toks[3:]))
                    for column, data_i in zip(columns, [gps_time, sat] + data):
                        data_map[column].append(data_i)
            rinex_dump = cls(columns=columns, data=data_map)
            rinex_dump.xyz = xyz
            rinex_dump.llh = llh
            rinex_dump.stn = stn
            rinex_dump.recv_type = recv_type
            rinex_dump.recv_p1c1 = recv_p1c1
            rinex_dump.p1c1_table = p1c1_table
            if p1c1:
                correct_p1c1(rinex_dump)
            return rinex_dump


def correct_p1c1(rinex_dump, replace_p1_with_c1=True):
    """
    Apply the P1-C1 code bias table to *rinex_dump* (receiver types
    1--3). When *replace_p1_with_c1*, fill missing P1 with C1.
    """
    if rinex_dump.recv_p1c1 not in [1, 2, 3]:
        raise ValueError('unknown receiver type {} (must be 1, 2, or 3)'.format(rinex_dump.recv_p1c1))
    for sat in sorted(set(rinex_dump.sat)):
        b = rinex_dump.p1c1_table[sat]
        if rinex_dump.recv_p1c1 == 1:
            rinex_dump.loc[rinex_dump.sat == sat, 'C1'] += b
            rinex_dump.loc[rinex_dump.sat == sat, 'P2'] += b
        elif rinex_dump.recv_p1c1 == 2:
            rinex_dump.loc[rinex_dump.sat == sat, 'C1'] += b
    if replace_p1_with_c1:
        I = pd.isnull(rinex_dump['P1'])
        rinex_dump.loc[I, 'P1'] = rinex_dump.loc[I, 'C1']
    return rinex_dump
