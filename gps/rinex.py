from __future__ import division

import os
import logging
from datetime import timedelta
from collections import OrderedDict, namedtuple
from cStringIO import StringIO

import sh
import scipy.constants as const

from constants import EPOCH, F_1, F_2, LAMBDA_1, LAMBDA_2

logger = logging.getLogger('pyrsss.gps.rinex')


try:
    GPSTK_BUILD_PATH = os.environ['GPSTK_BUILD']
except KeyError:
    raise RuntimeError('environment variable GPSTK_BUILD not set')


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


def dump_rinex(dump_fname,
               rinex_fname,
               nav_fname,
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
        reciever_position = get_receiver_position(rinex_fname)
    logger.info('dumping {} to {}'.format(rinex_fname,
                                          dump_fname))
    args = ['--nav', nav_fname,
            '--ref', ','.join(map(str, receiver_position)),
            rinex_fname] + RINDUMP_OBS_MAP.keys()
    rin_dump_command(*args,
                     _out=dump_fname,
                     _err=stderr_buffer)
    stderr = stderr_buffer.read()
    if len(stderr) > 0:
        raise RuntimeError('error dumping the contents of '
                           '{} with {} ({})'.format(rinex_fname,
                                                    rin_dump,
                                                    stderr))
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


LAMBDA_WL = const.c / (F_1 - F_2)
"""
???
"""

ALPHA = (F_1 / F_2)**2
"""
???
"""

MP_A = 1 + 2 / (ALPHA - 1)
"""
???
"""

MP_B = 2 / (ALPHA - 1)
"""
???
"""

MP_C = 2 * ALPHA / (ALPHA - 1)
"""
???
"""

MP_D = 2 * ALPHA / (ALPHA - 1) - 1
"""
???
"""


class Observation(namedtuple('Observation',
                             'C1 P1 P2 L1 L2 az el satx saty satz')):
    @property
    def L1m(self):
        """
        Return carrier 1 phase in [m].
        """
        return self.L1 * LAMBDA_1

    @property
    def L2m(self):
        """
        Return carrier 2 phase in [m].
        """
        return self.L2 * LAMBDA_2

    @property
    def N_WL(self):
        """
        """
        return self.L1 - self.L2 - (F_1 * self.P1 + F_2 * self.P2) / (LAMBDA_WL * (F_1 + F_2))

    @property
    def L_WL(self):
        """
        """
        return LAMBDA_WL * self.N_WL

    @property
    def P_I(self):
        """
        """
        return self.P2 - self.P1

    @property
    def L_I(self):
        """
        """
        return self.L1 - self.L2

    @property
    def L_Im(self):
        """
        """
        return self.L1m - self.L2m

    @property
    def MP1(self):
        """
        """
        return self.P1 - MP_A * self.L1m + MP_B * self.L2m

    @property
    def MP2(self):
        """
        """
        return self.P2 - MP_C * self.L1m + MP_D * self.L2m


class ObsTimeSeries(OrderedDict):
    def __setitem__(self, key, value):
        super(ObsTimeSeries, self).__setitem__(key,
                                               Observation(*value))


class ObsMap(dict):
    def __missing__(self, key):
        self[key] = ObsTimeSeries()
        return self[key]


def read_rindump(rindump_fname):
    """
    ???
    """
    obs_map = ObsMap()
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
            elif line.startswith('#'):
                # skip header lines
                continue
            else:
                cols = line.split()
                gps_week = int(cols[0])
                seconds = float(cols[1])
                sat = cols[2]
                dt = EPOCH + timedelta(days=7 * gps_week,
                                       seconds=seconds)
                obs_map[sat][dt] = reorder(map(float, cols[3:]))
    return obs_map


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)

    dump_rinex('/tmp/jplm0010.14o.dump',
               '/Users/butala/src/absolute_tec/jplm0010.14o',
               '/Users/butala/src/absolute_tec/jplm0010.14n',
               receiver_position=[-2493304.6796, -4655215.1032, 3565497.5918])

    obs_map = read_rindump('/tmp/jplm0010.14o.dump')

    print(obs_map.values()[0].items()[0])
