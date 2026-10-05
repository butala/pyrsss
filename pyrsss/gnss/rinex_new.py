import logging
from collections import OrderedDict, defaultdict

import pandas as pd

from .constants import week_sec2dt
from .preprocess import normalize_rinex
from .rinex import dump_rinex, RINDUMP_OBS_MAP


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


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    # rinex_fname = '/Users/butala/Documents/IARPA/absolute_tec/jplm0010.14o'
    normalized_rinex_fname = '/tmp/jplm0010.14o'

    # normalize_rinex(normalized_rinex_fname,
    #                 rinex_fname)

    dump_fname = '/tmp/jplm0010.14o.dump'
    nav_fname = '/Users/butala/Documents/IARPA/absolute_tec/jplm0010.14n'

    # dump_rinex(dump_fname,
               # normalized_rinex_fname,
               # nav_fname)

    dump = RinexDump.load(dump_fname)

    print(dump.stn)
    print(dump.xyz)
    print(dump.llh)
    # print(type(dump))
    # print(dump.recv_type)
    # print(dump.recv_p1c1)
    # print(dump.p1c1_table)
    dump.info()
