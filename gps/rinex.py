from __future__ import division

from datetime import timedelta
from collections import OrderedDict, namedtuple

import scipy.constants as const

from constants import EPOCH, F_1, F_2, LAMBDA_1, LAMBDA_2

"""
Ideally we would have a cython interface to the GPSTk RINEX reader
routines.
"""

RINDUMP_OBS_MAP = {'GC1C': 'C1',
                   'GC1P': 'P1',
                   'GL1C': 'L1',
                   'GC2W': 'P2',
                   'GL2X': 'L2',
                   'ELE':  'el'}


LAMBDA_WL = const.c / (F_1 - F_2)
""" ??? """

ALPHA = (F_1 / F_2)**2
""" ??? """

MP_A = 1 + 2 / (ALPHA - 1)
""" ??? """

MP_B = 2 / (ALPHA - 1)
""" ??? """

MP_C = 2 * ALPHA / (ALPHA - 1)
""" ??? """

MP_D = 2 * ALPHA / (ALPHA - 1) - 1
""" ??? """


class Observation(namedtuple('Observation', 'C1 P1 L1 P2 L2 el')):
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
    def MP1(self):
        return self.P1 - MP_A * self.L1m + MP_B * self.L2m

    @property
    def MP2(self):
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
                for x in ['C1', 'P1', 'L1', 'P2', 'L2', 'el']:
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
    obs_map = read_rindump('/tmp/test/jplm0010.14o.dump')

    print(len(obs_map))
    print(obs_map['G01'].values()[0])

    print(obs_map['G01'].values()[0].L1m)

    # import numpy as NP
    # print(sorted(set(NP.diff(data_map['G03'].keys()))))
    # print(data_map['G03'].values()[:3])
