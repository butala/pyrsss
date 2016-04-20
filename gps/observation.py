import logging
from datetime import timedelta
from collections import OrderedDict, namedtuple, Iterator

import scipy.constants as const
from tables import open_file, IsDescription, Time64Col, Float64Col

from ..util.date import UNIX_EPOCH
from constants import F_1, F_2, LAMBDA_1, LAMBDA_2

logger = logging.getLogger('pyrsss.gps.observation')



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


class ObsMapFlatIterator(Iterator):
    def __init__(self, obs_map):
        """ ??? """
        self.obs_map = obs_map
        sorted_sats = sorted(obs_map)
        self.obs_dts = OrderedDict([(x, self.obs_map[x].keys()) for x in sorted_sats])

    def next(self):
        """ ??? """
        front_dts = [None if len(x) == 0 else x[0] for x in self.obs_dts.itervalues()]
        if all([x is None for x in front_dts]):
            raise StopIteration
        # below is the argmin function that ignores entries that are
        # None
        I, min_dt = min(filter(lambda x: x[1] is not None,
                               enumerate(front_dts)),
                        key=lambda x: x[1])
        sat = self.obs_dts.keys()[I]
        self.obs_dts[sat].pop(0)
        return min_dt, sat, self.obs_map[sat][min_dt]


class ObsTimeSeries(OrderedDict):
    def __setitem__(self, key, value):
        """ ??? """
        super(ObsTimeSeries, self).__setitem__(key,
                                               Observation(*value))


class ObsMap(dict):
    def __init__(self, h5_fname=None):
        """ ??? """
        super(ObsMap, self).__init__()
        if h5_fname:
            self.undump(h5_fname)

    def __missing__(self, key):
        prn = int(key[1:])
        self[key] = ObsTimeSeries()
        return self[key]

    """ ??? """
    class Table(IsDescription):
        dt   = Time64Col()
        C1   = Float64Col()
        P1   = Float64Col()
        P2   = Float64Col()
        L1   = Float64Col()
        L2   = Float64Col()
        az   = Float64Col()
        el   = Float64Col()
        satx = Float64Col()
        saty = Float64Col()
        satz = Float64Col()

    def dump(self, h5_fname, title=''):
        """ ??? """
        h5file = open_file(h5_fname, mode='w', title=title)
        group = h5file.create_group('/', 'phase_arcs', 'Phase connected arcs')
        if hasattr(self, 'xyz'):
            group._v_attrs.xyz = self.xyz
        if hasattr(self, 'llh'):
            group._v_attrs.llh = self.llh
        for sat in sorted(self):
            assert sat[0] == 'G'
            table = h5file.create_table(group, sat, ObsMap.Table, 'GPS prn={} data'.format(sat[1:]))
            row = table.row
            for dt, obs in self[sat].iteritems():
                row['dt'] = (dt - UNIX_EPOCH).total_seconds()
                row['C1'] = obs.C1
                row['P1'] = obs.P1
                row['P2'] = obs.P2
                row['L1'] = obs.L1
                row['L2'] = obs.L2
                row['az'] = obs.az
                row['el'] = obs.el
                row['satx'] = obs.satx
                row['saty'] = obs.saty
                row['satz'] = obs.satz
                row.append()
            table.flush()
        return h5_fname

    def undump(self, h5_fname):
        """ ??? """
        h5file = open_file(h5_fname, mode='r')
        group = h5file.root.phase_arcs
        try:
            self.xyz = group._v_attrs.xyz
        except:
            logger.warning('{} does not contain XYZ position'.format(h5_fname))
        try:
            self.llh = group._v_attrs.llh
        except:
            logger.warning('{} does not contain LLH position'.format(h5_fname))
        for table in group:
            sat = table.name
            for row in table.iterrows():
                dt = UNIX_EPOCH + timedelta(seconds=row['dt'])
                obs = Observation(*[row[x] for x in ['C1',
                                                     'P1',
                                                     'P2',
                                                     'L1',
                                                     'L2',
                                                     'az',
                                                     'el',
                                                     'satx',
                                                     'saty',
                                                     'satz']])
                self[sat][dt] = obs
        return self
