import logging
from collections import namedtuple
from datetime import timedelta
from itertools import groupby

import numpy as NP

from rinex import ObsTimeSeries

logger = logging.getLogger('pyrsss.gps.level')


class Config(namedtuple('Config',
                        'minimum_elevation '
                        'minimum_arc_time '
                        'minimum_arc_points '
                        'scatter_factor '
                        'scatter_threshold')):
    pass


MINIMUM_ELEVATION = 10
"""
??? (in [deg])
"""

MINIMUM_ARC_TIME = 18 * 60
"""
??? (in [s])
"""

MINIMUM_ARC_POINTS = 3
"""
??? (in [#])
"""

SCATTER_FACTOR = 1.6
"""
??? (in [#])
"""

SCATTER_THRESHOLD = 20
"""
??? (in [???])
"""


DEFAULT_CONFIG = Config(MINIMUM_ELEVATION,
                        MINIMUM_ARC_TIME,
                        MINIMUM_ARC_POINTS,
                        SCATTER_FACTOR,
                        SCATTER_THRESHOLD)


def arc_iter(obs_time_series, gap_length):
    """
    """
    arc_label = []
    arc_index = 0
    for i, (dt, obs) in enumerate(obs_time_series.iteritems()):
        if i > 0 and dt - obs_time_series.keys()[i - 1] > gap_length:
            arc_index += 1
        arc_label.append(arc_index)
    for arc_index_i, group in groupby(zip(arc_label,
                                          obs_time_series.iteritems()),
                                      key=lambda x: x[0]):
        yield arc_index_i, ObsTimeSeries(zip(*group)[1])


def level_phase_to_code(obs_map,
                        gap_length=timedelta(minutes=10),
                        config=DEFAULT_CONFIG):
    """
    ???

    gap_length should come from phase_edit?
    """
    for sat in sorted(obs_map):
        for arc_index, obs_time_series in arc_iter(obs_map[sat], gap_length):
            arc_time_length = (obs_time_series.keys()[-1] -
                               obs_time_series.keys()[0]).total_seconds()
            if arc_time_length < config.minimum_arc_time:
                # reject short arc (time)
                logger.info('rejecting sat={} arc={} --- '
                            'begin={:%Y-%m-%d %H:%M:%S} '
                            'end={:%Y-%m-%d %H:%N:%S} '
                            'length={} [s] '
                            '< {} [s]'.format(sat,
                                              arc_index,
                                              obs_time_series.keys()[0],
                                              obs_time_series.keys()[-1],
                                              arc_time_length,
                                              config.minimum_arc_time))
                continue
            if len(obs_time_series) < config.minimum_arc_points:
                # reject short arc (number of epochs)
                logger.info('rejecting sat={} arc={} --- len(arc) = '
                            '{} < {}'.format(sat,
                                             arc_index,
                                             len(obs_time_series),
                                             config.minimum_arc_points))
                continue
            print(arc_index)
            print(len(obs_time_series))
            print(obs_time_series.items()[0])
            print(obs_time_series.items()[1])
            print(obs_time_series.items()[-2])
            print(obs_time_series.items()[-1])
            print('  ')
            print('  ')
        assert False



if __name__ == '__main__':
    from phase_edit import phase_edit, filter_obs_map
    from rinex import read_rindump

    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)

    rinex_fname = '/Users/butala/src/absolute_tec/jplm0010.14o'
    interval = 30

    (time_reject_map,
     phase_adjust_map) = phase_edit(rinex_fname, interval)

    rinex_dump = '/Users/butala/src/absolute_tec/jplm0010.14o.dump'

    obs_map = read_rindump(rinex_dump)

    obs_map = filter_obs_map(obs_map,
                             time_reject_map,
                             phase_adjust_map)

    arc_map = level_phase_to_code(obs_map)
