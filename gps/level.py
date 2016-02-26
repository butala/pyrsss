from __future__ import division

import logging
import math
from collections import namedtuple, defaultdict
from datetime import timedelta
from itertools import groupby

import numpy as NP

from ..util.stats import weighted_avg_and_std
from constants import TECU_TO_M, M_TO_TECU
from rinex import ObsTimeSeries
from rms_model import RMSModel

logger = logging.getLogger('pyrsss.gps.level')


class Config(namedtuple('Config',
                        'minimum_elevation '
                        'minimum_arc_time '
                        'minimum_arc_points '
                        'scatter_factor '
                        'scatter_threshold '
                        'p1p2_threshold')):
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

P1P2_THRESHOLD = 5e-6
"""
??? (in in [m])
"""


DEFAULT_CONFIG = Config(MINIMUM_ELEVATION,
                        MINIMUM_ARC_TIME,
                        MINIMUM_ARC_POINTS,
                        SCATTER_FACTOR,
                        SCATTER_THRESHOLD,
                        P1P2_THRESHOLD)


class LeveledArc(namedtuple('LeveledArc',
                            'dt '
                            'stec '
                            'sprn '
                            'az '
                            'el '
                            'satx '
                            'saty '
                            'satz '
                            'L '
                            'L_scatter')):
    pass


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
    arc_map = defaultdict(list)
    rms_model = RMSModel()
    for sat in sorted(obs_map):
        for arc_index, obs_time_series in arc_iter(obs_map[sat], gap_length):
            arc_time_length = (obs_time_series.keys()[-1] -
                               obs_time_series.keys()[0]).total_seconds()
            if arc_time_length < config.minimum_arc_time:
                # reject short arc (time)
                logger.info('rejecting sat={} arc={} --- '
                            'begin={:%Y-%m-%d %H:%M:%S} '
                            'end={:%Y-%m-%d %H:%M:%S} '
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
            # remove observations below minimum elevation limit
            el_filter = lambda x: x[1].el >= config.minimum_elevation
            # remove measurements with |p1 - p2| < threshold
            p1p2_filter = lambda x: abs(x[1].P1 - x[1].P2) > config.p1p2_threshold
            dts_obs = zip(*filter(lambda x: el_filter(x) and p1p2_filter(x),
                                  obs_time_series.iteritems()))
            if len(dts_obs) == 0:
                logger.info('rejecting sat={} arc={} --- el_filer or p1p2 '
                            'filter'.format(sat, arc_index))
                continue
            else:
                dts, obs = dts_obs
            P_I = NP.array([x.P_I for x in obs])
            L_Im = NP.array([x.L_Im for x in obs])
            diff = P_I - L_Im
            modeled_var = (NP.array(map(rms_model,
                                        [x.el for x in obs])) * TECU_TO_M)**2
            # compute level, level scatter, and modeled scatter
            N = len(diff)
            L, L_scatter = weighted_avg_and_std(diff, 1/modeled_var)
            sigma_scatter = NP.sqrt(NP.sum(modeled_var) / N)
            # check for excessive leveling uncertainty
            if L_scatter > config.scatter_factor * sigma_scatter:
                logger.info('rejecting sat={} arc={} --- L scatter={:.6f} '
                            '> {:.1f} * {:.6f}'.format(sat,
                                                       arc_index,
                                                       L_scatter,
                                                       config.scatter_factor,
                                                       sigma_scatter))
                continue
            if L_scatter / TECU_TO_M > config.scatter_threshold:
                logger.info('rejecting sat={} arc={} --- L uncertainty (in '
                            '[TECU])={:.1f} > '
                            '{:.1f}'.format(sat,
                                            arc_index,
                                            L_scatter * M_TO_TECU,
                                            config.scatter_threshold))
                continue
            # store information
            arc_map[sat].append(LeveledArc(dts,
                                           (L_Im + L) * M_TO_TECU,
                                           P_I * M_TO_TECU,
                                           [x.az for x in obs],
                                           [x.el for x in obs],
                                           [x.satx for x in obs],
                                           [x.saty for x in obs],
                                           [x.satz for x in obs],
                                           L * M_TO_TECU,
                                           L_scatter * M_TO_TECU))
    return arc_map


if __name__ == '__main__':
    from phase_edit import phase_edit, filter_obs_map
    from rinex import read_rindump, dump_rinex

    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)

    rinex_fname = '/Users/butala/src/absolute_tec/jplm0010.14o'
    nav_fname = '/Users/butala/src/absolute_tec/jplm0010.14n'
    interval = 30

    (time_reject_map,
     phase_adjust_map) = phase_edit(rinex_fname, interval)

    rinex_dump_fname = '/tmp/jplm0010.14o.dump'
    dump_rinex(rinex_dump_fname,
               rinex_fname,
               nav_fname)

    obs_map = read_rindump(rinex_dump_fname)

    obs_map = filter_obs_map(obs_map,
                             time_reject_map,
                             phase_adjust_map)

    arc_map = level_phase_to_code(obs_map)

    import matplotlib
    matplotlib.use('agg')
    import pylab as PL

    fig = PL.figure(figsize=(11, 8.5))

    for sat, arcs in arc_map.iteritems():
        fig.clf()
        for arc in arcs:
            PL.plot_date(arc.dt,
                         arc.sprn,
                         ls='None',
                         marker='x',
                         color='r')
            PL.plot_date(arc.dt,
                         arc.stec,
                         ls='-',
                         marker='None',
                         color='b')
        PL.savefig('/tmp/fig/{}.pdf'.format(sat),
                   bbox_inches='tight')
