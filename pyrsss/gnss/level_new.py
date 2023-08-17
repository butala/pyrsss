import logging

import numpy as np
import pandas as pd

from ..stats.stats import weighted_avg_and_std
from constants import LAMBDA_1, LAMBDA_2, TECU_TO_M, M_TO_TECU, glonass_lambda
from level import DEFAULT_CONFIG
from rms_model import RMSModel

logger = logging.getLogger('pyrsss.gps.level_new')


class LeveledArc(pd.DataFrame):
    _metadata = ['xyz',
                 'llh',
                 'stn',
                 'recv_type',
                 'stn',
                 'sat',
                 'L',
                 'L_scatter']

    @property
    def _constructor(self):
        return LeveledArc


def convert_phase_m(df_arc, sat):
    """
    """
    print(df_arc.shape)
    if sat[0] == 'G':
        return (df_arc.L1 * LAMBDA_1,
                df_arc.L2 * LAMBDA_2)
    elif sat[0] == 'R':
        dt = df_arc.iloc[0].gps_time
        slot = int(sat[1:])
        lambda1, lambda2 = glonass_lambda(slot, dt)
        return (df_arc.L1 * lambda1,
                df_arc.L2 * lambda2)
    else:
        raise ValueError('cannot convert phase to [m] for {}'.format(sat))
    assert False


def level(rinex_dump,
          config=DEFAULT_CONFIG):
    """
    """
    rms_model = RMSModel()
    leveled_arcs = []
    for arc_index, arc in enumerate(sorted(set(rinex_dump.arc))):
        df_arc = rinex_dump[rinex_dump.arc == arc]
        sat = df_arc.iloc[0].sat
        delta = df_arc.iloc[-1].gps_time - df_arc.iloc[0].gps_time
        arc_time_length = delta.total_seconds()
        if arc_time_length < config.minimum_arc_time:
            # reject short arc (time)
            logger.info('rejecting arc={} --- '
                        'begin={:%Y-%m-%d %H:%M:%S} '
                        'end={:%Y-%m-%d %H:%M:%S} '
                        'length={} [s] '
                        '< {} [s]'.format(arc,
                                          df_arc.iloc[0].gps_time,
                                          df_arc.iloc[-1].gps_time,
                                          arc_time_length,
                                          config.minimum_arc_time))
            continue
        if df_arc.shape[0] < config.minimum_arc_points:
            # reject short arc (number of epochs)
            logger.info('rejecting arc={} --- len(arc) = '
                        '{} < {}'.format(sat,
                                         arc,
                                         df_arc.shape[0],
                                         config.minimum_arc_points))
            continue
        # remove observations below minimum elevation limit
        I = df_arc.el >= config.minimum_elevation
        # remove observations for which P1, P2, L1, or L2 are nan
        I &= df_arc.P1.notnull()
        I &= df_arc.P2.notnull()
        I &= df_arc.L1.notnull()
        I &= df_arc.L2.notnull()
        # remove measurements with |p1 - p2| < threshold
        I &= abs(df_arc.P1 - df_arc.P2) > config.p1p2_threshold
        # compute geometry free combinations
        df_arc = df_arc.loc[I, :]
        if df_arc.shape[0] == 0:
            continue
        P_I = df_arc.P2 - df_arc.P1
        L1m, L2m = convert_phase_m(df_arc, sat)
        # THIS IS ONLY TRUE FOR GPS!!!
        # L1m = df_arc.L1 * LAMBDA_1
        # L2m = df_arc.L2 * LAMBDA_2
        L_Im = L1m - L2m
        diff = P_I - L_Im
        modeled_var = (np.array(map(rms_model,
                                    df_arc.el.values)) * TECU_TO_M)**2
        # compute level, level scatter, and modeled scatter
        N = len(diff)
        if N == 0:
            continue
        L, L_scatter = weighted_avg_and_std(diff, 1 / modeled_var)
        sigma_scatter = np.sqrt(np.sum(modeled_var) / N)
        # check for excessive leveling uncertainty
        if L_scatter > config.scatter_factor * sigma_scatter:
            logger.info('rejecting arc={} --- L scatter={:.6f} '
                        '> {:.1f} * {:.6f}'.format(arc_index,
                                                   L_scatter,
                                                   config.scatter_factor,
                                                   sigma_scatter))
            continue
        if L_scatter / TECU_TO_M > config.scatter_threshold:
            logger.info('rejecting arc={} --- L uncertainty (in '
                        '[TECU])={:.1f} > '
                        '{:.1f}'.format(arc_index,
                                        L_scatter * M_TO_TECU,
                                        config.scatter_threshold))
            continue
        # store information
        data_map = {'gps_time': df_arc.gps_time.values,
                    'az': df_arc.az.values,
                    'el': df_arc.el.values,
                    'satx': df_arc.satx.values,
                    'saty': df_arc.saty.values,
                    'satz': df_arc.satz.values,
                    'P_I': P_I * M_TO_TECU,
                    'L_I': (L_Im + L) * M_TO_TECU}
        leveled_arc = LeveledArc(data_map)
        leveled_arc.xyz = df_arc.xyz
        leveled_arc.llh = df_arc.llh
        leveled_arc.stn = df_arc.stn
        leveled_arc.recv_type = df_arc.recv_type
        leveled_arc.sat = sat
        leveled_arc.L = L * M_TO_TECU
        leveled_arc.L_scatter = L_scatter * M_TO_TECU
        leveled_arcs.append(leveled_arc)
    return leveled_arcs


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    pkl_fname = '/tmp/jplm0010.14o.pkl'
    rinex_dump = pd.read_pickle(pkl_fname)

    leveled_arcs = level(rinex_dump)

    # print(leveled_arcs[0])
    # print(leveled_arcs[10])

    import pylab as PL

    fig = PL.figure()
    for i, arc in enumerate(leveled_arcs):
        fig.clf()
        PL.plot_date(arc.gps_time, arc.P_I, marker='x', ls='None', color='r')
        PL.plot_date(arc.gps_time, arc.L_I, marker='None', ls='-', color='b')
        PL.title('arc = {}  sat = {}'.format(i, arc.sat))
        PL.savefig('/tmp/test/arc_{:03d}.pdf'.format(i), bbox_inches='tight')
