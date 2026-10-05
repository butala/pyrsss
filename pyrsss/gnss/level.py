"""
Phase-level carrier phase to code for phase-connected arcs.
"""
import logging
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import namedtuple

import numpy as np
import pandas as pd

from ..stats.stats import weighted_avg_and_std
from .constants import (LAMBDA_1, LAMBDA_2, M_TO_TECU, TECU_TO_M,
                        glonass_lambda)
from .rms_model import RMSModel

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
"""Minimum elevation angle cut-off [deg]."""

MINIMUM_ARC_TIME = 18 * 60
"""Minimum arc time span [s]."""

MINIMUM_ARC_POINTS = 3
"""Minimum number of points per arc [#]."""

SCATTER_FACTOR = 1.6
"""Allowed leveling scatter in multiples of the modeled scatter [#]."""

SCATTER_THRESHOLD = 20
"""Maximum leveling uncertainty [TECU]."""

P1P2_THRESHOLD = 5e-6
"""Minimum |P1 - P2| used to detect bad code observations [m]."""


CONFIG_UNITS = {'minimum_elevation': 'deg',
                'minimum_arc_time': 's',
                'minimum_arc_points': '#',
                'scatter_factor': '#',
                'scatter_threshold': 'TECU',
                'p1p2_threshold': 'm'}


DEFAULT_CONFIG = Config(MINIMUM_ELEVATION,
                        MINIMUM_ARC_TIME,
                        MINIMUM_ARC_POINTS,
                        SCATTER_FACTOR,
                        SCATTER_THRESHOLD,
                        P1P2_THRESHOLD)


class LeveledArc(pd.DataFrame):
    """
    One phase-leveled arc: columns gps_time, az, el, satx/y/z, P_I, L_I
    ([TECU]) with station/satellite metadata and leveling results (L,
    L_scatter in [TECU]) attached as attributes.
    """
    _metadata = ['xyz',
                 'llh',
                 'stn',
                 'recv_type',
                 'sat',
                 'L',
                 'L_scatter']

    @property
    def _constructor(self):
        return LeveledArc


def convert_phase_m(df_arc, sat):
    """
    Return the (L1, L2) carrier phase of *df_arc* in [m]. Resolve
    satellite *sat* wavelengths (GLONASS frequencies vary by slot).
    """
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


def level(rinex_dump,
          config=DEFAULT_CONFIG):
    """
    Phase-level the arcs labeled in *rinex_dump* (see
    :func:`phase_edit.label_phase_arcs`) to code. Return the list of
    accepted :class:`LeveledArc` according to the *config* rejection
    rules.
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
        L_Im = L1m - L2m
        diff = P_I - L_Im
        modeled_var = (np.array([rms_model(el) for el in df_arc.el.values]) * TECU_TO_M)**2
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


def parse_override(config_overrides, config):
    """
    Return a :class:`Config` with the *name*=*value* entries of
    *config_overrides* applied to *config*.
    """
    config_map = config._asdict().copy()
    for token in config_overrides:
        name, value = token.split('=')
        if name not in config._fields:
            raise ValueError('unrecognized leveling configuration option {}'.format(name))
        if name in ['scatter_factor', 'p1p2_threshold']:
            value = float(value)
        else:
            value = int(value)
        config_map[name] = value
    return Config(*config_map.values())


def get_epilog(config=DEFAULT_CONFIG,
               config_units=CONFIG_UNITS):
    """
    Return the help epilog documenting the leveling configuration.
    """
    output = 'Default configuration names, values, and units:\n'
    for name, value in config._asdict().items():
        output += '\t{}={}\t[{}]\n'.format(name, value, config_units[name])
    return output


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Level GPS phase to code.',
                            formatter_class=RawDescriptionHelpFormatter,
                            epilog=get_epilog())
    parser.add_argument('output_fname',
                        type=str,
                        help='output pickle file containing leveled phase arcs')
    parser.add_argument('rinex_dump_fname',
                        type=str,
                        help='input pickle file containing an edited, arc-labeled RinexDump')
    parser.add_argument('--config',
                        '-c',
                        type=str,
                        nargs='+',
                        default=[],
                        help='leveling configuration overrides (specify as, e.g., minimum_elevation=15)')
    args = parser.parse_args(argv[1:])

    config = parse_override(args.config, DEFAULT_CONFIG)
    rinex_dump = pd.read_pickle(args.rinex_dump_fname)
    leveled_arcs = level(rinex_dump, config=config)
    pd.to_pickle(leveled_arcs, args.output_fname)
    return args.output_fname


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)
    sys.exit(main())
