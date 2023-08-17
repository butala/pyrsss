import logging
import sys
import math
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from datetime import timedelta
from collections import namedtuple, OrderedDict, Iterator
from datetime import timedelta
from itertools import groupby

import numpy as np
from tables import open_file, IsDescription, Time64Col, Float64Col
from more_itertools import peekable

from ..stats.stats import weighted_avg_and_std
from ..util.date import UNIX_EPOCH
from .constants import TECU_TO_M, M_TO_TECU
from .rms_model import RMSModel
from .observation import ObsMap, ObsTimeSeries

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
??? (in [TECU])
"""

P1P2_THRESHOLD = 5e-6
"""
??? (in in [m])
"""


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
    def timeiter(self, skip_fields=['L', 'L_scatter']):
        """
        ???
        """
        ArcPoint = namedtuple('ArcPoint', filter(lambda x: x not in skip_fields, self._fields))
        I_fields, fields = zip(*[(i, x) for i, x in enumerate(self._fields) if x not in skip_fields])
        for i in range(len(self.dt)):
            yield ArcPoint(*[self[i_field][i] for i_field in I_fields])


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


class ArcMapFlatIterator(Iterator):
    def __init__(self, arc_map):
        """ ??? """
        self.arc_map = arc_map
        sorted_sats = sorted(arc_map)
        self.flat_iters = [peekable(arc_map[x].flat) for x in sorted_sats]

    def next(self):
        """ ??? """
        front_entries = [x.peek(None) for x in self.flat_iters]
        if all([x is None for x in front_entries]):
            raise StopIteration
        # below is the argmin function that ignores entries that are
        # None
        I, _ = min(filter(lambda x: x[1] is not None,
                          enumerate(front_entries)),
                   key=lambda x: x[1].dt)
        return self.flat_iters[I].next()


"""
???
"""
class ArcList(list):
    @property
    def flat(self):
        for arc in self:
            for x in arc.timeiter():
                yield x


class ArcMap(OrderedDict):
    def __init__(self, h5_fname=None):
        """ ??? """
        super(ArcMap, self).__init__()
        if h5_fname:
            self.undump(h5_fname)

    def __missing__(self, key):
        """ ??? """
        self[key] = ArcList()
        return self[key]

    @property
    def flat(self):
        return ArcMapFlatIterator(self)

    """ ??? """
    class Table(IsDescription):
        dt   = Time64Col()
        stec = Float64Col()
        sprn = Float64Col()
        az   = Float64Col()
        el   = Float64Col()
        satx = Float64Col()
        saty = Float64Col()
        satz = Float64Col()

    def dump(self, h5_fname):
        """ ??? """
        h5file = open_file(h5_fname, mode='w', title='pyrsss.gps.level output')
        leveled_phase_arcs_group = h5file.create_group('/',
                                                       'leveled_phase_arcs',
                                                       'Leveled phase connected arcs')
        if hasattr(self, 'xyz'):
            leveled_phase_arcs_group._v_attrs.xyz = self.xyz
        if hasattr(self, 'llh'):
            leveled_phase_arcs_group._v_attrs.llh = self.llh
        for sat in sorted(self):
            assert sat[0] == 'G'
            sat_group = h5file.create_group(leveled_phase_arcs_group,
                                            sat,
                                            'Leveled phase connected arcs for {}'.format(sat))
            for i, leveled_arc in enumerate(self[sat]):
                table = h5file.create_table(sat_group,
                                            'arc' + str(i),
                                            ArcMap.Table,
                                            'GPS prn={} arc={} data'.format(sat, i))
                table.attrs.L = leveled_arc.L
                table.attrs.L_scatter = leveled_arc.L_scatter
                row = table.row
                for j in range(len(leveled_arc.dt)):
                    row['dt'] = (leveled_arc.dt[j] - UNIX_EPOCH).total_seconds()
                    row['stec'] = leveled_arc.stec[j]
                    row['sprn'] = leveled_arc.sprn[j]
                    row['az'] = leveled_arc.az[j]
                    row['el'] = leveled_arc.el[j]
                    row['satx'] = leveled_arc.satx[j]
                    row['saty'] = leveled_arc.saty[j]
                    row['satz'] = leveled_arc.satz[j]
                    row.append()
                table.flush()
        h5file.close()
        return h5_fname

    def undump(self, h5_fname):
        """ ??? """
        h5file = open_file(h5_fname, mode='r')
        leveled_phase_arcs_group = h5file.root.leveled_phase_arcs
        try:
            self.xyz = leveled_phase_arcs_group._v_attrs.xyz
        except:
            logger.warning('{} does not contain XYZ position'.format(h5_fname))
        try:
            self.llh = leveled_phase_arcs_group._v_attrs.llh
        except:
            logger.warning('{} does not contain LLH position'.format(h5_fname))
        for sat_group in leveled_phase_arcs_group:
            sat = sat_group._v_name
            for arc_table in sat_group:
                dt = []
                stec = []
                sprn = []
                az = []
                el = []
                satx = []
                saty = []
                satz = []
                for row in arc_table.iterrows():
                    dt.append(UNIX_EPOCH + timedelta(seconds=row['dt']))
                    stec.append(row['stec'])
                    sprn.append(row['sprn'])
                    az.append(row['az'])
                    el.append(row['el'])
                    satx.append(row['satx'])
                    saty.append(row['saty'])
                    satz.append(row['satz'])
                self[sat].append(LeveledArc(dt,
                                            stec,
                                            sprn,
                                            az,
                                            el,
                                            satx,
                                            saty,
                                            satz,
                                            arc_table.attrs.L,
                                            arc_table.attrs.L_scatter))
        h5file.close()
        return self


def level_phase_to_code(obs_map,
                        gap_length=timedelta(minutes=10),
                        config=DEFAULT_CONFIG):
    """
    ???

    gap_length should come from phase_edit?

    THIS FUNCTION IS TOO LONG --- BREAK INTO COMPONENTS
    """
    arc_map = ArcMap()
    arc_map.xyz = obs_map.xyz
    arc_map.llh = obs_map.llh
    rms_model = RMSModel()
    for sat in sorted(obs_map):
        # NO NO NO!!! Breaking arcs by gap length is wrong! Instead,
        # parse the log generated by DiscFix to determine the phase
        # breaks. The output of phase_edit should be phase connected
        # arcs suitable for leveling --- there is no reason to further
        # break up arcs at this point!
        for arc_index, obs_time_series in arc_iter(obs_map[sat], gap_length):
            logger.info('processing sat={} arc={}'.format(sat, arc_index))
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
            # remove observations for which P1, P2, L1, or L2 are nan
            valid_filter = lambda x: not any(np.isnan([getattr(x[1], attr) for attr in ['P1', 'P2', 'L1', 'L2']]))
            # remove measurements with |p1 - p2| < threshold
            p1p2_filter = lambda x: abs(x[1].P1 - x[1].P2) > config.p1p2_threshold
            filter_map = OrderedDict(el_filter=el_filter,
                                     valid_filter=valid_filter,
                                     p1p2_filter=p1p2_filter)
            dts_obs = obs_time_series.items()
            for name, obs_filter in filter_map.iteritems():
                dts_obs = filter(obs_filter, dts_obs)
                if len(dts_obs) == 0:
                    break
            if len(dts_obs) == 0:
                logger.info('rejecting sat={} arc={} after {} pass'.format(sat,
                                                                           arc_index,
                                                                           name))
                continue

            dts, obs = zip(*dts_obs)

            P_I = np.array([x.P_I for x in obs])
            L_Im = np.array([x.L_Im for x in obs])
            diff = P_I - L_Im
            modeled_var = (np.array(map(rms_model,
                                        [x.el for x in obs])) * TECU_TO_M)**2
            # compute level, level scatter, and modeled scatter
            N = len(diff)
            L, L_scatter = weighted_avg_and_std(diff, 1/modeled_var)
            sigma_scatter = np.sqrt(np.sum(modeled_var) / N)
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


def parse_override(config_overrides, config):
    """
    ???
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


def level_process(output_h5_fname,
                  input_h5_fname,
                  config_overrides=[],
                  config=DEFAULT_CONFIG):
    """
    """
    logger.info('reading phase connected arcs from {}'.format(input_h5_fname))
    obs_map = ObsMap(input_h5_fname)
    logger.info('beginning level phase to code process')
    config = parse_override(config_overrides, config)
    logger.info('{}'.format(config))
    arc_map = level_phase_to_code(obs_map, config=config)
    logger.info('storing leveled phase arcs to {}'.format(output_h5_fname))
    arc_map.dump(output_h5_fname)
    return output_h5_fname


def get_epilog(config=DEFAULT_CONFIG,
               config_units=CONFIG_UNITS):
    """
    """
    output = 'Default configuration names, values, and units:\n'
    for name, value in config._asdict().iteritems():
        output += '\t{}={}\t[{}]\n'.format(name, value, config_units[name])
    return output



def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Level GPS phase to code.',
                            formatter_class=RawDescriptionHelpFormatter,
                            epilog=get_epilog())
    parser.add_argument('leveled_arc_h5_fname',
                        type=str,
                        help='output H5 file containing leveled phase arcs')
    parser.add_argument('phase_edit_h5_fname',
                        type=str,
                        help='input H5 file generated by pyrsss.gps.phase_edit')
    parser.add_argument('--config',
                        '-c',
                        type=str,
                        nargs='+',
                        default=[],
                        help='leveling configuration overrides (specify as, e.g., minimum_elevation=15)')
    args = parser.parse_args(argv[1:])

    level_process(args.leveled_arc_h5_fname,
                  args.phase_edit_h5_fname,
                  config_overrides=args.config)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)
    sys.exit(main())
