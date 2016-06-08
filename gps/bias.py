import sys
import logging
import os
import posixpath
from collections import namedtuple, defaultdict
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import datetime

import numpy as NP
from scipy.interpolate import RectBivariateSpline

from constants import SHELL_HEIGHT, TECU_TO_NS
from level import LeveledArc, ArcMap
from util import shell_mapping
from ipp import cnv_azel2latlon
from teqc import rinex_info
from sideshow import update_sideshow_file
from ..ionex.read_ionex import interpolate2D_temporal
from ..util.search import find_le
from ..util.path import SmartTempDir

logger = logging.getLogger('pyrsss.gps.bias')

"""
Instead, call this unbias or debias?

Remove transmitter and receiver biases from phase-connected arcs. Use
Attila's method to estimate the receiver bias (using IGS IONEX records
containing satellite biases and modeled VTEC maps).
"""


class AugLeveledArc(namedtuple('AugLeveledArc',
                               ' '.join(LeveledArc._fields + ('el_map',
                                                              'ipp_lat',
                                                              'ipp_lon')))):
    def __new__(cls,
                leveled_arc,
                site_lat,
                site_lon,
                height=SHELL_HEIGHT):
        fields = leveled_arc._asdict()
        fields['el_map'] = [shell_mapping(el_i, h=height) for el_i in fields['el']]
        ipp_lat, ipp_lon = zip(*[cnv_azel2latlon(az_i,
                                                 el_i,
                                                 (site_lat, site_lon),
                                                 ht=height) for (az_i, el_i) in zip(fields['az'],
                                                                                    fields['el'])])
        fields['ipp_lat'] = ipp_lat
        fields['ipp_lon'] = ipp_lon
        return cls._make(fields.values())


def augmented_arc_map(arc_map,
                      site_lat,
                      site_lon,
                      height=SHELL_HEIGHT):
    """
    """
    augmented_map = defaultdict(list)
    for key, arc_list in arc_map.iteritems():
        augmented_map[key] = [AugLeveledArc(x,
                                            site_lat,
                                            site_lon,
                                            height=height) for x in arc_list]
    return augmented_map


def get_dt_list(arc_map):
    """
    """
    dt_set = set()
    for arc_list in arc_map.itervalues():
        for arc in arc_list:
            dt_set.update(arc.dt)
    return sorted(dt_set)


def ionex_stec_map(ionex_fname,
                   augmented_arc_map):
    """
    """
    logger.info('computing interpolated and mapped STEC from {}'.format(ionex_fname))
    # compute temporally interpolated IONEX VTEC maps
    dt_list = get_dt_list(augmented_arc_map)
    (grid_lon, grid_lat, vtec,
     _, sat_biases, _) = interpolate2D_temporal(ionex_fname,
                                                dt_list)
    # compute spline interpolators
    bbox = [-180, 180, -90, 90]
    # check that latitude grid is in decreasing order
    assert grid_lat[1] < grid_lat[0]
    interp_map = {dt: RectBivariateSpline(grid_lon,
                                          grid_lat[::-1],
                                          vtec[:, ::-1, i],
                                          bbox=bbox) for i, dt in enumerate(dt_list)}
    # compute interpolated stec
    stec_map = defaultdict(list)
    for key, arc_list in augmented_arc_map.iteritems():
        for arc in arc_list:
            ionex_stec = []
            for (dt_i, ipp_lat_i, ipp_lon_i, el_map_i) in zip(arc.dt,
                                                              arc.ipp_lat,
                                                              arc.ipp_lon,
                                                              arc.el_map):
                i, _ = find_le(dt_list, dt_i)
                interpolator = interp_map[dt_list[i]]
                vtec_i = float(interpolator.ev(ipp_lon_i,
                                               ipp_lat_i))
                ionex_stec.append(vtec_i * el_map_i)
            stec_map[key].append(ionex_stec)
    return stec_map, sat_biases


def estimate_receiver_bias(arc_map,
                           model_stec_map,
                           sat_biases,
                           n_std=3):
    """
    """
    logger.info('estimating receiver bias')
    # gather vectors
    stec_minus_sat_bias = []
    el = []
    model_sobs = []
    for sat, arcs in arc_map.iteritems():
        if not sat.startswith('G'):
            raise NotImplementedError('only GPS satellites are currently '
                                      'supported')
        # IONEX DCBs are given in [ns] --- convert to [TECU]
        sat_bias = -sat_biases['GPS'][int(sat[1:])][0] / TECU_TO_NS
        for i, arc in enumerate(arcs):
            stec_minus_sat_bias.extend([x - sat_bias for x in arc.stec])
            el.extend(arc.el)
            model_sobs.extend(model_stec_map[sat][i])
    # flag outliers
    srgim = NP.ma.masked_invalid(NP.array(model_sobs) -
                                 NP.array(stec_minus_sat_bias))
    mean = NP.mean(srgim)
    sigma = NP.std(srgim)
    I = NP.nonzero(NP.abs(srgim[~srgim.mask]) > abs(mean) + n_std * sigma)
    srgim[I] = NP.ma.masked
    bias = -NP.sum(srgim * el) / NP.sum(el)
    sigma = NP.std(srgim)
    return bias, sigma


"""
???
"""
JPLH_TEMPLATE = '/pub/iono_daily/IONEX_rapid/JPLH{date:%j}0.{date:%y}I.gz'


"""
???
"""
JPLH_ARCHIVE_TEMPLATE = '/pub/iono_daily/IONEX_rapid/archive/JPLH{date:%j}0.{date:%y}I.gz'


def fetch_sideshow_ionex(path,
                         date,
                         work_path=None,
                         templates=[JPLH_TEMPLATE, JPLH_ARCHIVE_TEMPLATE]):
    """
    ???
    """
    with SmartTempDir(work_path) as work_path:
        for template in templates:
            server_fname = template.format(date=date)
            local_fname = os.path.join(path, posixpath.basename(server_fname)[:-3])
            try:
                update_sideshow_file(local_fname,
                                     server_fname)
                return local_fname
            except:
                logger.info('could not download {}'.format(server_fname))
                continue
    raise RuntimeError('could not download IONEX file from sideshow for {:%Y-%m-%d}'.format(date))


def bias_process(output_h5_fname,
                 leveled_arc_h5_fname,
                 ionex_fname):
    """
    ???
    """
    arc_map = ArcMap(leveled_arc_h5_fname)
    return output_h5_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Estimate receiver bias from phase-leveled data.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_h5_fname',
                        type=str,
                        help='???')
    parser.add_argument('leveled_arc_h5_fname',
                        type=str,
                        help='input H5 file containing leveled phase arcs')
    parser.add_argument('--work-path',
                        '-w',
                        type=str,
                        default=None,
                        help='path to store intermediate files (use an '
                             'automatically cleaned up area if not specified)')
    calibration_group = parser.add_mutually_exclusive_group(required=True)
    calibration_group.add_argument('--ionex_fname',
                                   '-i',
                                   type=str,
                                   help='calibrate using the given IONEX file for satellite biases and ionospheric delay (fetch from JPL sideshow if not specified)')
    calibration_group.add_argument('--date',
                                   '-d',
                                   type=lambda x: datetime.strptime(x, '%Y-%m-%d'),
                                   help='fetch IONEX for the given date')
    args = parser.parse_args(argv[1:])


    with SmartTempDir(args.work_path) as work_path:
        if args.ionex_fname is None:
            ionex_fname = fetch_sideshow_ionex(work_path, args.date)
        else:
            ionex_fname = args.ionex_fname
        bias_process(args.output_h5_fname,
                     args.leveled_arc_h5_fname,
                     ionex_fname)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)
    sys.exit(main())


# if __name__ == '__main__':
#     from phase_edit import phase_edit, filter_obs_map
#     from rinex import read_rindump, dump_rinex
#     from level import level_phase_to_code

#     logging.basicConfig(level=logging.INFO)
#     logging.getLogger('sh').setLevel(logging.WARNING)
#     stn = 'artu'
#     rinex_fname = '/Users/butala/Documents/IARPA/absolute_tec/{}0010.14o'.format(stn)
#     nav_fname = '/Users/butala/Documents/IARPA/absolute_tec/{}0010.14n'.format(stn)
#     interval = 30

#     ionex_fname = '/Users/butala/Documents/IARPA/absolute_tec/JPLH0010.14I'

#     (time_reject_map,
#      phase_adjust_map) = phase_edit(rinex_fname, interval)

#     import os
#     rinex_dump_fname = os.path.join('/tmp', os.path.basename(rinex_fname) + '.dump')
#     dump_rinex(rinex_dump_fname,
#                rinex_fname,
#                nav_fname)

#     obs_map = read_rindump(rinex_dump_fname)

#     obs_map = filter_obs_map(obs_map,
#                              time_reject_map,
#                              phase_adjust_map)

#     arc_map = level_phase_to_code(obs_map)

#     info = rinex_info(rinex_fname,
#                       nav_fname)

#     aug_arc_map = augmented_arc_map(arc_map,
#                                     info['lat'],
#                                     info['lon'])

#     (stec_map,
#      sat_biases) = ionex_stec_map(ionex_fname,
#                                   aug_arc_map)

#     bias, sigma = estimate_receiver_bias(aug_arc_map,
#                                          stec_map,
#                                          sat_biases)

#     print('estimated receiver bias = {:.1f} [TECU]'.format(bias))
