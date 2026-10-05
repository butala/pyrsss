"""
Remove transmitter and receiver biases from phase-connected arcs. Use
Attila's method to estimate the receiver bias (using IGS IONEX records
containing satellite biases and modeled VTEC maps) or apply the IONEX
DCB tables for both satellite and station biases directly.
"""
import logging
import os
import posixpath
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.interpolate import RectBivariateSpline

from ..ionex.read_ionex import interpolate2D_temporal
from ..util.angle import convert_lon
from ..util.path import SmartTempDir
from ..util.search import find_le
from .constants import SHELL_HEIGHT, TECU_TO_NS, NS_TO_TECU
from .level import LeveledArc
from .sideshow import update_sideshow_file
from .util import shell_mapping

logger = logging.getLogger('pyrsss.gps.bias')


class CalibratedArc(pd.DataFrame):
    """
    Absolutely calibrated arc: columns gps_time, az, el, satx/y/z, sobs,
    sprn ([TECU]) and, when augmented (Attila's method), el_map, ipp_lat,
    ipp_lon. Bias metadata is attached as attributes.
    """
    _metadata = ['xyz',
                 'llh',
                 'stn',
                 'recv_type',
                 'sat',
                 'L',
                 'L_scatter',
                 'sat_bias',
                 'stn_bias',
                 'stn_bias_sigma']

    @property
    def _constructor(self):
        return CalibratedArc


def augment_arc(arc, stn_pos=None, shell_height=SHELL_HEIGHT):
    """
    Return a copy of the :class:`LeveledArc` *arc* with the mapped shell
    height (*el_map*) and ionospheric pierce point (*ipp_lat*,
    *ipp_lon* [deg]) columns added. *stn_pos* is a :class:`PyPosition`
    built from *arc.llh* when not given.
    """
    # the gnsstk extension is only required for IPP computation
    from ..gnsstk import PyPosition
    from .ipp import ipp_from_azel

    if stn_pos is None:
        stn_pos = PyPosition(*arc.llh,
                             s=PyPosition.CoordinateSystem['geodetic'])
    ipp_pos = [ipp_from_azel(stn_pos, az_i, el_i)
               for az_i, el_i in zip(arc.az, arc.el)]
    aug = arc.copy()
    aug['el_map'] = [shell_mapping(el_i, h=shell_height) for el_i in arc.el]
    aug['ipp_lat'] = [x.geodeticLatitude for x in ipp_pos]
    aug['ipp_lon'] = [convert_lon(x.longitude) for x in ipp_pos]
    return aug


def ionex_stec_map(ionex_fname,
                   arcs):
    """
    Return the tuple (*model_stec*, *sat_biases*) where *model_stec* is a
    list (one entry per :class:`LeveledArc` of *arcs*) of IONEX-modeled
    STEC [TECU] along the line of sight and *sat_biases* is the IONEX
    DCB table.
    """
    logger.info('computing interpolated and mapped STEC from {}'.format(ionex_fname))
    # compute temporally interpolated IONEX VTEC maps
    dt_list = sorted({dt for arc in arcs for dt in arc.gps_time})
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
    model_stec = []
    for arc in arcs:
        ionex_stec = []
        for (dt_i, ipp_lat_i, ipp_lon_i, el_map_i) in zip(arc.gps_time,
                                                          arc.ipp_lat,
                                                          arc.ipp_lon,
                                                          arc.el_map):
            i, _ = find_le(dt_list, dt_i)
            interpolator = interp_map[dt_list[i]]
            vtec_i = float(interpolator.ev(ipp_lon_i,
                                           ipp_lat_i))
            ionex_stec.append(vtec_i * el_map_i)
        model_stec.append(ionex_stec)
    return model_stec, sat_biases


def estimate_receiver_bias(arcs,
                           model_stec,
                           sat_biases,
                           n_std=3):
    """
    Estimate the receiver bias [TECU] and its uncertainty given the
    augmented *arcs* and the IONEX-modeled STEC *model_stec* (see
    :func:`ionex_stec_map`). Return the (bias, sigma) tuple.
    """
    logger.info('estimating receiver bias')
    # gather vectors
    stec_minus_sat_bias = []
    el = []
    model_sobs = []
    for arc, model in zip(arcs, model_stec):
        sat = arc.sat
        if not sat.startswith('G'):
            raise NotImplementedError('only GPS satellites are currently '
                                      'supported')
        # IONEX DCBs are given in [ns] --- convert to [TECU]
        sat_bias = -sat_biases['GPS'][int(sat[1:])][0] / TECU_TO_NS
        stec_minus_sat_bias.extend([x - sat_bias for x in arc.L_I])
        el.extend(arc.el)
        model_sobs.extend(model)
    # flag outliers
    srgim = np.ma.masked_invalid(np.array(model_sobs) -
                                 np.array(stec_minus_sat_bias))
    mean = np.mean(srgim)
    sigma = np.std(srgim)
    I = np.nonzero(np.abs(srgim[~srgim.mask]) > abs(mean) + n_std * sigma)
    srgim[I] = np.ma.masked
    bias = -np.sum(srgim * el) / np.sum(el)
    sigma = np.std(srgim)
    return bias, sigma


def calibrate_arcs(arcs,
                   sat_biases,
                   stn_bias,
                   stn_bias_sigma=None):
    """
    Return the list of :class:`CalibratedArc` obtained by removing the
    IONEX satellite DCB and the estimated receiver bias *stn_bias*
    [TECU] from the leveled phase (L_I) of each arc of *arcs* (Attila's
    method).
    """
    calibrated_arcs = []
    for arc in arcs:
        assert arc.sat.startswith('G')
        # IONEX DCBs are given in [ns] --- convert to [TECU]
        sat_bias = -sat_biases['GPS'][int(arc.sat[1:])][0] / TECU_TO_NS
        data_map = {'gps_time': arc.gps_time.values,
                    'az': arc.az.values,
                    'el': arc.el.values,
                    'satx': arc.satx.values,
                    'saty': arc.saty.values,
                    'satz': arc.satz.values,
                    'sobs': arc.L_I.values - (sat_bias + stn_bias),
                    'sprn': arc.P_I.values}
        for column in ('el_map', 'ipp_lat', 'ipp_lon'):
            if column in arc.columns:
                data_map[column] = arc[column].values
        calibrated_arc = CalibratedArc(data_map)
        for attr in LeveledArc._metadata:
            setattr(calibrated_arc, attr, getattr(arc, attr))
        calibrated_arc.sat_bias = sat_bias
        calibrated_arc.stn_bias = stn_bias
        calibrated_arc.stn_bias_sigma = stn_bias_sigma
        calibrated_arcs.append(calibrated_arc)
    return calibrated_arcs


def calibrate_dcb(arcs,
                  sat_biases,
                  stn_biases):
    """
    Return the list of :class:`CalibratedArc` obtained by applying the
    IONEX DCB tables *sat_biases* and *stn_biases* to both the leveled
    phase (L_I) and code (P_I) of each arc of *arcs*.
    """
    calibrated_arcs = []
    for arc in arcs:
        if arc.sat[0] == 'G':
            sat_bias = sat_biases['GPS'][int(arc.sat[1:])][0] * NS_TO_TECU
            stn_bias = stn_biases['GPS'][arc.stn.upper()][0] * NS_TO_TECU
        elif arc.sat[0] == 'R':
            sat_bias = sat_biases['GLONASS'][int(arc.sat[1:])][0] * NS_TO_TECU
            stn_bias = stn_biases['GLONASS'][arc.stn.upper()][0] * NS_TO_TECU
        else:
            raise ValueError('Satellite bias for {} not found'.format(arc.sat))

        data_map = {'gps_time': arc.gps_time.values,
                    'az': arc.az.values,
                    'el': arc.el.values,
                    'satx': arc.satx.values,
                    'saty': arc.saty.values,
                    'satz': arc.satz.values,
                    'sobs': arc.L_I.values + sat_bias + stn_bias,
                    'sprn': arc.P_I.values + sat_bias + stn_bias}
        calibrated_arc = CalibratedArc(data_map)
        for attr in LeveledArc._metadata:
            setattr(calibrated_arc, attr, getattr(arc, attr))
        calibrated_arc.sat_bias = sat_bias
        calibrated_arc.stn_bias = stn_bias
        calibrated_arcs.append(calibrated_arc)
    return calibrated_arcs


JPLH_TEMPLATE = '/pub/iono_daily/IONEX_rapid/JPLH{date:%j}0.{date:%y}I.gz'
"""Sideshow path template for rapid JPL IONEX records."""


JPLH_ARCHIVE_TEMPLATE = '/pub/iono_daily/IONEX_rapid/archive/JPLH{date:%j}0.{date:%y}I.gz'
"""Sideshow path template for archived JPL IONEX records."""


def fetch_sideshow_ionex(path,
                         date,
                         work_path=None,
                         templates=[JPLH_TEMPLATE, JPLH_ARCHIVE_TEMPLATE]):
    """
    Fetch the JPL IONEX record for *date* into *path* from the sideshow
    FTP server. Return the local file name.
    """
    with SmartTempDir(work_path) as work_path:
        for template in templates:
            server_fname = template.format(date=date)
            local_fname = os.path.join(path, posixpath.basename(server_fname)[:-3])
            try:
                update_sideshow_file(local_fname,
                                     server_fname)
                return local_fname
            except Exception:
                logger.info('could not download {}'.format(server_fname))
                continue
    raise RuntimeError('could not download IONEX file from sideshow for {:%Y-%m-%d}'.format(date))


def bias_process(leveled_arcs,
                 ionex_fname):
    """
    Run Attila's method end-to-end: augment the *leveled_arcs* with IPP
    information, estimate the receiver bias against the IONEX record
    *ionex_fname*, and return the list of :class:`CalibratedArc`.
    """
    logger.info('computing IPPs')
    aug_arcs = [augment_arc(arc) for arc in leveled_arcs]
    logger.info('computing STEC from IONEX')
    (model_stec,
     sat_biases) = ionex_stec_map(ionex_fname,
                                  aug_arcs)
    logger.info('least squares estimate of receiver bias')
    stn_bias, stn_bias_sigma = estimate_receiver_bias(aug_arcs,
                                                      model_stec,
                                                      sat_biases)
    return calibrate_arcs(aug_arcs,
                          sat_biases,
                          stn_bias,
                          stn_bias_sigma)


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Estimate receiver bias from phase-leveled data.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_fname',
                        type=str,
                        help='output pickle file containing calibrated, leveled phase arcs')
    parser.add_argument('leveled_arcs_fname',
                        type=str,
                        help='input pickle file containing leveled phase arcs')
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
        leveled_arcs = pd.read_pickle(args.leveled_arcs_fname)
        calibrated_arcs = bias_process(leveled_arcs, ionex_fname)
        pd.to_pickle(calibrated_arcs, args.output_fname)
    return


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)
    sys.exit(main())
