import sys
import logging
import os
import posixpath
from collections import namedtuple, defaultdict, OrderedDict
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import datetime, timedelta

import numpy as NP
from scipy.interpolate import RectBivariateSpline
from tables import open_file, IsDescription, Time64Col, Float64Col

from constants import SHELL_HEIGHT, TECU_TO_NS
from level import LeveledArc, ArcMap
from util import shell_mapping
from ipp import ipp_from_azel
from teqc import rinex_info
from sideshow import update_sideshow_file
from ..gpstk import PyPosition
from ..ionex.read_ionex import interpolate2D_temporal
from ..util.search import find_le
from ..util.path import SmartTempDir
from ..util.date import UNIX_EPOCH
from ..util.angle import convert_lon

logger = logging.getLogger('pyrsss.gps.bias')

"""
Remove transmitter and receiver biases from phase-connected arcs. Use
Attila's method to estimate the receiver bias (using IGS IONEX records
containing satellite biases and modeled VTEC maps).
"""


class AugLeveledArc(namedtuple('AugLeveledArc',
                               ' '.join(LeveledArc._fields + ('el_map',
                                                              'ipp_lat',
                                                              'ipp_lon'))),
                    LeveledArc):
    def __new__(cls,
                leveled_arc,
                stn_pos,
                shell_height=SHELL_HEIGHT):
        fields = leveled_arc._asdict()
        fields['el_map'] = [shell_mapping(el_i, h=shell_height) for el_i in fields['el']]
        ipp_pos = [ipp_from_azel(stn_pos, az_i, el_i) for az_i, el_i in zip(fields['az'],
                                                                            fields['el'])]
        fields['ipp_lat'] = [x.geodeticLatitude for x in ipp_pos]
        fields['ipp_lon'] = [convert_lon(x.longitude) for x in ipp_pos]
        return cls._make(fields.values())


class AugmentedArcMap(OrderedDict):
    def __init__(self, arc_map, shell_height=SHELL_HEIGHT):
        super(AugmentedArcMap, self).__init__()
        stn_pos = PyPosition(*arc_map.llh,
                             s=PyPosition.CoordinateSystem['geodetic'])
        for key, arc_list in arc_map.iteritems():
            self[key] = [AugLeveledArc(x,
                                       stn_pos,
                                       shell_height=shell_height) for x in arc_list]
        self.xyz = arc_map.xyz
        self.llh = arc_map.llh
        self.shell_height = shell_height


    def __missing__(self, key):
        """ ??? """
        self[key] = list()
        return self[key]


class CalibratedArc(namedtuple('CalibratedArc',
                               ' '.join(AugLeveledArc._fields).replace('stec',
                                                                       'sobs')),
                    AugLeveledArc):
    @classmethod
    def from_aug_leveled_arc(cls,
                             aug_leveled_arc,
                             sat_bias,
                             stn_bias):
        fields = aug_leveled_arc._asdict()
        fields['stec'] -= sat_bias + stn_bias
        return cls._make(fields.values())


class CalibratedArcMap(OrderedDict):
    def __missing__(self, key):
        """ ??? """
        self[key] = list()
        return self[key]

    @classmethod
    def from_aug_arc_map(cls,
                         aug_arc_map,
                         sat_biases,
                         stn_bias,
                         stn_bias_sigma):
        """ ??? """
        calibrated_arc_map = cls()
        for sat, aug_leveled_arcs in aug_arc_map.iteritems():
             assert sat.startswith('G')
             sat_bias = -sat_biases['GPS'][int(sat[1:])][0] / TECU_TO_NS
             calibrated_arc_map[sat] = [CalibratedArc.from_aug_leveled_arc(x, sat_bias, stn_bias) for x in aug_leveled_arcs]
        calibrated_arc_map.llh = aug_arc_map.llh
        calibrated_arc_map.xyz = aug_arc_map.xyz
        calibrated_arc_map.sat_biases = sat_biases
        calibrated_arc_map.stn_bias = stn_bias
        calibrated_arc_map.stn_bias_sigma = stn_bias_sigma
        return calibrated_arc_map

    """ ??? """
    class Table(IsDescription):
        dt   = Time64Col()
        sobs = Float64Col()
        sprn = Float64Col()
        az   = Float64Col()
        el   = Float64Col()
        satx = Float64Col()
        saty = Float64Col()
        satz = Float64Col()
        el_map = Float64Col()
        ipp_lat = Float64Col()
        ipp_lon = Float64Col()

    def dump(self, h5_fname):
        """ ??? """
        h5file = open_file(h5_fname, mode='w', title='pyrsss.gps.bias output')
        calibrated_phase_arcs_group = h5file.create_group('/',
                                                        'calibrated_phase_arcs',
                                                        'Calibrated phase connected arcs')
        calibrated_phase_arcs_group._v_attrs.xyz = self.xyz
        calibrated_phase_arcs_group._v_attrs.llh = self.llh
        calibrated_phase_arcs_group._v_attrs.stn_bias = self.stn_bias
        calibrated_phase_arcs_group._v_attrs.stn_bias_sigma = self.stn_bias_sigma
        for prn in sorted(self.sat_biases['GPS']):
            bias, rms = self.sat_biases['GPS'][prn]
            setattr(calibrated_phase_arcs_group._v_attrs, 'G{:02d}'.format(prn), -bias / TECU_TO_NS)
            setattr(calibrated_phase_arcs_group._v_attrs, 'G{:02d}_rms'.format(prn), rms / TECU_TO_NS)
        for sat in sorted(self):
            assert sat[0] == 'G'
            sat_group = h5file.create_group(calibrated_phase_arcs_group,
                                            sat,
                                            'Calibrated phase connected arcs for {}'.format(sat))
            for i, calibrated_arc in enumerate(self[sat]):
                table = h5file.create_table(sat_group,
                                            'arc' + str(i),
                                            CalibratedArcMap.Table,
                                            'GPS prn={} arc={} data'.format(sat, i))
                table.attrs.L = calibrated_arc.L
                table.attrs.L_scatter = calibrated_arc.L_scatter
                row = table.row
                for j in range(len(calibrated_arc.dt)):
                    row['dt'] = (calibrated_arc.dt[j] - UNIX_EPOCH).total_seconds()
                    row['sobs'] = calibrated_arc.sobs[j]
                    row['sprn'] = calibrated_arc.sprn[j]
                    row['az'] = calibrated_arc.az[j]
                    row['el'] = calibrated_arc.el[j]
                    row['satx'] = calibrated_arc.satx[j]
                    row['saty'] = calibrated_arc.saty[j]
                    row['satz'] = calibrated_arc.satz[j]
                    row['el_map'] = calibrated_arc.el_map[j]
                    row['ipp_lat'] = calibrated_arc.ipp_lat[j]
                    row['ipp_lon'] = calibrated_arc.ipp_lon[j]
                    row.append()
                table.flush()
        h5file.close()
        return h5_fname

    @classmethod
    def undump(cls, h5_fname):
        calibrated_arc_map = cls()
        h5file = open_file(h5_fname, mode='r')
        calibrated_phase_arcs_group = h5file.root.calibrated_phase_arcs
        for attr in calibrated_phase_arcs_group._v_attrs._f_list():
            setattr(calibrated_arc_map, attr, getattr(calibrated_phase_arcs_group._v_attrs, attr))
        for sat_group in calibrated_phase_arcs_group:
            sat = sat_group._v_name
            for arc_table in sat_group:
                dt = []
                sobs = []
                sprn = []
                az = []
                el = []
                satx = []
                saty = []
                satz = []
                el_map = []
                ipp_lat = []
                ipp_lon = []
                for row in arc_table.iterrows():
                    dt.append(UNIX_EPOCH + timedelta(seconds=row['dt']))
                    sobs.append(row['sobs'])
                    sprn.append(row['sprn'])
                    az.append(row['az'])
                    el.append(row['el'])
                    satx.append(row['satx'])
                    saty.append(row['saty'])
                    satz.append(row['satz'])
                    el_map.append(row['el_map'])
                    ipp_lat.append(row['ipp_lat'])
                    ipp_lon.append(row['ipp_lon'])
                calibrated_arc_map[sat].append(CalibratedArc(dt,
                                                         sobs,
                                                         sprn,
                                                         az,
                                                         el,
                                                         satx,
                                                         saty,
                                                         satz,
                                                         arc_table.attrs.L,
                                                         arc_table.attrs.L_scatter,
                                                         el_map,
                                                         ipp_lat,
                                                         ipp_lon))
        h5file.close()
        return calibrated_arc_map


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
    # load arc map
    arc_map = ArcMap(leveled_arc_h5_fname)
    # compute IPPs
    logger.info('computing IPPs')
    aug_arc_map = AugmentedArcMap(arc_map)
    # compute VTEC mapped to STEC for arc map lines of site
    logger.info('computing STEC from IONEX')
    (stec_map,
     sat_biases) = ionex_stec_map(ionex_fname,
                                  aug_arc_map)
    # estimate receiver bias and uncertainty
    logger.info('least squares estimate of receiver bias')
    stn_bias, stn_bias_sigma = estimate_receiver_bias(aug_arc_map,
                                                      stec_map,
                                                      sat_biases)
    # store output
    calibrated_arc_map = CalibratedArcMap.from_aug_arc_map(aug_arc_map,
                                                       sat_biases,
                                                       stn_bias,
                                                       stn_bias_sigma)
    calibrated_arc_map.dump(output_h5_fname)
    return output_h5_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Estimate receiver bias from phase-leveled data.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_h5_fname',
                        type=str,
                        help='output H5 files containing calibrated, leveled phase arcs')
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
