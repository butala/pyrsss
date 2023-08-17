import logging

import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from constants import NS_TO_TECU
from preprocess import normalize_rinex
from rinex import dump_rinex, get_receiver_position
from phase_edit_new import phase_edit, apply_phase_adjustments, apply_rejections, label_phase_arcs, parse_discfix_log
from rinex_new import week_sec2dt, RinexDump
from level_new import level
from ..ionex.read_ionex import read_header


class CalibratedArc(pd.DataFrame):
    _metadata = ['xyz',
                 'llh',
                 'stn',
                 'recv_type',
                 'stn',
                 'sat',
                 'L',
                 'L_scatter',
                 'stn_bias',
                 'sat_bias']

    @property
    def _constructor(self):
        return CalibratedArc


def calibrate(leveled_arcs, sat_biases, stn_biases):
    """
    ???
    """
    calibrated_arcs = []
    for arc in leveled_arcs:
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
                    'sobs': arc.L_I + sat_bias + stn_bias,
                    'sprn': arc.P_I + sat_bias + stn_bias}
        calibrated_arc = CalibratedArc(data_map)
        calibrated_arc.xyz = arc.xyz
        calibrated_arc.llh = arc.llh
        calibrated_arc.stn = arc.stn
        calibrated_arc.recv_type = arc.recv_type
        calibrated_arc.sat = arc.sat
        calibrated_arc.L = arc.L
        calibrated_arc.L_scatter = arc.L_scatter
        calibrated_arc.sat_bias = sat_bias
        calibrated_arc.stn_bias = stn_bias
        calibrated_arcs.append(calibrated_arc)
    return calibrated_arcs


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    FRONTEND = True

    # pkl_fname = '/tmp/jplm0010.14o.gps.pkl'
    pkl_fname = '/tmp/algo0010.14o.gps.pkl'

    ionex_jpl_fname = '/Users/butala/Documents/IARPA/absolute_tec/JPLH0010.14I'
    ionex_code_fname = '/Users/butala/Documents/IARPA/absolute_tec/codg0010.14i'

    if FRONTEND:
        # preprocess
        # rinex_fname = '/Users/butala/Documents/IARPA/absolute_tec/jplm0010.14o'
        rinex_fname = '/Users/butala/Documents/IARPA/absolute_tec/algo0010.14o'
        # normalized_rinex_fname = '/tmp/jplm0010.14o'
        normalized_rinex_fname = '/tmp/algo0010.14o'

        normalize_rinex(normalized_rinex_fname,
                        rinex_fname)

        # dump GPS
        # dump_fname = '/tmp/jplm0010.14o.gps.dump'
        dump_fname = '/tmp/algo0010.14o.gps.dump'
        # nav_fname = '/Users/butala/Documents/IARPA/absolute_tec/jplm0010.14n'
        nav_fname = '/Users/butala/Documents/IARPA/absolute_tec/algo0010.14n'

        receiver_position = get_receiver_position(normalized_rinex_fname,
                                                  nav_fname)

        dump_rinex(dump_fname,
                   normalized_rinex_fname,
                   nav_fname,
                   receiver_position=receiver_position)

        # apply GPSTk phase editor
        (time_reject_map,
         phase_adjust_map) = phase_edit(normalized_rinex_fname,
                                        work_path='/tmp/work')

        rinex_dump = RinexDump.load(dump_fname)
        # rinex_dump = RinexDump.load(dump_fname, p1c1=False)

        # log_fname = '/tmp/work/jplm0010.14o.df.log'
        log_fname = '/tmp/work/algo0010.14o.df.log'
        phase_breaks = parse_discfix_log(log_fname)
        label_phase_arcs(rinex_dump, phase_breaks)

        apply_phase_adjustments(rinex_dump, phase_adjust_map)
        apply_rejections(rinex_dump, time_reject_map)

        rinex_dump.to_pickle(pkl_fname)
    else:
        rinex_dump = pd.read_pickle(pkl_fname)

    leveled_arcs = level(rinex_dump)

    with open(ionex_jpl_fname) as fid:
        _, sat_biases_jpl, stn_biases_jpl = read_header(fid)

    with open(ionex_code_fname) as fid:
        _, sat_biases_code, stn_biases_code = read_header(fid)

    # print(leveled_arcs[0].columns)

    # print(stn_biases)
    # print(sat_biases)

    calibrated_arcs_jpl = calibrate(leveled_arcs, sat_biases_jpl, stn_biases_jpl)
    calibrated_arcs_code = calibrate(leveled_arcs, sat_biases_code, stn_biases_code)

    fig = plt.figure(figsize=(11, 8.5))
    for arc in calibrated_arcs_jpl:
        plt.plot_date(arc.gps_time.values,
                      arc.sobs.values,
                      marker='.',
                      ls='-',
                      c='b')
    for arc in calibrated_arcs_code:
        plt.plot_date(arc.gps_time.values,
                      arc.sobs.values,
                      marker='.',
                      ls='-',
                      c='r')
    plt.savefig('/tmp/gps_test.pdf',
                bbox_inches='tight')
