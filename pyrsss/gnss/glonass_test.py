import logging

import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from .preprocess import normalize_rinex
from .rinex import dump_rinex, get_receiver_position
from .phase_edit_new import phase_edit, apply_phase_adjustments, apply_rejections, label_phase_arcs, parse_discfix_log
from .level_new import level
from .rinex_new import RinexDump
from ..ionex.read_ionex import read_header

from .gps_test import calibrate



if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    """
    # preprocess
    rinex_fname = '/Users/butala/Documents/IARPA/absolute_tec/jplm0010.14o'
    normalized_rinex_fname = '/tmp/jplm0010.14o'

    normalize_rinex(normalized_rinex_fname,
                    rinex_fname,
                    gps=True,
                    glonass=True)

    # dump GPS
    dump_fname = '/tmp/jplm0010.14o.gps.dump'
    nav_fname = '/Users/butala/Documents/IARPA/absolute_tec/jplm0010.14n'

    receiver_position = get_receiver_position(normalized_rinex_fname,
                                              nav_fname)

    dump_rinex(dump_fname,
               normalized_rinex_fname,
               nav_fname,
               receiver_position=receiver_position)

    # dump GLONASS
    dump_fname = '/tmp/jplm0010.14o.glo.dump'
    gnav_fname = '/Users/butala/Documents/IARPA/absolute_tec/jplm0010.14g'
    data_keys = ['RC1C', 'RC1P', 'RC2P', 'RL1C', 'RL2C']

    dump_rinex(dump_fname,
               normalized_rinex_fname,
               gnav_fname,
               data_keys=data_keys,
               receiver_position=receiver_position)

    # apply GNSSTk phase editor
    (time_reject_map,
     phase_adjust_map) = phase_edit(normalized_rinex_fname,
                                    work_path='/tmp/work',
                                    glonass=True)
    """

    FRONTEND = False

    # pkl_fname = '/tmp/jplm0010.14o.gps.pkl'
    pkl_fname = '/tmp/algo0010.14o.gps.pkl'

    ionex_jpl_fname = '/Users/butala/Documents/IARPA/absolute_tec/JPLH0010.14I'
    ionex_code_fname = '/Users/butala/Documents/IARPA/absolute_tec/codg0010.14i'

    if FRONTEND:
        # preprocess
        rinex_fname = '/Users/butala/Documents/IARPA/absolute_tec/algo0010.14o'
        normalized_rinex_fname = '/tmp/algo0010.14o'

        normalize_rinex(normalized_rinex_fname,
                        rinex_fname,
                        gps=True,
                        glonass=True)

        # dump GLONASS
        dump_fname = '/tmp/algo0010.14o.glonass.dump'
        nav_fname = '/Users/butala/Documents/IARPA/absolute_tec/algo0010.14n'
        gnav_fname = '/Users/butala/Documents/IARPA/absolute_tec/algo0010.14g'

        receiver_position = get_receiver_position(normalized_rinex_fname,
                                                  nav_fname)

        data_keys = ['RC1C', 'RC1P', 'RC2P', 'RL1C', 'RL2C', 'AZI', 'ELE', 'SVX', 'SVY', 'SVZ']

        dump_rinex(dump_fname,
                   normalized_rinex_fname,
                   gnav_fname,
                   #nav_fname,
                   data_keys=data_keys,
                   receiver_position=receiver_position)

        # apply GNSSTk phase editor
        (time_reject_map,
         phase_adjust_map) = phase_edit(normalized_rinex_fname,
                                        work_path='/tmp/work',
                                        glonass=True)

        rinex_dump = RinexDump.load(dump_fname, p1c1=False)

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

    with open(ionex_code_fname) as fid:
        _, sat_biases_code, stn_biases_code = read_header(fid)

    # print(leveled_arcs[0].columns)

    # print(stn_biases)
    # print(sat_biases)

    calibrated_arcs_code = calibrate(leveled_arcs, sat_biases_code, stn_biases_code)

    fig = plt.figure(figsize=(11, 8.5))
    for arc in calibrated_arcs_code:
        plt.plot_date(arc.gps_time.values,
                      arc.sobs.values,
                      marker='.',
                      ls='-',
                      c='r')
    plt.savefig('/tmp/glonass_test.pdf',
                bbox_inches='tight')
