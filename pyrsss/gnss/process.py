import logging
import os
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import pandas as pd

from ..util.path import SmartTempDir, replace_path
from .bias import bias_process, fetch_sideshow_ionex
from .level import DEFAULT_CONFIG, level, parse_override
from .phase_edit import phase_edit_rinex
from .rinex import fname2date

logger = logging.getLogger('pyrsss.gps.process')


def process(path,
            rinex_fnames,
            nav_fname,
            work_path=None,
            discfix_args=[],
            leveling_config_overrides=[],
            ionex_fname=None):
    """
    End-to-end processing of the RINEX files *rinex_fnames* (with the
    navigation file *nav_fname*) to absolutely calibrated phase arcs
    written as pickles under *path*. Return the list of output file
    names.
    """
    calibrated_fnames = []
    ionex_map = {}
    config = parse_override(leveling_config_overrides, DEFAULT_CONFIG)
    with SmartTempDir(work_path) as work_path:
        for rinex_fname in rinex_fnames:
            # phase edit + observable dump
            logger.info('editing {}'.format(rinex_fname))
            try:
                rinex_dump = phase_edit_rinex(rinex_fname,
                                              nav_fname,
                                              work_path=work_path,
                                              discfix_args=discfix_args)
            except Exception as e:
                logger.warning('phase edit step failed for {} ({}) --- '
                               'skipping'.format(rinex_fname, e))
                continue
            # level phase to code
            logger.info('leveling {}'.format(rinex_fname))
            try:
                leveled_arcs = level(rinex_dump, config=config)
            except Exception as e:
                logger.warning('level step failed for {} ({}) --- '
                               'skipping'.format(rinex_fname, e))
                continue
            # receiver bias estimation and subtraction
            date = fname2date(rinex_fname)
            if ionex_fname:
                ionex_fname_date = ionex_fname
            else:
                if date not in ionex_map:
                    logger.info('fetching IONEX for {:%Y-%m-%d}'.format(date))
                    ionex_map[date] = fetch_sideshow_ionex(work_path, date)
                ionex_fname_date = ionex_map[date]
            logger.info('calibrating {}'.format(rinex_fname))
            try:
                calibrated_arcs = bias_process(leveled_arcs,
                                               ionex_fname_date)
            except Exception as e:
                logger.warning('bias calibration step failed for {} ({}) --- '
                               'skipping'.format(rinex_fname, e))
                continue
            output_fname = replace_path(path, rinex_fname + '.pkl')
            pd.to_pickle(calibrated_arcs, output_fname)
            calibrated_fnames.append(output_fname)
        return calibrated_fnames


def add_dashes(s):
    """
    Prefix each token of *s* with dashes suitable for the DiscFix
    command line (single character tokens get one dash).
    """
    output = []
    for x in s.split():
        if len(x) == 1:
            output.append('-' + x)
        else:
            output.append('--' + x)
    return output


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('End-to-end processing of RINEX to absolutely calibrated arcs.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('path',
                        type=str,
                        help='output path')
    parser.add_argument('nav_fname',
                        type=str,
                        help='input RINEX navigation file')
    parser.add_argument('rinex_fnames',
                        type=str,
                        nargs='+',
                        metavar='rinex_fname',
                        help='input RINEX file')
    parser.add_argument('--work-path',
                        '-w',
                        type=str,
                        default=None,
                        help='path to store intermediate files (use an '
                             'automatically cleaned up area if not specified)')
    parser.add_argument('--discfix-options',
                        '-d',
                        type=add_dashes,
                        default=[],
                        help='options to pass to the GNSSTk discontinuity fixer (see help message for pyrsss.gps.phase_edit) --- do not include the dashes')
    parser.add_argument('--leveling-config-overrides',
                        '-l',
                        metavar='leveling_config_override',
                        type=str,
                        nargs='+',
                        default=[],
                        help='overrides to default leveling configuration (see the help message for pyrsss.gps.level for the possibilities)')
    parser.add_argument('--ionex-fname',
                        '-i',
                        type=str,
                        default=None,
                        help='use the specified IONEX record for satellite biases and VTEC (if not specified, download automatically from JPL sideshow)')
    args = parser.parse_args(argv[1:])

    process(args.path,
            args.rinex_fnames,
            args.nav_fname,
            work_path=args.work_path,
            discfix_args=args.discfix_options,
            leveling_config_overrides=args.leveling_config_overrides,
            ionex_fname=args.ionex_fname)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)
    sys.exit(main())
