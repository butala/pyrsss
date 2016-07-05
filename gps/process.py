import logging
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from phase_edit import phase_edit_process
from level import level_process
from bias import fetch_sideshow_ionex, bias_process
from rinex import fname2date
from ..util.path import SmartTempDir, replace_path

logger = logging.getLogger('pyrsss.gps.process')


def process(path,
            rinex_fnames,
            nav_fname,
            work_path=None,
            discfix_args=[],
            leveling_config_overrides=[],
            ionex_fname=None):
    """
    ???
    """
    with SmartTempDir(work_path) as work_path:
        # phase edit
        phase_edit_h5 = []
        for rinex_fname in rinex_fnames:
            logger.info('editing {}'.format(rinex_fname))
            try:
                phase_edit_h5.append(
                    phase_edit_process(replace_path(work_path,
                                                    rinex_fname + '.phase_edit.h5'),
                                       rinex_fname,
                                       nav_fname,
                                       work_path=work_path,
                                       discfix_args=discfix_args))
            except Exception as e:
                logger.warning('phase edit step failed for {} ({}) --- '
                               'skipping'.format(rinex_fname, e))
                continue
        # level phase to code
        level_h5 = []
        for phase_edit_h5_i in phase_edit_h5:
            logger.info('leveling {}'.format(phase_edit_h5_i))
            try:
                level_h5.append(
                    level_process(replace_path(work_path,
                                               rinex_fname + '.level.h5'),
                                  phase_edit_h5_i,
                                  config_overrides=leveling_config_overrides))
            except Exception as e:
                logger.warning('level step failed for {} ({}) --- '
                               'skipping'.format(phase_edit_h5_i, e))
                continue
        # receiver bias estimation and subtraction
        calibrated_h5 = []
        ionex_map = {}
        for level_h5_i, rinex_fname in zip(level_h5, rinex_fnames):
            date = fname2date(rinex_fname)
            if ionex_fname:
                ionex_fname_date = ionex_fname
            else:
                if date not in ionex_map:
                    logger.info('fetching IONEX for {:%Y-%m-%d}'.format(date))
                    ionex_map[date] = fetch_sideshow_ionex(work_path, date)
                ionex_fname_date = ionex_map[date]
            logger.info('calibrating {}'.format(level_h5_i))
            try:
                calibrated_h5.append(
                    bias_process(replace_path(path,
                                              rinex_fname + '.h5'),
                                 level_h5_i,
                                 ionex_fname_date))
            except Exception as e:
                logger.warning('bias calibration step failed for {} ({}) --- '
                               'skipping'.format(level_h5_i, e))
                continue
        return calibrated_h5


def add_dashes(s):
    """
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
                        metavar='rinex_fanme',
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
                        help='options to pass to the GPSTk discontinuity fixer (see help message for pyrsss.gps.phase_edit) --- do not include the dashes')
    parser.add_argument('--leveling-config-overrides',
                        '-l',
                        metavar='leveling_config_override',
                        type=str,
                        nargs='+',
                        default=[],
                        help='overrides to default leveling configuration (see the help message fro pyrsss.gps.level for the possibilities)')
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
