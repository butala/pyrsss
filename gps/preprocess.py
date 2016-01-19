import os
import sys
import shutil
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import sh

from .. util.path import SmartTempDir

LOGGER = logging.getLogger('pyrsss.gps.preprocess')


def normalize_rinex(rinex_fname,
                    gps=True,
                    glonass=False,
                    sbas=False,
                    galileo=False,
                    compass=False,
                    qzss=False,
                    clockprep=True,
                    obs_types = 'L1L2P1P2C1',
                    work_path=None):
    """
    Sanitize *rinex_fname* (warning: in place) with teqc (which is
    expected to be found in the path). Satellite system data are
    excluded if *gps*, *glonass*, *sbas*, *galileo*, *compass*, or
    *qzss* are `False`. If *clockprep*, smooth time tags (useful for
    receivers with steered clocks). Only output the observations givne
    in *obs_types*. Store intermediate files in *work_path*. Return
    *rinex_fname*.
    """
    with SmartTempDir(work_path) as work_path:
        rinex_copy_fname = os.path.join(work_path, os.path.basename(rinex_fname))
        LOGGER.info('copying {} to {}'.format(rinex_fname,
                                              rinex_copy_fname))
        shutil.copyfile(rinex_fname, rinex_copy_fname)
        args = ['+C2']
        args += ['+G'] if gps else ['-G']
        args += ['+R'] if glonass else ['-R']
        args += ['+S'] if sbas else ['-S']
        args += ['+E'] if galileo else ['-E']
        args += ['+C'] if compass else ['-C']
        args += ['+J'] if qzss else ['-J']
        args += ['+smtt'] if clockprep else ['-smtt']
        args += ['-O.obs_types', obs_types]
        args += [rinex_copy_fname]
        LOGGER.info('calling teqc with arguments: {}'.format(' '.join(args)))
        sh.teqc(args,
                _out=rinex_fname,
                _err=sys.stderr)
    return rinex_fname


def main(args):
    parser = ArgumentParser('Normalize a given RINEX file.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('rinex_fnames',
                        nargs='+',
                        type=str,
                        metavar='rinex_fname',
                        help='RINEX file to normalize')
    parser.add_argument('--gps',
                        '-g',
                        action='store_true',
                        help='exclude GPS measurements')
    parser.add_argument('--glonass',
                        '-r',
                        action='store_true',
                        help='exclude GLONASS measurements')
    parser.add_argument('--sbas',
                        '-s',
                        action='store_true',
                        help='exclude SBAS measurements')
    parser.add_argument('--galileo',
                        '-e',
                        action='store_true',
                        help='exclude Galileo measurements')
    parser.add_argument('--compass',
                        '-c',
                        action='store_true',
                        help='exclude Compass measurements')
    parser.add_argument('--qzss',
                        '-j',
                        action='store_true',
                        help='exclude QZSS measurements')
    parser.add_argument('--work-path',
                        '-w',
                        type=str,
                        default=None,
                        help='location to store intermediate files')
    args = parser.parse_args(args)

    for rinex_fname in args.rinex_fnames:
        normalize_rinex(rinex_fname,
                        gps=not args.gps,
                        glonass=not args.glonass,
                        sbas=not args.sbas,
                        galileo=not args.galileo,
                        compass=not args.compass,
                        qzss=not args.qzss,
                        work_path=args.work_path)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)
    main(sys.argv[1:])
