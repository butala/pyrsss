import os
import sys
import shutil
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import sh

from .. util.path import replace_path

logger = logging.getLogger('pyrsss.gps.preprocess')


def normalize_rinex(output_rinex_fname,
                    rinex_fname,
                    gps=True,
                    glonass=False,
                    sbas=False,
                    galileo=False,
                    compass=False,
                    qzss=False,
                    decimate=None,
                    clockprep=True,
                    obs_types = 'L1L2P1P2C1'):
    """
    Sanitize *rinex_fname* with teqc (which is expected to be found in
    the path) and store output in *output_rinex_fname*. Satellite
    system data are excluded if *gps*, *glonass*, *sbas*, *galileo*,
    *compass*, or *qzss* are `False`. If *decimate* is given, decimate
    the output to the given interval (in [s]). If *clockprep*, smooth
    time tags (useful for receivers with steered clocks). Only output
    the observations givne in *obs_types*. Return
    *output_rinex_fname*.
    """
    if not os.path.isfile(rinex_fname):
        raise ValueError('RINEX observation file {} not found'.format(rinex_fname))
    args = ['+C2']
    args += ['+G'] if gps else ['-G']
    args += ['+R'] if glonass else ['-R']
    args += ['+S'] if sbas else ['-S']
    args += ['+E'] if galileo else ['-E']
    args += ['+C'] if compass else ['-C']
    args += ['+J'] if qzss else ['-J']
    args += ['+smtt'] if clockprep else ['-smtt']
    args += ['-O.obs_types', obs_types]
    if decimate is not None:
        args += ['-O.dec', str(decimate)]
    args += [rinex_fname]
    logger.info('calling teqc with arguments: {}'.format(' '.join(args)))
    sh.teqc(args,
            _out=output_rinex_fname,
            _err=sys.stderr)
    return output_rinex_fname


def normalize_to_path(path,
                      rinex_fnames,
                      **kwds):
    """ ??? """
    output_rinex_fnames = []
    for rinex_fname in rinex_fnames:
        output_rinex_fnames.append(replace_path(path, rinex_fname))
        normalize_rinex(output_rinex_fnames[-1],
                        rinex_fname,
                        **kwds)
    return output_rinex_fnames


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Normalize a given RINEX file.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('path',
                        type=str,
                        help='output normalized RINEX files to path')
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
    parser.add_argument('--decimate',
                        '-d',
                        type=int,
                        default=None,
                        help='decimate output to given interval in [s] (and adjust header accordingly)')
    args = parser.parse_args(argv[1:])

    normalize_to_path(args.path,
                      args.rinex_fnames,
                      gps=not args.gps,
                      glonass=not args.glonass,
                      sbas=not args.sbas,
                      galileo=not args.galileo,
                      compass=not args.compass,
                      qzss=not args.qzss,
                      decimate=args.decimate)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)
    sys.exit(main())
