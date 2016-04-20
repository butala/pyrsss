import logging
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from itertools import repeat, izip
from multiprocessing import Pool, cpu_count
from collections import OrderedDict, defaultdict
from cPickle import dump

from ..gpstk import PyPosition
from ..util.path import SmartTempDir, replace_path
from ..iri.iri_stec import iri_stec
from rinex import dump_preprocessed_rinex, read_rindump
from observation import ObsMapFlatIterator


def iri_stec_wrap(x):
    """ ??? """
    stn_xyz, (dt, sat, obs) = x
    stn_pos = PyPosition(*stn_xyz)
    sat_pos = PyPosition(obs.satx,
                         obs.saty,
                         obs.satz)
    return sat, dt, iri_stec(dt, stn_pos, sat_pos)


def rinex_iri_stec(obs_fname,
                   nav_fname,
                   work_path=None,
                   decimate=None,
                   processes=cpu_count()):
    """ ??? """
    with SmartTempDir(work_path) as work_path:
        dump_fname = replace_path(work_path, obs_fname + '.dump')
        dump_preprocessed_rinex(dump_fname,
                                obs_fname,
                                nav_fname,
                                work_path=work_path,
                                decimate=decimate)
        obs_map = read_rindump(dump_fname)
        obs_map_iter = ObsMapFlatIterator(obs_map)

        pool = Pool(processes)
        stec_output = pool.map(iri_stec_wrap, izip(repeat(obs_map.xyz), obs_map_iter))
        pool.close()
        pool.join()

        stec_map = defaultdict(OrderedDict)
        for sat, dt, stec in stec_output:
            stec_map[sat][dt] = stec
    return stec_map


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Compute IRI slant TEC along lines of sight found in given RINEX file.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_pkl_fname',
                        type=str,
                        help='output pickle file')
    parser.add_argument('rinex_obs_fname',
                        type=str,
                        help='input RINEX observation file')
    parser.add_argument('rinex_nav_fname',
                        type=str,
                        help='input RINEX navigation file')
    parser.add_argument('--work-path',
                        '-w',
                        type=str,
                        default=None,
                        help='path to store intermediate files (if not specified, use an automatically cleaned up temporary area)')
    parser.add_argument('--decimate',
                        '-d',
                        type=int,
                        default=None,
                        help='decimate to time interval in [s]')
    parser.add_argument('--processes',
                        '-p',
                        type=int,
                        default=cpu_count(),
                        help='use the given number of processes')
    args = parser.parse_args(argv[1:])

    stec_map = rinex_iri_stec(args.rinex_obs_fname,
                              args.rinex_nav_fname,
                              work_path=args.work_path,
                              decimate=args.decimate,
                              processes=args.processes)

    with open(args.output_pkl_fname, 'w') as fid:
        dump(stec_map, fid, -1)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)
    sys.exit(main())
