import logging
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from itertools import repeat, izip
from multiprocessing import Pool, cpu_count
from collections import OrderedDict, namedtuple

from tables import open_file, IsDescription, Time64Col, Float64Col

from ..gpstk import PyPosition
from ..util.path import SmartTempDir, replace_path
from ..util.date import UNIX_EPOCH
from ..iri.iri_stec import iri_stec
from rinex import dump_preprocessed_rinex, read_rindump
from observation import ObsMapFlatIterator


class STecInfo(namedtuple('StecInfo',
                          'stec az el satx saty satz')):
    pass


class STecMap(dict):
    def __missing__(self, key):
        self[key] = OrderedDict()
        return self[key]


def iri_stec_wrap(x):
    """ ??? """
    stn_xyz, (dt, sat, obs) = x
    stn_pos = PyPosition(*stn_xyz)
    sat_pos = PyPosition(obs.satx,
                         obs.saty,
                         obs.satz)
    return (sat,
            dt,
            STecInfo(iri_stec(dt, stn_pos, sat_pos),
                     obs.az,
                     obs.el,
                     obs.satx,
                     obs.saty,
                     obs.satz))


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

        stec_map = STecMap()
        stec_map.xyz = obs_map.xyz
        stec_map.llh = obs_map.llh
        for sat, dt, stec_info in stec_output:
            stec_map[sat][dt] = stec_info
    return stec_map


class STecTable(IsDescription):
    dt   = Time64Col()
    stec = Float64Col()
    az   = Float64Col()
    el   = Float64Col()
    satx = Float64Col()
    saty = Float64Col()
    satz = Float64Col()


def dump_stec_map(h5_fname, stec_map):
    """ ??? """
    h5file = open_file(h5_fname, mode='w', title='IRI simulated slant TEC')
    group = h5file.create_group('/', 'phase_arcs', 'Phase connected arcs')
    if hasattr(stec_map, 'xyz'):
        group._v_attrs.xyz = stec_map.xyz
    if hasattr(stec_map, 'llh'):
        group._v_attrs.llh = stec_map.llh
    for sat in sorted(stec_map):
        assert sat[0] == 'G'
        table = h5file.create_table(group, sat, STecTable, 'GPS prn={} data'.format(sat[1:]))
        row = table.row
        for dt, stec_info in stec_map[sat].iteritems():
            row['dt'] = (dt - UNIX_EPOCH).total_seconds()
            row['stec'] = stec_info.stec
            row['az'] = stec_info.az
            row['el'] = stec_info.el
            row['satx'] = stec_info.satx
            row['saty'] = stec_info.saty
            row['satz'] = stec_info.satz
            row.append()
        table.flush()
    h5file.close()
    return h5_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Compute IRI slant TEC along lines of sight found in given RINEX file.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_h5_fname',
                        type=str,
                        help='output HDF5 file')
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

    dump_stec_map(args.output_h5_fname,
                  stec_map)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)
    sys.exit(main())
