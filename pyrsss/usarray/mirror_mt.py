import os
import sys
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict

from iris import fetch
from info import info_map
from ..util.path import touch_path

logger = logging.getLogger('pyrsss.usarray.mirror_mt')


def mirror(path):
    """
    Fetch all USArray MT data and store records, one per station, at
    *path*. Return the mapping between station identifiers and the
    file name of the associated stored data record.
    """
    touch_path(path)
    info = info_map()
    file_map = OrderedDict()
    for row in info.itertuples():
        stn = row.Index
        start = row.start.to_pydatetime()
        stop = row.end.to_pydatetime()
        hdf_fname = os.path.join(path, '{}.hdf'.format(stn))
        if os.path.isfile(hdf_fname):
            logger.info('Data record {} already exists --- skipping'.format(hdf_fname))
            continue
        logger.info('Fetch data for {}, {} --- {}'.format(stn, start, stop))
        try:
            df = fetch(stn, start, stop)
        except Exception as e:
            logger.error('Error fetching data ({}) --- skipping'.format(e))
        logger.info('Writing record to {}'.format(hdf_fname))
        df.to_hdf(hdf_fname, 'level0')
        file_map[stn] = hdf_fname
        break
    return file_map


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Mirror the USArray MT data repository.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('path',
                        type=str,
                        help='location to store data records')
    args = parser.parse_args(argv[1:])

    mirror(args.path)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
