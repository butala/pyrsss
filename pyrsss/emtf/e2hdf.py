import sys
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import numpy as np
import pandas as pd

from .calc_e import apply_transfer_function as tf_1D
from .calc_e_3d import apply_transfer_function as tf_3D
from ..mag.iaga2hdf import read_hdf, write_hdf
from ..usarray_emtf.index import get_index
from .usgs_regions import get_region


logger = logging.getLogger('pyrsss.emtf.e2hdf')


def find_1D(header):
    """
    Find and return the key associated with the Fernberg physiographic
    region at the location of the magnetometer as given in *header*.
    """
    lat = header['geodetic_latitude']
    lon = header['geodetic_longitude']
    region = get_region(lat, lon)
    if region:
        region = region.replace('-', '_')
    return region


def find_3D(index,
            header,
            max_distance,
            quality=5):
    """
    Return the list of USArray 3-D EMTF keys (for repository
    :class:`Index` *index*) that are less than *max_distance* (in km)
    from the magnetometer location specified in *header*. Only include
    *quality* or greater sites.
    """
    lat = header['geodetic_latitude']
    lon = header['geodetic_longitude']
    return index.quality_subset(min_quality=quality).by_distance(lat, lon, max_distance)


def apply_emtf(df_E,
               df_B,
               emtf_key,
               index,
               extrapolate0=True):
    """
    Apply the EMTF associated with *emtf_key* to magnetometer data
    found in *df_B* and store result to *df_E*. Use USArray .xml
    repository information :class:`Index` to process the 3-D EMTFs.
    """
    logger.info('applying transfer function {}'.format(emtf_key))
    interval = np.diff(df_B.index.values[:2])[0] / np.timedelta64(1, 's')
    Bx = df_B.B_X.values
    By = df_B.B_Y.values
    if emtf_key.startswith('USArray'):
        xml_fname = index[emtf_key][1]
        Ex, Ey = tf_3D(Bx, By, interval, xml_fname, extrapolate0=extrapolate0)
    else:
        Ex, Ey = tf_1D(Bx, By, interval, emtf_key)
    df_E[emtf_key + '_X'] = Ex
    df_E[emtf_key + '_Y'] = Ey
    return df_E


def e2hdf(hdf_fname,
          source_key='B',
          key='E',
          replace=False,
          include=[],
          exclude=[],
          _3D=None,
          _1D=False,
          quality=5,
          verbose=True):
    """
    Add modeled E columns to the :class:`DataFrame` record stored at
    *hdf_fname*. The input to the process (processed magnetometer B_X
    and B_Y) are found at *source_key* and the output (E_X and E_Y for
    various models) is stored associated with *key*. If *replace*,
    replace the *key* data record otherwise add to the data
    record. The following control which EMTFs are applied:

    - *include*: list of keys to always include

    - *exclude*: list of keys to never include

    - *_3D*: tuple of two values: 1) the maximum range (in km) from
             the measurement site to include USArray EMTF and 2) the
             path containing the repository of USArray EMTF .xml files

    - *_1D*: if True, include E generated by the 1-D Fernberg model
             for the physiographic region enclosing the magnetometer
             measurement point

    - *quality*: exclude USArray 3-D EMTFs with flagged at a quality
                 index less than the specified values (5 being the
                 highest)

    If *verbose*, then report the set of applied EMTFs to the logging
    facilities.
    """
    # setup target DataFrame
    df, header = read_hdf(hdf_fname, source_key)
    def empty_record():
        return pd.DataFrame(index=df.index)
    if replace:
        logger.info('creating new E record')
        df_e = empty_record()
    else:
        try:
            df_e, _ = read_hdf(hdf_fname, key)
            logger.info('appending to existing E record')
        except KeyError:
            logger.info('creating new E record')
            df_e = empty_record()
    # determine which EMTFs to use
    emtf_set = set(include) - set(exclude)
    if _1D:
        emtf_1D = find_1D(header)
        if emtf_1D is not None:
            emtf_set.add(emtf_1D)
    if _3D is not None:
        d_km, repository_path = _3D
        index = get_index(repository_path)
        for emtf_3D in find_3D(index, header, d_km):
            emtf_set.add(emtf_3D)
    else:
        index = None
    # apply EMTFs
    if verbose:
        logger.info('Applying EMTFs: {}'.format(', '.join(sorted(emtf_set))))
    for emtf_key in sorted(emtf_set):
        df_e = apply_emtf(df_e, df, emtf_key, index)
    # output DataFrame
    write_hdf(hdf_fname, df_e, key, header)
    return hdf_fname



def e2hdf_3D(hdf_fname,
             keys_3D,
             repository_path,
             source_key='B',
             key='E',
             replace=False):
    """
    Add modeled E columns to the :class:`DataFrame` record stored at
    *hdf_fname*. The input to the process (processed magnetometer B_X
    and B_Y) are found at *source_key* and the output (E_X and E_Y for
    various models) is stored associated with *key*. If *replace*,
    replace the *key* data record otherwise add to the data
    record. Apply the 3-D transfer functions identified by the list
    *keys_3D*. Search for EMTF data records at *repository_path*.
    """
    # setup target DataFrame
    df, header = read_hdf(hdf_fname, source_key)
    def empty_record():
        return pd.DataFrame(index=df.index)
    if replace:
        logger.info('creating new E record')
        df_e = empty_record()
    else:
        try:
            df_e, _ = read_hdf(hdf_fname, key)
            logger.info('appending to existing E record')
        except KeyError:
            logger.info('creating new E record')
            df_e = empty_record()
    # determine which EMTFs to use
    emtf_list = []
    index = get_index(repository_path)
    for key_3D in keys_3D:
        candidates = [x for x in index if key_3D in x]
        if len(candidates) == 1:
            emtf_list.append(candidates[0])
        elif len(candidates) == 0:
            raise KeyError('could not find {} in index.pkl at {}'.format(key_3D,
                                                                         repository_path))
        else:
            raise KeyError('could not unambiguously resolve {} in index.pkl at {} ({} are all candidates)'.format(key_3D,
                                                                                                                  repository_path,
                                                                                                                  ', '.join(candidates)))
    # apply EMTFs
    logger.info('Applying EMTFs: {}'.format(', '.join(sorted(emtf_list))))
    for emtf_key in emtf_list:
        df_e = apply_emtf(df_e, df, emtf_key, index)
    # output DataFrame
    write_hdf(hdf_fname, df_e, key, header)
    return hdf_fname


def float_or_str(x):
    """
    Return *x* converted to float or *x* if that fails.
    """
    try:
        return float(x)
    except:
        return x


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Add modeled E field records to HDF file containing processed B field data.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('hdf_fnames',
                        type=str,
                        nargs='*',
                        metavar='hdf_fname',
                        help='HDF file record to process')
    parser.add_argument('--source-key',
                        '-s',
                        type=str,
                        default='B',
                        help='')
    parser.add_argument('--key',
                        '-k',
                        type=str,
                        default='E',
                        help='key to associate with the processed records')
    parser.add_argument('--replace',
                        '-r',
                        action='store_true',
                        help='replace modeled E field record (otherwise, append to the existing record)')
    parser.add_argument('--include',
                        '-i',
                        nargs='+',
                        type=str,
                        default=[],
                        help='EMTFs to include')
    parser.add_argument('--exclude',
                        '-e',
                        nargs='+',
                        type=str,
                        default=[],
                        help='EMTFs to exclude')
    parser.add_argument('--1D',
                        action='store_true',
                        help='include Fernberg 1-D model result for the physiographic region at the measurement location (if interior to a physiographic region)')
    parser.add_argument('--3D',
                        type=float_or_str,
                        nargs=2,
                        help='two arguments: 1) include USArray 3-D model results within the specified geodetic distance (in km) from the measurement location and 2) the path to the USArray EMTF .xml repository')
    parser.add_argument('--quality',
                        '-q',
                        choices=range(6),
                        type=int,
                        default=5,
                        help='minimum acceptable quality USArray EMTF (i.e., 0 means use all and 5 means use only highest flagged transfer functions)')
    args = parser.parse_args(argv[1:])

    for hdf_fname in args.hdf_fnames:
        e2hdf(hdf_fname,
              source_key=args.source_key,
              key=args.key,
              replace=args.replace,
              include=args.include,
              exclude=args.exclude,
              _3D=getattr(args, '3D'),
              _1D=getattr(args, '1D'),
              quality=args.quality)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
