import sys
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict

import numpy as np
import pandas as pd

from ..usarray_emtf.index import get_index
from .calc_e_3d import apply_transfer_function

logger = logging.getLogger('pyrsss.emtf.e3d2hdf')


"""
Mappings between coordinate frame identifiers and northward and eastward measurement suffixes.
"""
COORDINATE_FRAME_SUFFIX = {'geo': {'north': 'X',
                                   'east':  'Y'},
                           'mag': {'north': 'N',
                                   'east':  'E'}}


def coordinate_frame_help(coord_map):
    """
    Return a help message explaining the contents of the coordinate
    frame measurement suffix mapping *coord_map*.
    """
    message = 'mappings between coordinate systems and measurement suffixes: '
    for frame, mapping in coord_map.items():
        message += ' [{}: north -> _{}, east -> {}]'.format(frame,
                                                            mapping['north'],
                                                            mapping['east'])
    return message


def get_xml_map(repository_path, station_names):
    """
    Find XML 3-D EMTF record for each identifier in *station_names* at
    *repository_path*. Return a mapping been the identifier and the
    tuple containing related information and XML file name.
    """
    index = get_index(repository_path)
    xml_map = OrderedDict()
    for station_name in station_names:
        candidates = [x for x in index.keys() if station_name in x]
        if len(candidates) == 1:
            xml_map[station_name] = index[candidates[0]]
        elif len(candidates) == 0:
            raise ValueError('no XML file associated with {} found at {}'.format(station_name,
                                                                                 repository_path))
        else:
            raise ValueError('multiple XML files associated with {} found at {}: {}'.format(station_name,
                                                                                            repository_path,
                                                                                            candidates))
    return xml_map


def e3d2hdf(hdf_fname,
            station_names,
            input_key='B',
            N_suffix='_X',
            E_suffix='_Y',
            output_hdf_fname=None,
            repository_path='.'):
    """
    Compute the 3-D EMTF modeled E-field for the B-field data found in
    *hdf_fname* for each station identifier given in
    *station_names*. Use the B-field data associated with record
    *input_key*. Write the data to a data frame record with key
    prefixed by *output_key_prefix*. Northward data (B and E) use
    suffix *N_suffix* and eastward data use suffix *E_suffix*. Output
    records to *output_hdf_fname* (and use *hdf_fname* by
    default). Use *repository_path* to find XML 3-D EMTF
    records. Return *output_hdf_fname*.
    """
    df_B = pd.read_hdf(hdf_fname, input_key)
    if output_hdf_fname is None:
        output_hdf_fname = hdf_fname
    interval = np.diff(df_B.index.values[:2])[0] / np.timedelta64(1, 's')
    logger.info('coordinate frame: north suffix = {}, east suffix = {}'.format(N_suffix, E_suffix))
    for station_name, (xml_info, xml_fname) in get_xml_map(repository_path, station_names).items():
        logger.info('applying transfer function {}'.format(station_name))
        logger.info('XML record = {}'.format(xml_fname))
        df_E = pd.DataFrame(index=df_B.index)
        B_north = df_B['B_' + N_suffix].values
        B_east = df_B['B_' + E_suffix].values
        E_north, E_east = apply_transfer_function(B_north,
                                                  B_east,
                                                  interval,
                                                  xml_fname,
                                                  extrapolate0=True)
        E_north_key = '{}_{}'.format(output_key_prefix, N_suffix)
        E_east_key = '{}_{}'.format(output_key_prefix, E_suffix)
        logger.info('writing columns {} and {} to {} with key {}'.format(E_north_key,
                                                                         E_east_key,
                                                                         output_hdf_fname,
                                                                         station_name))
        df_E[E_north_key] = E_north
        df_E[E_east_key] = E_east
        df_E.to_hdf(output_hdf_fname, station_name)
    return output_hdf_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Add 3-D modeled E field records to an HDF file containing B field data.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('hdf_fname',
                        type=str,
                        help='HDF file record to process')
    parser.add_argument('station_names',
                        type=str,
                        nargs='+',
                        metavar='station_name',
                        help='station identifiers for the 3-D EMTFs')
    parser.add_argument('--input-key',
                        '-i',
                        type=str,
                        default='B',
                        help='')
    parser.add_argument('--output-key-prefix',
                        '-o',
                        type=str,
                        default='E',
                        help='key prefix to associate with the processed records ({prefix}_{station_name})')
    parser.add_argument('--coordinate_frame',
                        '-c',
                        type=str,
                        choices=COORDINATE_FRAME_SUFFIX.keys(),
                        default='geo',
                        help=coordinate_frame_help(COORDINATE_FRAME_SUFFIX))
    parser.add_argument('--output_hdf_fname',
                        type=str,
                        default=None,
                        help='output HDF file name (if not specified, use the input HDF file)')
    parser.add_argument('--repository_path',
                        '-r',
                        type=str,
                        default='.',
                        help='path to search for EMTF repository')
    args = parser.parse_args(argv[1:])


    e3d2hdf(args.hdf_fname,
            args.station_names,
            input_key=args.input_key,
            output_key_prefix=args.output_key_prefix,
            #coordinate_frame=args.coordinate_frame,
            N_suffix=COORDINATE_FRAME_SUFFIX[args.coordinate_frame]['north'],
            E_suffix=COORDINATE_FRAME_SUFFIX[args.coordinate_frame]['east'],
            output_hdf_fname=args.output_hdf_fname,
            repository_path=args.repository_path)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
