from __future__ import division

import sys
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import numpy as NP
import pandas as PD
from geomagio.StreamConverter import get_obs_from_geo
from obspy.core.stream import Stream
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.trace import Trace

from fgm2iaga import parse
from iaga2hdf import get_dec_tenths_arcminute, write_hdf

logger = logging.getLogger('pyrsss.mat.fgm2hdf')


def build_header(data_list,
                 keys=[('geodetic_latitude', 'lat'),
                       ('geodetic_longitude', 'lon'),
                       ('station', 'siteid')],
                 elevation=None,
                 baseline_declination=None):
    """
    Build and return meta information mapping based on information in
    *data_list* with *keys*, *elevation*, and *baseline_declination*
    (the baseline declination, determined from IGRF if not specified).
    """
    header = {}
    for key, search_key in keys:
        values = set([getattr(x, search_key) for x in data_list])
        if len(values) > 1:
            raise ValueError('multiple values for {} encountered'.format(search_key))
        if len(values) == 0:
            logger.warning('{} not found'.format(key))
            continue
        value = values.pop()
        header[key] = value
    if 'station' in header:
        header['station'] = header['station'][:3]
    d1 = PD.to_datetime(data_list[0].index.values[0]).to_pydatetime()
    d2 = PD.to_datetime(data_list[-1].index.values[-1]).to_pydatetime()
    d1_obj = UTCDateTime('{:%Y-%m-%d %H:%H:%S}'.format(d1))
    d2_obj = UTCDateTime('{:%Y-%m-%d %H:%H:%S}'.format(d2))
    header['starttime'] = d1_obj
    header['endtime'] = d2_obj
    if elevation is None:
        logger.warning('no elevation found --- using default of 0')
        header['elevation'] = 0
    else:
        header['elevation'] = elevation
    delta = NP.diff(data_list[0].index.values[:2])[0] / NP.timedelta64(1, 's')
    fs = 1 / delta
    header['sampling_rate'] = fs
    if baseline_declination is None:
        d = {'starttime': header['starttime'],
             'Geodetic Latitude': header['geodetic_latitude'],
             'Geodetic Longitude': header['geodetic_longitude'],
             'Elevation': header['elevation']}
        baseline_declination = get_dec_tenths_arcminute(d, d1)
    header['declination_base'] = baseline_declination
    header['npts'] = sum(map(len, data_list))
    return header


def fgm2hdf(hdf_fname,
            fgm_fnames,
            he=False,
            elevation=0,
            key='B_raw'):
    """
    Convert data found in FGM files *fgm_fnames* to an HDF record at
    *hdf_fname*. Write to the HDF record associated with *key*. If
    *he*, store the h (mag north) and e (mag east) components. Use
    *elevation* in specifying the measurement location.
    """
    data_list = []
    for fgm_fname in fgm_fnames:
        logger.info('reading {}'.format(fgm_fname))
        data_list.append(parse(fgm_fname))
    header = build_header(data_list,
                          elevation=elevation)

    df = PD.concat([x[['x', 'y', 'z', 'f']] for x in data_list])
    df.rename(columns={'x': 'B_X',
                       'y': 'B_Y',
                       'z': 'B_Z',
                       'f': 'B_F'},
              inplace=True)
    write_hdf(hdf_fname, df, key, header)
    return hdf_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Convert FGM format data to HDF.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('hdf_fname',
                        type=str,
                        help='output HDF file name')
    parser.add_argument('fgm_fnames',
                        type=str,
                        metavar='fgm_fname',
                        nargs='*',
                        help='input FGM file (in time order)')

    args = parser.parse_args(argv[1:])

    fgm2hdf(args.hdf_fname,
            args.fgm_fnames)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
