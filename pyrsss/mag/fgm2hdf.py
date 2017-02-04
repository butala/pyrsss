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
from iaga2hdf import get_dec_tenths_arcminute, stream2df

logger = logging.getLogger('pyrsss.mat.fgm2hdf')


def build_header(data_list,
                 keys=[('geodetic_latitude', 'lat'),
                       ('geodetic_longitude', 'lon'),
                       ('station', 'siteid')],
                 elevation=None,
                 baseline_declination=None):
    """
    ???
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
             'Geodetic_Latitude': header['geodetic_latitude'],
             'Geodetic_Longitude': header['geodetic_longitude'],
             'Elevation': header['elevation']}
        baseline_declination = get_dec_tenths_arcminute(d)
    header['declination_base'] = baseline_declination
    header['npts'] = sum(map(len, data_list))
    return header


def build_stream(data_list,
                 network='NT',
                 location='R0',
                 he=False,
                 elevation=None):
    """
    """
    header = build_header(data_list,
                          elevation=elevation)
    traces = []
    for channel in ['x', 'y', 'z', 'f']:
        vals = []
        for data_i in data_list:
            vals.extend(data_i[channel])
        header_i = header.copy()
        header_i['channel'] = channel.upper()
        header_i['network'] = network
        header_i['location'] = location
        traces.append(Trace(data = NP.array(vals),
                            header = header_i))
    return Stream(traces=traces)


def fgm2hdf(hdf_fname,
            fgm_fnames,
            he=False):
    """
    ???
    """
    data_list = []
    for fgm_fname in fgm_fnames:
        logger.info('reading {}'.format(fgm_fname))
        data_list.append(parse(fgm_fname))
    geo = build_stream(data_list)
    obs = get_obs_from_geo(geo) if he else False
    df = stream2df(geo, he=obs)
    df.to_hdf(hdf_fname, 'B', mode='w')
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
    parser.add_argument('--he',
                        action='store_true',
                        help='include data in HE coordinate')

    args = parser.parse_args(argv[1:])

    fgm2hdf(args.hdf_fname,
            args.fgm_fnames,
            he=args.he)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
