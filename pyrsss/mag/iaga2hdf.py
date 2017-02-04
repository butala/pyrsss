import sys
import math
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import numpy as NP
import pandas as PD
from pyglow.pyglow import Point
from geomagio.StreamConverter import get_obs_from_geo
from obspy.core.stream import Stream
from obspy.core.utcdatetime import UTCDateTime
from obspy.core.trace import Trace

from iaga2002 import parse
from ..util.angle import deg2tenths_of_arcminute

logger = logging.getLogger('pyrsss.mag.iaga2hdf')


def find_value(key, headers):
    """
    Search *headers* for *key* and return the common value shared
    across each element of *headers* and raise an exception if
    multiple values are encountered.
    """
    values = set([getattr(x, key, None) for x in headers if getattr(x, key, None) is not None])
    if len(values) == 0:
        raise KeyError('{} not found'.format(key))
    elif len(values) > 1:
        raise ValueError('multiple values of {} detected ({})'.format(key,
                                                                      ', '.join(map(str, values))))
    return values.pop()


def reduce_headers(headers,
                   keys=['IAGA_CODE',
                         'Geodetic_Latitude',
                         'Geodetic_Longitude',
                         'Elevation',
                         'Reported',
                         'Comment']):
    """
    Search *headers* for values associated with each of *keys* and
    confirm those values are equal. Return a mapping between *keys*
    and each unique value.
    """
    return {k: find_value(k, headers) for k in keys}


def fix_sign(x, N=360 * 60 * 10):
    """
    Convert negative tenths of arcminutes *x* to positive by checking
    bounds and taking the modulus N (360 degrees * 60 minutes per
    degree * 10 tenths per 1).
    """
    if x < 0:
        assert x > -N
        x += N
    assert x < N
    return x % N


def build_stream(header,
                 data_maps,
                 dec_tenths_arcminute,
                 network='NT',
                 location='R0',
                 default_elevation=0,
                 he=False):
    """
    Build obspy :class:`Stream` from *header* information and the data
    arranged in *data_maps*. Use *dec_tenths_arcminute* (local
    magnetic declination in determining magnetic north and east or XYZ
    in unite of tenths of arcminutes).
    """
    glon = header['Geodetic_Longitude']
    if glon < 0:
        glon += 360
    delta = NP.diff(data_maps[0].keys()[:2])[0].total_seconds()
    fs = 1 / delta
    d1 = data_maps[0].keys()[0]
    d2 = data_maps[-1].keys()[-1]
    d1_obj = UTCDateTime('{:%Y-%m-%d %H:%H:%S}'.format(d1))
    d2_obj = UTCDateTime('{:%Y-%m-%d %H:%H:%S}'.format(d2))
    N = sum(map(len, data_maps))
    stream_header = {'geodetic_latitude': header['Geodetic_Latitude'],
                     'geodetic_longitude': glon,
                     'station': header['IAGA_CODE'],
                     'sampling_rate': fs,
                     'starttime': d1_obj,
                     'endtime': d2_obj,
                     'declination_base': dec_tenths_arcminute,
                     'npts': N}
    try:
        stream_header['elevation'] = header['elevation']
    except KeyError:
        stream_header['elevation'] = default_elevation
        logger.warning('elevation is unknown --- inserting {}'.format(default_elevation))
    traces = []
    for channel in ['x', 'y', 'z', 'f']:
        vals = [getattr(x, channel) for dm in data_maps for x in dm.itervalues()]
        header_i = stream_header.copy()
        header_i['channel'] = channel.upper()
        header_i['network'] = network
        header_i['location'] = location
        traces.append(Trace(data = NP.array(vals),
                            header = header_i))
    return Stream(traces=traces)


def stream2dt(stream):
    """
    Verify time indices of all *stream* traces are equal and return a
    list of sample of time suitable for a :class:`DataFrame` index.
    """
    dt_match = None
    for t in stream.traces:
        if dt_match is None:
            dt_match = [(t.stats.starttime + x).datetime for x in t.times()]
            continue
        dt_i = [(t.stats.starttime + x).datetime for x in t.times()]
        if dt_match != dt_i:
            raise ValueError('timing discrepancy detected in stream elements')
    return dt_match


def stream2df(geo, he=False):
    """
    Build and return :class:`DataFrame` from *geo* stream elements. By
    default, these are traces x, y, z, and f. If *he*, also include
    the mag north, h, and mag east, e, components.
    """
    if len(set(map(len, geo.traces))) != 1:
        raise ValueError('size mismatch detected (dt and all geo traces should be the same length)')
    data = {'B' + k: geo.select(channel=k).traces[0].data for k in ['x', 'y', 'z', 'f']}
    if he:
        for k in ['h', 'e']:
            data['B' + k] = he.select(channel=k).traces[0].data
    return PD.DataFrame(index=stream2dt(geo),
                        data=data)



def find_decbas(header):
    """
    Search IAGA *header* Comment for DECBAS. Return either this
    quantity (declination in tenths of arcmin) or `None` if not found.
    """
    if 'Comment' in header:
        if 'DECBAS' in header['Comment']:
            for line in header['Comment'].split('\n'):
                if line.startswith('DECBAS'):
                    _, decbas_str, _ = line.split(None, 2)
                    return float(decbas_str)
    return None


def iaga2hdf(hdf_fname,
             iaga2002_fnames,
             he=False,
             strict=False):
    """
    Convert data found in IAGA 2002 files *iaga2002_fname* to an HDF
    record. If *he*, store the h (mag north) and e (mag east)
    components. If *strict*, use strict conformity checks when parsing
    the IAGA 2002 records.
    """
    headers = []
    data_maps = []
    for fname in iaga2002_fnames:
        logger.info('reading {}'.format(fname))
        header_i, data_map_i = parse(fname, strict=strict)
        headers.append(header_i)
        data_maps.append(data_map_i)
    header = reduce_headers(headers)
    if header['Reported'] != 'XYZF':
        raise NotImplementedError('only XYZF implemented so far')
    if he and header['Reported'] == 'XYZF':
        decbas = find_decbas(header)
        if decbas:
            dec_tenths_arcminute = decbas
        else:
            # transform to obs coordinates
            date = data_maps[0].keys()[0]
            point = Point(date,
                          header['Geodetic_Latitude'],
                          header['Geodetic_Longitude'],
                          header['Elevation'])
            point.run_igrf()
            dec_deg = point.dec
            dec_tenths_arcminute = fix_sign(deg2tenths_of_arcminute(dec_deg))
    geo = build_stream(header, data_maps, dec_tenths_arcminute, he=he)
    obs = get_obs_from_geo(geo) if he else False
    df = stream2df(geo, he=obs)
    df.to_hdf(hdf_fname, 'B', mode='w')
    return hdf_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Convert IAGA2002 magnetometer data to HDF.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('hdf_fname',
                        type=str,
                        help='output HDF file name')
    parser.add_argument('iaga2002_fnames',
                        type=str,
                        metavar='iaga2002_fname',
                        nargs='*',
                        help='input IAGA2002 file (in time order)')
    parser.add_argument('--he',
                        action='store_true',
                        help='include data in HE coordinate')
    args = parser.parse_args(argv[1:])

    iaga2hdf(args.hdf_fname,
             args.iaga2002_fnames,
             he=args.he)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
