import sys
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import numpy as np
import pandas as pd
from .hdfio import (get_dec_tenths_arcminute, hdx_to_xyz, read_hdf,
                    write_hdf, xyz_to_hdx)
from .iaga2002 import iaga2df

logger = logging.getLogger('pyrsss.mag.iaga2hdf')


def find_value(key, headers):
    """
    Search *headers* for *key* and return the common value shared
    across each element of *headers* and raise an exception if
    multiple values are encountered.
    """
    values = set([x[key] for x in headers])
    if len(values) == 0:
        raise KeyError('{} not found'.format(key))
    elif len(values) > 1:
        raise ValueError('multiple values of {} detected ({})'.format(key,
                                                                      ', '.join(map(str, values))))
    return values.pop()


def reduce_headers(headers,
                   keys=['IAGA CODE',
                         'Geodetic Latitude',
                         'Geodetic Longitude',
                         'Elevation',
                         'Reported',
                         'decbas']):
    """
    Search *headers* for values associated with each of *keys* and
    confirm those values are equal. Return a mapping between *keys*
    and each unique value.
    """
    reduced_header = {}
    for key in keys:
        try:
            value = find_value(key, headers)
        except KeyError:
            logger.warning('Entry for {} not found in IAGA headers --- skipping'.format(key))
            continue
        reduced_header[key] = value
    return reduced_header


def df2stream(df,
              header,
              network='NT',
              location='R0',
              radians=True,
              default_elevation=0):
    # obspy is an `usarray`-extra dependency and only this function needs
    # it; importing here keeps the module loadable without it.
    from obspy.core.stream import Stream
    from obspy.core.trace import Trace
    from obspy.core.utcdatetime import UTCDateTime
    """
    Build and return obspy :class:`Stream` from *header* information
    and the :class:`DataFrame` *df*. Use *dec_tenths_arcminute* (local
    magnetic declination in determining magnetic north and east or XYZ
    in units of tenths of arcminutes). If *radians*, angles in D are
    given in degrees and must be converted to radians. If the site
    elevation is not included in the header, use
    *default_elevation*. The *network* and *location* identifiers are
    used in forming the trace names in the resultant stream.
    """
    glon = header['Geodetic Longitude']
    if glon < 0:
        glon += 360
    delta = (df.index[1] - df.index[0]).total_seconds()
    fs = 1 / delta
    d1 = df.index[0]
    d2 = df.index[-1]
    d1_obj = UTCDateTime('{:%Y-%m-%d %H:%H:%S}'.format(d1))
    d2_obj = UTCDateTime('{:%Y-%m-%d %H:%H:%S}'.format(d2))
    dec_tenths_arcminute = header.get('decbas',
                                      get_dec_tenths_arcminute(header,
                                                               d1.to_pydatetime()))
    logger.info('using declination baseline = {:.1f} (tenths of arcminutes)'.format(dec_tenths_arcminute))
    N = df.shape[0]
    stream_header = {'geodetic_latitude': header['Geodetic Latitude'],
                     'geodetic_longitude': glon,
                     'station': header['IAGA CODE'],
                     'sampling_rate': fs,
                     'starttime': d1_obj,
                     'endtime': d2_obj,
                     'declination_base': dec_tenths_arcminute,
                     'npts': N}
    try:
        stream_header['elevation'] = header['Elevation']
    except KeyError:
        stream_header['elevation'] = default_elevation
        logger.warning('elevation is unknown --- inserting {}'.format(default_elevation))
    traces = []
    for column in df.columns:
        channel = column[-1]
        header_i = stream_header.copy()
        header_i['channel'] = channel.upper()
        header_i['network'] = network
        header_i['location'] = location
        vals = df[column].values
        if channel == 'D' and radians:
            vals = np.radians(vals)
        traces.append(Trace(data=vals,
                            header=header_i))
    return Stream(traces=traces)


def combine_iaga(iaga2002_fnames):
    """
    Load one or more IAGA-2002 data records *iaga_fnames* and
    concatenate the data into a single :class:`DataFrame` record.
    """
    df_list = []
    header_list = []
    for iaga2002_fname in iaga2002_fnames:
        df_i, header_i = iaga2df(iaga2002_fname)
        df_list.append(df_i)
        header_list.append(header_i)
    return pd.concat(df_list), reduce_headers(header_list)


def _dec_degrees(df, header):
    """
    The station declination driving the component rotation, in degrees.

    The IAGA header's ``decbas`` baseline (tenths of arcminutes) wins when
    present; otherwise IGRF at the header's site and the record's first
    epoch (see ``hdfio.get_dec_tenths_arcminute``).
    """
    return get_dec_tenths_arcminute(
        header, pd.Timestamp(df.index[0]).to_pydatetime()) / (60.0 * 10.0)


def xy2df(df, header):
    """
    Add `B_X` and `B_Y` (surface magnetic field in geographic
    coordinates, X is north and Y is east) to the :class:`DataFrame`
    *df* by the declination rotation, and return. The record *header*
    supplies the declination (its ``decbas`` or IGRF).

    Needs ``B_H`` and ``B_E``: this is the HDZ -> XYZ rotation of the
    horizontal components (geomagio's ``get_geo_from_obs``, replaced by
    ``hdfio.hdx_to_xyz``).
    """
    for col in ('B_H', 'B_E', 'B_Z'):
        if col not in df.columns:
            raise ValueError(
                f'xy2df needs {col}; the HDZ->XYZ rotation takes B_H, B_E '
                f'and B_Z (columns found: {list(df.columns)})')
    dec_deg = _dec_degrees(df, header)
    X, Y, _ = hdx_to_xyz(df['B_H'].values, df['B_E'].values,
                         df['B_Z'].values, dec_deg)
    return df.assign(B_X=X, B_Y=Y)


def he2df(df, header):
    """
    Add `B_H` and `B_E` (surface magnetic field in local geomagnetic
    coordinates, H is local north and E is local east) to the
    :class:`DataFrame` *df* by the declination rotation, and return.
    The record *header* supplies the declination.

    Needs ``B_X`` and ``B_Y``: this is the XYZ -> HDZ rotation
    (geomagio's ``get_obs_from_geo``, replaced by ``hdfio.xyz_to_hdx``).
    """
    for col in ('B_X', 'B_Y', 'B_Z'):
        if col not in df.columns:
            raise ValueError(
                f'he2df needs {col}; the XYZ->HDZ rotation takes B_X, B_Y '
                f'and B_Z (columns found: {list(df.columns)})')
    if ('B_F' not in df.columns):
        values = zip(df['B_X'].values,
                     df['B_Y'].values,
                     df['B_Z'].values)
        df = df.assign(B_F=[np.linalg.norm([x, y, z]) for x, y, z in values])
    dec_deg = _dec_degrees(df, header)
    H, E, _ = xyz_to_hdx(df['B_X'].values, df['B_Y'].values,
                         df['B_Z'].values, dec_deg)
    return df.assign(B_H=H, B_E=E)


def add_columns(df, header, xy, he):
    """
    Add columns to *df* as needed and return a new :class:`DataFrame`,
    reporting magnetic field in additional coordinates systems. If
    *xy*, ensure data are provided in geographic XYZ coordinates. If
    *he*, ensure data are provided in local geomagnetic HEZ
    coordinates.
    """
    assert 'B_Z' in df.columns
    if xy:
        if not all([x in df.columns for x in ['B_X', 'B_Y']]):
            df = xy2df(df, header)
    if he:
        if not all([x in df.columns for x in ['B_H', 'B_E']]):
            df = he2df(df, header)
    return df


def iaga2hdf(hdf_fname,
             iaga2002_fnames,
             xy=True,
             he=False,
             key='B_raw'):
    """
    Convert data found in IAGA 2002 files *iaga2002_fnames* to an HDF
    record at *hdf_fname*. Write to the HDF record associated with
    *key*. If *he*, store the H (mag north) and E (mag east)
    components.
    """
    df, header = combine_iaga(iaga2002_fnames)
    df = add_columns(df, header, xy, he)
    # change header names
    if 'Geodetic Latitude' in header:
        header['geodetic_latitude'] = header.pop('Geodetic Latitude')
    if 'Geodetic Longitude' in header:
        header['geodetic_longitude'] = header.pop('Geodetic Longitude')
    write_hdf(hdf_fname, df, key, {k.lower(): v for k, v in header.items()})
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
