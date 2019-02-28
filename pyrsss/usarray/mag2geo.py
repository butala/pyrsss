import sys
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import numpy as np
import pandas as pd

from pyglow import Point

from .info import get_info_map

logger = logging.getLogger('pyrsss.usarray.mag2geo')


def compute_mag2geo(df, dec_deg):
    """
    """
    cos_dec = np.cos(np.deg2rad(dec_deg))
    sin_dec = np.sin(np.deg2rad(dec_deg))
    # rotate magnetic field measurements
    B_N = df.B_N.values
    B_E = df.B_E.values
    B_X = cos_dec * B_N - sin_dec * B_E
    B_Y = sin_dec * B_N + cos_dec * B_E
    # rotate electric field measurements
    E_N = df.E_N.values
    E_E = df.E_E.values
    E_X = cos_dec * E_N - sin_dec * E_E
    E_Y = sin_dec * E_N + cos_dec * E_E
    return pd.DataFrame(index=df.index,
                        data={'B_X': B_X,
                              'B_Y': B_Y,
                              'E_X': E_X,
                              'E_Y': E_Y})


def get_declination(station_id):
    """
    Retrieve the magnetic declination angle at USArray site
    *station_id* using IGRF. Use the start date as the IGRF
    epoch. Return the angle in degrees.
    """
    info_map = get_info_map()
    station_info = info_map.loc[station_id]
    point = Point(station_info.start.to_pydatetime(),
                  station_info.lat,
                  station_info.lon,
                  station_info.elev / 1e3)
    point.run_igrf()
    dec = point.dec
    logger.info('IGRF declination angle at site {} on {:%Y-%m-%d} is {:.2f} [deg]'.format(station_id,
                                                                                          station_info.start,
                                                                                          dec))
    return dec


def mag2geo(hdf_fname,
            input_key,
            output_key,
            declination=None,
            station_id=None,
            output_fname=None):
    """
    """
    df = pd.read_hdf(hdf_fname, input_key)

    if declination is None:
        assert station_id is not None
        dec_deg = get_declination(station_id)
    else:
        dec_deg = declination

    if output_fname is None:
        output_fname = hdf_fname

    df_geo = compute_mag2geo(df, dec_deg)
    df_geo.to_hdf(output_fname, output_key)
    return output_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Rotate geomagnetic coordinate data to geographic',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('hdf_fname',
                        type=str,
                        help='HDF file')
    dec_group = parser.add_mutually_exclusive_group(required=True)
    dec_group.add_argument('--declination',
                           '-d',
                           type=float,
                           default=None,
                           help='declination angle (degrees)')
    dec_group.add_argument('--station_id',
                           '-s',
                           type=str,
                           default=None,
                           help='station identifier (compute declination at station location)')
    parser.add_argument('--input_key',
                        type=str,
                        default='iris',
                        help='key to search for B_N, B_E, E_N, and E_E')
    parser.add_argument('--output_key',
                        type=str,
                        default='geo',
                        help='key to write B_X, B_Y, E_X, and E_Y')
    parser.add_argument('--output_fname',
                        type=str,
                        default=None,
                        help='output HDF file (use the input file if not specified)')
    args = parser.parse_args(argv[1:])

    mag2geo(args.hdf_fname,
            args.input_key,
            args.output_key,
            declination=args.declination,
            station_id=args.station_id,
            output_fname=args.output_fname)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
