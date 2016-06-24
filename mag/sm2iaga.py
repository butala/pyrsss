from __future__ import division

import sys
import os
import logging
import math
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import datetime, timedelta
from collections import defaultdict, OrderedDict

import pandas as PD

from pyglow.pyglow import Point


"""
NOTE: This code assumes that the input and output data are 1-minute interval.
"""


HEADER_TEMPLATE = """\
 Format                 IAGA-2002                                    |
 Source of Data         SuperMAG                                     |
 IAGA CODE              {stn}                                          |
 Geodetic Latitude      {lat:<8.3f}                                     |
 Geodetic Longitude     {lon:<8.3f}                                     |
 Elevation              {el}                                          |
 Reported               {reported}                                         |
 Digital Sampling       1-Minute                                     |
DATE       TIME         DOY     {stn}{C1}       {stn}{C2}     {stn}Z      {stn}F   |
"""


def read_sm_csv(csv_fname):
    """
    """
    df = PD.read_csv(csv_fname,
                     header=0,
                     parse_dates=[0],
                     date_parser=lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S'))
    return {name: group for name, group in df.groupby('IAGA')}


def dataframe2iaga(fname,
                   stn,
                   date,
                   df,
                   lat,
                   lon,
                   elevation=None,
                   nez=False,
                   header_template=HEADER_TEMPLATE):
    """
    ASSUMES DF SPANS ONLY A SINGLE DAY
    """
    #reported, coord_map = get_reported(df)
    # build map of data records
    data_map = {}
    for _, row in df.iterrows():
        N = row['N']
        E = row['E']
        Z = row['Z']
        F = math.hypot(math.hypot(N, E), Z)
        data_map[row.Date_UTC] = (N, E, Z, F)
    # fill missing records
    filled_data = []
    dts = list(PD.date_range(start=date,
                             end=date + timedelta(days=1),
                             freq='min'))[:-1]
    for dt in dts:
        filled_data.append(data_map.get(dt,
                                        (88888,) * 4))
    if not nez:
        # get magnetic declination angle for local magnetic (NEZ) to
        # geographic (XYZ) conversion
        point = Point(dts[0], lat, lon, elevation / 1e3 if elevation else 0)
        point.run_igrf()
        dec_deg = point.dec
        dec_rad = math.radians(dec_deg)
        cos_dec = math.cos(dec_rad)
        sin_dec = math.sin(dec_rad)

    # write data records to file
    with open(fname, 'w') as fid:
        fid.write(header_template.format(stn=stn,
                                         reported='NEZF' if nez else 'XYZF',
                                         lat=lat,
                                         lon=lon,
                                         el='XXX' if elevation is None else elevation,
                                         C1='N' if nez else 'X',
                                         C2='E' if nez else 'Y'))
        for dt, (N, E, Z, F) in zip(dts, filled_data):
            if nez:
                C1 = N
                C2 = E
            else:
                C1 = cos_dec * N - sin_dec * E
                C2 = sin_dec * N + cos_dec * E
            fid.write('{date:%Y-%m-%d %H:%M:%S.000} {date:%j}'
                      '    {C1:>9.2f} {C2:>9.2f} {Z:>9.2f} {F:>9.2f}\n'.format(date=dt,
                                                                               C1=C1,
                                                                               C2=C2,
                                                                               Z=Z,
                                                                               F=F))
    return fname


def iaga_fname(path, stn, date, cutoff_date=datetime(2007, 1, 1)):
    """
    """
    if date >= cutoff_date.date():
        ext = 'dmin.min'
    else:
        ext = 'd.min'
    return os.path.join(path, '{stn}{date:%Y%m%d}{ext}'.format(stn=stn.lower(),
                                                               date=date,
                                                               ext=ext))


def df_map2iaga(path, df_map, lat, lon, elevation=None, nez=False):
    """
    """
    output_map = defaultdict(dict)
    for stn, df in df_map.iteritems():
        for date, df_date in df.groupby(df['Date_UTC'].map(lambda x: x.date())):
            output_fname = iaga_fname(path, stn, date)
            dataframe2iaga(output_fname,
                           stn,
                           date,
                           df_date,
                           lat,
                           lon,
                           elevation=elevation,
                           nez=nez)
            output_map[stn][date] = output_fname
    return output_map


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Convert SuperMAG CSV format data to daily IAGA2002 format.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_path',
                        type=str,
                        help='path to store daily IAGA2002 format files')
    parser.add_argument('csv_fname',
                        type=str,
                        help='SuperMAG CSV file')
    parser.add_argument('latitude',
                        type=float,
                        help='site geodetic latitude [deg]')
    parser.add_argument('longitude',
                        type=float,
                        help='site longitude [deg]')
    parser.add_argument('--elevation',
                        '-e',
                        type=float,
                        default=None,
                        help='site elevation [m]')
    parser.add_argument('--nez',
                        action='store_true',
                        help='store (raw) HEZ components (aligned to local magnetic field) instead of XYZ components')
    args = parser.parse_args(argv[1:])

    df_map = read_sm_csv(args.csv_fname)
    df_map2iaga(args.output_path,
                df_map,
                args.latitude,
                args.longitude,
                elevation=args.elevation,
                nez=args.nez)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
