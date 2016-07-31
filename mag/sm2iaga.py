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

from sm_stations import STATION_MAP


"""
NOTE: This code assumes that the input and output data are 1-minute interval.
"""


HEADER_TEMPLATE = """\
 Format                 IAGA-2002                                    |
 Source of Data         SuperMAG                                     |
 IAGA CODE              {stn}                                          |
 Geodetic Latitude      {lat:<8.3f}                                     |
 Geodetic Longitude     {lon:<8.3f}                                     |
 Elevation              {el:<8.3f}                                     |
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


def ne2xy_converter(stn, dt, nan=True, elevation=None):
    """ ??? """
    station_info = STATION_MAP[stn.upper()]
    lat = station_info.glat
    lon = station_info.glon
    point = Point(dt, lat, lon, elevation / 1e3 if elevation else 0)
    point.run_igrf()
    dec_deg = point.dec
    dec_rad = math.radians(dec_deg)
    cos_dec = math.cos(dec_rad)
    sin_dec = math.sin(dec_rad)
    def ne2xy(n, e):
        if n in [88888, 99999] or e in [88888, 99999]:
            if nan:
                return NP.nan, NP.nan
            else:
                return 88888, 88888
        x = cos_dec * n - sin_dec * e
        y = sin_dec * n + cos_dec * e
        return x, y
    return ne2xy



def fill_data(df,
              start_date,
              end_date,
              nan=True):
    """
    Fill data in missing rows of *df* in the time interval
    *start_date* to *end_date* (one minute sample period, closed on
    the left and open on the right)). If *nan* will with not a
    numbers, otherwise fill with 88888 (i.e., the missing value in
    IAGA2002 data records). Return the tuple of the list of date/times
    for each data record and the list of tuple values (containing the
    N, E, Z, and F in [nT] in that order).
    """
    data_map = {}
    for row in df.itertuples():
        N = row.N
        E = row.E
        Z = row.Z
        F = math.hypot(math.hypot(N, E), Z)
        data_map[row.Date_UTC] = (N, E, Z, F)
    # fill missing records
    filled_data = []
    dts = list(PD.date_range(start=start_date,
                             end=end_date,
                             freq='min'))[:-1]
    if nan:
        fill_value = (float('nan'),) * 4
    else:
        fill_value = (88888,) * 4
    for dt in dts:
        filled_data.append(data_map.get(dt,
                                        fill_value))
    return dts, filled_data


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
    dts, filled_data = fill_data(df,
                                 date,
                                 date + timedelta(days=1),
                                 nan=False)
    if not nez:
        ne2xy = ne2xy_converter(stn, dts[0], nan=False)

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
                C1, C2 = ne2xy(N, E)
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


def df_map2iaga(path, df_map, nez=False):
    """
    """
    output_map = defaultdict(dict)
    for stn, df in df_map.iteritems():
        station_info = STATION_MAP[stn.upper()]
        lat = station_info.glat
        lon = station_info.glon
        for date, df_date in df.groupby(df['Date_UTC'].map(lambda x: x.date())):
            output_fname = iaga_fname(path, stn, date)
            dataframe2iaga(output_fname,
                           stn,
                           date,
                           df_date,
                           lat,
                           lon,
                           elevation=0,
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
    parser.add_argument('--nez',
                        action='store_true',
                        help='store (raw) HEZ components (aligned to local magnetic field) instead of XYZ components')
    args = parser.parse_args(argv[1:])

    df_map = read_sm_csv(args.csv_fname)
    df_map2iaga(args.output_path,
                df_map,
                nez=args.nez)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
