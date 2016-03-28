import sys
import os
import logging
import math
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import datetime, timedelta
from collections import defaultdict, OrderedDict

import pandas as PD


"""
NOTE: This code assumes that the input and output data are 1-minute interval.
"""


HEADER_TEMPLATE = """\
 Format                 IAGA-2002                                    |
 Source of Data         SuperMAG                                     |
 IAGA CODE              {stn}                                          |
 Geodetic Latitude      XXX                                          |
 Geodetic Longitude     XXX                                          |
 Elevation              XXX                                          |
 Reported               {reported}                                         |
DATE       TIME         DOY     {stn}{C1}      {stn}{C2}      {stn}{C3}      {stn}F   |
"""


def read_sm_csv(csv_fname):
    """
    """
    df = PD.read_csv(csv_fname,
                     header=0,
                     parse_dates=[0],
                     date_parser=lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S'))
    return {name: group for name, group in df.groupby('IAGA')}


def get_reported(df):
    """
    """
    if all([x in df.columns for x in ['N', 'E', 'Z']]):
        return 'HDZF', {'C1': ('H', 'N'),
                        'C2': ('D', 'E'),
                        'C3': ('Z', 'Z')}
    raise ValueError('could not determine coordinate system from columns: ' + ', '.join(df.columns))


def dataframe2iaga(fname, stn, date, df, header_template=HEADER_TEMPLATE):
    """
    ASSUMES DF SPANS ONLY A SINGLE DAY
    """
    reported, coord_map = get_reported(df)
    # build map of data records
    data_map = {}
    for _, row in df.iterrows():
        C1 = row[coord_map['C1'][1]]
        C2 = row[coord_map['C2'][1]]
        C3 = row[coord_map['C3'][1]]
        F = math.hypot(math.hypot(C1, C2), C3)
        data_map[row.Date_UTC] = (C1, C2, C3, F)
    # fill missing records
    filled_data = []
    dts = list(PD.date_range(start=date,
                             end=date + timedelta(days=1),
                             freq='min'))[:-1]
    for dt in dts:
        filled_data.append(data_map.get(dt,
                                        (88888,) * 4))
    # write data records to file
    with open(fname, 'w') as fid:
        fid.write(header_template.format(stn=stn,
                                         reported=reported,
                                         C1=coord_map['C1'][0],
                                         C2=coord_map['C2'][0],
                                         C3=coord_map['C3'][0]))
        for dt, (C1, C2, C3, F) in zip(dts, filled_data):
            fid.write('{date:%Y-%m-%d %H:%M:%S.000} {date:%j}'
                      '    {C1:>9.2f} {C2:>9.2f} {C3:>9.2f} {F:>9.2f}\n'.format(date=dt,
                                                                                C1=C1,
                                                                                C2=C2,
                                                                                C3=C3,
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


def df_map2iaga(path, df_map):
    """
    """
    output_map = defaultdict(dict)
    for stn, df in df_map.iteritems():
        for date, df_date in df.groupby(df['Date_UTC'].map(lambda x: x.date())):
            output_fname = iaga_fname(path, stn, date)
            dataframe2iaga(output_fname, stn, date, df_date)
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
    args = parser.parse_args(argv[1:])

    df_map = read_sm_csv(args.csv_fname)
    df_map2iaga(args.output_path,
                df_map)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
