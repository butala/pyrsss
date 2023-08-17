import sys
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import datetime

import pandas as pd


def read_sm_csv(csv_fname):
    """
    Parse the SuperMAG CSV format data record *csv_fname*. For each
    station, store the information in pandas
    :class:`DataFrame`. Return a mapping between the station
    identifier and data frame.
    """
    df = pd.read_csv(csv_fname,
                     header=0,
                     parse_dates=[0],
                     date_parser=lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S'),
                     index_col=0)
    df_map = {name: group for name, group in df.groupby('IAGA')}
    for df in df_map.values():
        del df['IAGA']
        df.rename(columns={'N': 'B_N',
                           'E': 'B_E',
                           'Z': 'B_Z'},
                  inplace=True)
    return df_map


def sm2hdf(hdf_fname, csv_fname):
    """
    Parse the SuperMAG CSV format data record *csv_fname* and store
    the data in the HDF record *hdf_fname*. If data from only one
    station are found in *csv_fname* then the key "raw" is
    used. Otherwise, use key "{station ID}_raw".
    """
    df_map = read_sm_csv(csv_fname)
    if len(df_map) == 1:
        next(iter(df_map.values())).to_hdf(hdf_fname,
                                           key='raw',
                                           mode='w')
    else:
        for i, (stn, df_stn) in enumerate(df_map.items()):
            df_stn.to_hdf(hdf_fname,
                          key='{}_raw'.format(stn),
                          mode='w' if i == 0 else 'a')
    return hdf_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Convert SuperMAG CSV data record to HDF.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('hdf_fname',
                        type=str,
                        help='HDF file to store data')
    parser.add_argument('csv_fname',
                        type=str,
                        help='SuperMAG CSV data record')
    args = parser.parse_args(argv[1:])

    sm2hdf(args.hdf_fname,
           args.csv_fname)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
