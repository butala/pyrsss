import logging
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict

import pandas as pd

from ..stats.stats import despike, max_consectuive_true
from ..mag.process_hdf import fill_nans, second_interval_filter, minute_interval_filter
from ..signal.lfilter import lp_fir_filter

logger = logging.getLogger('pyrsss.usarray.preprocess')


def nan_interpolate(df):
    """
    """
    I = df.isna()
    cumsum, consecutive = max_consectuive_true(I)
    for column, cumsum_i, consecutive_i in zip(df.columns, cumsum, consecutive):
        if cumsum_i > 0:
            logger.info('column {}: {} points ({:.2f}%) removed --- max run removed = {} points'.format(column, cumsum_i, cumsum_i / len(df) * 100, consecutive_i))
    if any(cumsum):
        df = df.interpolate(inplace=True)
    return df


"""
"""
BANDPASS_MAP = {'1to100mHz': {1: second_interval_filter,
                              60: minute_interval_filter}}


def preprocess(hdf_fname,
               input_key='iris',
               output_key='filtered',
               columns='all',
               output_hdf_fname=None,
               apply_despike=True,
               apply_median_subtraction=True,
               bandpass_filter='1to100mHz'):
    """
    """
    logger.info('preprocessing {} key={}'.format(hdf_fname, input_key))
    # read data
    df = pd.read_hdf(hdf_fname, input_key)
    # select columns
    if columns == 'all':
        columns = [x for x in df.columns]
    logger.info('processing columns {}'.format(columns))
    df_col = df[columns]
    # outlier removal
    if apply_despike:
        logger.info('applying data despike')
        df_col = despike(df_col)
    # fill gaps
    logger.info('filling gaps with NaNs')
    df_col, delta = fill_nans(df_col)
    # interpolate across gaps
    logger.info('gap interpolation')
    df_col = nan_interpolate(df_col)
    # apply bandpass filter
    logger.info('apply bandpass filter')
    h = BANDPASS_MAP[bandpass_filter][delta.seconds]()
    data = OrderedDict()
    dt = df_col.index
    for i, column in enumerate(columns):
        if i == 0:
            (col_filtered,
             dt_filtered) = lp_fir_filter(h,
                                          df_col[column].values,
                                          mode='valid',
                                          index=dt)
        else:
            col_filtered = lp_fir_filter(h,
                                         df_col[column].values,
                                         mode='valid')
        data[column] = col_filtered
    df_filtered = pd.DataFrame(index=dt_filtered,
                               data=data)
    # median subtraction
    if apply_median_subtraction:
        logger.info('median subtraction')
        df_filtered = df_filtered - df_filtered.median()
    # output result
    if output_hdf_fname is None:
        output_hdf_fname = hdf_fname
    logger.info('storing preprocessed data to {} with key={}'.format(output_hdf_fname, output_key))
    df_filtered.to_hdf(output_hdf_fname, output_key)
    return output_hdf_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Preprocess data frame columns (time series).',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('hdf_fname',
                        type=str,
                        help='HDF file record to process')
    parser.add_argument('--input-key',
                        '-i',
                        type=str,
                        default='iris',
                        help='HDF record key for the data to preprocess')
    parser.add_argument('--output-key',
                        '-o',
                        type=str,
                        default='filtered',
                        help='HDF record key to store the preprocessed data')
    parser.add_argument('--columns',
                        '-c',
                        type=str,
                        nargs='*',
                        default='all',
                        help='columns to preprocess (all columns by default)')
    parser.add_argument('--output-hdf-fname',
                        type=str,
                        default=None,
                        help='output HDF file name (if not specified, use the input HDF file)')
    parser.add_argument('--no-despike',
                        action='store_false',
                        help='disable the despike step')
    parser.add_argument('--no-median-subtraction',
                        action='store_false',
                        help='disable the median subtraction step')
    parser.add_argument('--bandpass_filter',
                        '-b',
                        choices=BANDPASS_MAP.keys(),
                        default='1to100mHz',
                        help='choice of bandpass filter')
    args = parser.parse_args(argv[1:])

    preprocess(args.hdf_fname,
               input_key=args.input_key,
               output_key=args.output_key,
               columns=args.columns,
               output_hdf_fname=args.output_hdf_fname,
               apply_despike=args.no_despike,
               apply_median_subtraction=args.no_median_subtraction,
               bandpass_filter=args.bandpass_filter)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
