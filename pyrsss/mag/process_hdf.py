import sys
import logging
import warnings
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import timedelta
from itertools import groupby
from collections import OrderedDict

import numpy as np
import pandas as pd
import scipy.signal

from .iaga2hdf import read_hdf, write_hdf
from .filters import minute_interval_filter, second_interval_filter
from ..util.nan import nan_interp
from ..signal.lfilter import lp_fir_filter
from ..stats.stats import despike

logger = logging.getLogger('pyrsss.mag.process_hdf')


"""
A NOTE ON THE COORDINATE FRAME

The following describes the XYZ convention
(from https://github.com/usgs/geomag-algorithms/blob/master/docs/algorithms/XYZ.md):

- X is the magnitude of the geographic north pole component of the H vector;
- Y is the magnitude of the east component of the H vector;
- Z is the downward component of the geomagnetic field, same as before.
"""

def consecutive_nans(x, y):
    """
    Return the maximum number of nans encountered in the logical or of
    the *x* and *y* time series.
    """
    nan = np.isnan(x) | np.isnan(y)
    lengths = []
    for key, group in groupby(nan):
        if key == True:
            lengths.append(len(list(group)))
    if lengths:
        return max(lengths)
    return 0


def process_timeseries(dt,
                       Bx,
                       By,
                       c1='B_X',
                       c2='B_Y',
                       despike_data=True,
                       remove_mean=True):
    """
    Process surface magnetic field measurement time series with
    indices *dt* and components *Bx* and *By*. Output a
    :class:`DataFrame` with columns *c1* and *c2* associated with the
    processed output time series. If *despike_data*, remove outliers
    prior to filtering. If *remove_mean*, remove the mean from the
    output time series.
    """
    warnings.warn('use process_df instead',
                  PendingDeprecationWarning)
    n = consecutive_nans(Bx, By)
    interval = (dt[1] - dt[0]).total_seconds()
    if n > 0:
        logger.warning('longest contiguous gap = {:.2f} minutes'.format(n * interval / 60))
    # fill data gaps via linear interpolation
    Bx = nan_interp(Bx)
    By = nan_interp(By)
    if despike_data:
        # remove outliers
        df = pd.DataFrame(index=dt,
                          data={'Bx': Bx,
                                'By': By})
        df = despike(df)
        dt = df.index.to_pydatetime()
        Bx = df.Bx.values
        By = df.By.values
    # apply 1 - 100 mHz bandpass filter
    if interval == 1.0:
        h = second_interval_filter()
    elif interval == 60.0:
        h = minute_interval_filter()
    else:
        raise ValueError('1 to 100 mHz filter not yet synthesized for {} s interval data'.format(interval))
    Bx_filtered, dt_filtered = lp_fir_filter(h, Bx, mode='valid', index=dt)
    By_filtered = lp_fir_filter(h, By, mode='valid')
    # remove mean
    if remove_mean:
        Bx_filtered -= np.mean(Bx_filtered)
        By_filtered -= np.mean(By_filtered)
    # build DataFrame and store to disk
    return pd.DataFrame(index=dt_filtered,
                        data={c1: Bx_filtered,
                              c2: By_filtered})


def process(hdf_fname,
            source_key='B_raw',
            key='B',
            he=False,
            despike_data=True,
            remove_mean=True):
    """
    Process the magnetic field columns of *hdf_fname*, applying
    pre-processing (nan interpolation) and a band-pass filter. Look
    for input at *source_key* and store output under identifier
    *key*. If *he*, process the H and E magnetic field components. If
    *remove_mean*, remove the mean from each column.
    """
    logger.info('processing {}'.format(hdf_fname))
    df_raw, header = read_hdf(hdf_fname, source_key)
    dt = df_raw.index.to_pydatetime()
    Bx_raw = df_raw['B_X'].values * 1e-9
    By_raw = df_raw['B_Y'].values * 1e-9
    df_filtered = process_timeseries(dt,
                                     Bx_raw,
                                     By_raw,
                                     despike_data=despike_data,
                                     remove_mean=remove_mean)
    if he:
        Bh_raw = df_raw['B_H'].values * 1e-9
        Be_raw = df_raw['B_E'].values * 1e-9
        df_he_filtered = process_timeseries(dt,
                                            Bh_raw,
                                            Be_raw,
                                            c1='B_H',
                                            c2='B_E',
                                            despike_data=despike_data,
                                            remove_mean=remove_mean)
        df_filtered = df_filtered.join(df_he_filtered)
    write_hdf(hdf_fname, df_filtered, key, header)
    return hdf_fname


def fill_nans(df, delta=None):
    """
    """
    if not delta:
        dt_diff = np.diff(df.index.values)
        delta_timedelta64 = min(dt_diff)
        delta_seconds = delta_timedelta64 / np.timedelta64(1, 's')
        delta = timedelta(seconds=delta_seconds)
    logger.info('using delta = {} (s)'.format(delta.total_seconds()))
    index_new = pd.date_range(start=df.index[0],
                              end=df.index[-1],
                              freq=delta)
    missing = sorted(set(index_new) - set(df.index))
    if missing:
        logger.warning('Missing time indices (filled by NaNs):')
        for x in missing:
            logger.warning(x)
        df = df.reindex(index_new, copy=False)
    return df, delta


def nan_interpolate(df):
    """

    Reference:
    https://stackoverflow.com/questions/29007830/identifying-consecutive-nans-with-pandas
    """
    sum_nan = df.isnull().sum()
    df_null_int = df.isnull().astype(int)
    for col in df.columns:
        max_run = df[col].isnull().astype(int).groupby(df[col].notnull().astype(int).cumsum()).sum()
        if sum_nan[col]:
            # BELOW IS BROKEN!!!
            pass
            # logger.warning('column {} has {} NaNs ({} max consecutive run)'.format(col,
                                                                                   # sum_nan[col],
                                                                                   # max_run))
    df.interpolate(inplace=True)
    return df


def process_df(df,
               delta=None,
               despike_data=True,
               subtract_median=True):
    """
    """
    if despike_data:
        logger.info('despike')
        df = despike(df)
    logger.info('Fill gaps')
    df, delta = fill_nans(df, delta=delta)
    logger.info('Gap interpolation')
    df = nan_interpolate(df)
    # apply 1 - 100 mHz bandpass filter
    interval = delta.total_seconds()
    if interval == 1.0:
        h = second_interval_filter()
    elif interval == 60.0:
        h = minute_interval_filter()
    else:
        raise ValueError('1 to 100 mHz filter not yet synthesized for {} s interval data'.format(interval))
    data = OrderedDict()
    dt = df.index.to_pydatetime()
    for i, col in enumerate(df.columns):
        logger.info('Band-pass filter {}'.format(col))
        if i == 0:
            col_filtered, dt_filtered = lp_fir_filter(h, df[col].values, mode='valid', index=dt)
        else:
            col_filtered = lp_fir_filter(h, df[col].values, mode='valid')
        data[col] = col_filtered
    df_filtered = pd.DataFrame(index=dt_filtered,
                               data=data)
    # remove median
    if subtract_median:
        logger.info('Subtract median')
        df_filtered = df_filtered.sub(df.median(axis=1), axis=0)
    return df_filtered


def process_new(hdf_fname,
                source_key='B_raw',
                key='B',
                despike_data=True,
                subtract_median=True):
    """
    """
    logger.info('processing {}'.format(hdf_fname))
    df_raw, header = read_hdf(hdf_fname, source_key)
    df = df_raw[['B_X', 'B_Y']] * 1e-9
    df_filtered = process_df(df,
                             despike_data=despike_data,
                             subtract_median=subtract_median)
    write_hdf(hdf_fname, df_filtered, key, header)
    return hdf_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Apply preprocessing steps to raw magnetometer data.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('hdf_fnames',
                        type=str,
                        nargs='*',
                        metavar='hdf_fname',
                        help='HDF file record to process')
    parser.add_argument('--source-key',
                        '-s',
                        type=str,
                        default='B_raw',
                        help='')
    parser.add_argument('--key',
                        '-k',
                        type=str,
                        default='B',
                        help='key to associate with the processed records')
    parser.add_argument('--he',
                        action='store_true',
                        help='include results in HE coordinate')
    args = parser.parse_args(argv[1:])

    for hdf_fname in args.hdf_fnames:
        process(hdf_fname,
                source_key=args.source_key,
                key=args.key,
                he=args.he)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
