from __future__ import division

import sys
import logging
import warnings
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import timedelta
from itertools import groupby
from collections import OrderedDict

import numpy as NP
import pandas as PD
import scipy.signal

from iaga2hdf import read_hdf, write_hdf
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


def minute_interval_filter_firwin():
    """
    Synthesize and return impulse response of minute interval filter
    for analysis of magnetometer data. Design a filter using the
    window method. The filter is designed to stop from 0 to 0.75 mHz
    and pass from 1 mHz to Nyquist.

    For reference: this is the filter used to process the data for the
    Space Weather 2017 paper.

    fs_1m = 1./60
    bands_1m = [0, 1e-3 - 0.25e-3, 1e-3, fs_1m/2]
    desired_1m = [0, 1]

    _ = fir_response(h_1m, bands_1m, desired_1m, Hz=fs_1m)

    FIR filter is linear phase type 1

    band 0 abs deviations from 0:
    min = 4.466e-05  (db=-87.000846)
    med = -4.218e-03
    max = 3.027e-02 (db=-30.379491

    band 1 abs deviations from 1:
    min = 4.412e-06  (db=-107.106525)
    med = 2.398e-03
    std = 4.059e-03
    max = 2.982e-02 (db=-30.509428
    """
    Ts = 60
    fs = 1. / Ts
    fn = fs / 2
    width = 0.25e-3
    cutoff = 1e-3
    ripple = 30
    numtaps, beta = scipy.signal.kaiserord(ripple, width / fn)
    if numtaps % 2 == 0:
        numtaps += 1
    return scipy.signal.firwin(numtaps,
                               cutoff - width/2,
                               width=width,
                               pass_zero=False,
                               nyq=fn)


def minute_interval_filter(N_remez=201):
    """
    Synthesize and return impulse response of minute interval filter
    for analysis of magnetometer data. Design a length *N_remez*
    min-max optimal filter. The filter is designed to stop from 0 to
    0.75 mHz and pass from 1 mHz to Nyquist.

    fs_1m = 1./60
    bands_1m = [0, 1e-3 - 0.25e-3, 1e-3, fs_1m/2]
    desired_1m = [0, 1]

    _ = fir_response(h_1m, bands_1m, desired_1m, Hz=fs_1m)

    FIR filter is linear phase type 1

    band 0 abs deviations from 0:
    min = 1.487e-05  (db=-96.553225)
    med = -1.102e-03
    max = 1.569e-03 (db=-56.089238

    band 1 abs deviations from 1:
    min = 1.018e-06  (db=-119.848775)
    med = -1.066e-05
    std = 4.820e-04
    max = 1.571e-03 (db=-56.074809
    """
    return scipy.signal.remez(N_remez,
                              [0, 0.75e-3, 1e-3, 1./60/2],
                              [0, 1],
                              [1, 1],
                              Hz=1./60)


def second_interval_filter(N_remez=4001):
    """
    Synthesize and return impulse response of second interval filter
    for analysis of magnetometer data. Design a length *N_remez*
    min-max optimal filter. The filter is designed to stop from 0 to
    0.25 mHz, pass from 1 to 100 mHz, and stop from 101 mHz to
    Nyquist.


    fs_1m = 1.
    bands_1s = [0, 1e-3 - 0.25e-3,
                1e-3, 100e-3,
                101e-3, fs_1s/2]
    desired_1s = [0, 1, 0]

    _ = fir_response(h_1m, bands_1m, desired_1m, Hz=fs_1m)

    FIR filter is linear phase type 1

    band 0: 0 Hz -- 750 \muHz
    abs deviations from 0 statistics:
    min = 6.872e-05  (db=-83.257739)
    med = -4.906e-02
    max = 7.486e-01 (db=-2.514808)

    band 1: 1 mHz -- 100 mHz
    abs deviations from 1 statistics:
    min = 2.323e-07  (db=-132.678419)
    med = -6.093e-07
    std = 5.871e-04
    max = 2.349e-03 (db=-52.584106)

    band 2: 101 mHz -- 500 mHz
    abs deviations from 0 statistics:
    min = 1.430e-07  (db=-136.893586)
    med = -1.353e-03
    max = 2.853e-03 (db=-50.895031)
    """
    return scipy.signal.remez(N_remez,
                              [0, 0.25e-3, 1e-3, 100e-3, 101e-3, 0.5],
                              [0, 1, 0],
                              [1, 1, 1],
                              Hz=1)


def consecutive_nans(x, y):
    """
    Return the maximum number of nans encountered in the logical or of
    the *x* and *y* time series.
    """
    nan = NP.isnan(x) | NP.isnan(y)
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
        df = PD.DataFrame(index=dt,
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
        Bx_filtered -= NP.mean(Bx_filtered)
        By_filtered -= NP.mean(By_filtered)
    # build DataFrame and store to disk
    return PD.DataFrame(index=dt_filtered,
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
        dt_diff = NP.diff(df.index.values)
        delta_timedelta64 = min(dt_diff)
        delta_seconds = delta_timedelta64 / NP.timedelta64(1, 's')
        delta = timedelta(seconds=delta_seconds)
    logger.info('Using delta = {} (s)'.format(delta.total_seconds()))
    index_new = PD.date_range(start=df.index[0],
                              end=df.index[-1],
                              freq=delta)
    missing = sorted(set(index_new) - set(df.index))
    if missing:
        logger.warning('Missing time indices (filled by NaNs):')
        for x in missing:
            logger.warning(x)
    return df.reindex(index_new, copy=False), delta


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
    df_filtered = PD.DataFrame(index=dt_filtered,
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
