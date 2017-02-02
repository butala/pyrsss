from __future__ import division

import sys
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from itertools import groupby

import numpy as NP
import pandas as PD
import scipy.signal

from iaga2002 import parse
from ..util.nan import nan_interp
from ..util.signal import lp_fir_filter

logger = logging.getLogger('pyrsss.mag.process_iaga')


def minute_interval_filter():
    """ ??? """
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


def second_interval_filter():
    """ ??? """
    Ts = 1
    fs = 1./ Ts
    fn = fs / 2
    width = 0.1e-3

    freq = [0, 1e-3 - width, 1e-3, 100e-3, 100e-3 + width, fn]
    gain = [0, 0, 1, 1, 0, 0]

    width_1sec = 0.25e-3
    fn_1sec = 1. / 60 / 2
    ripple = 50
    numtaps, beta = scipy.signal.kaiserord(ripple, width_1sec / fn_1sec)

    # want type 3 filter (0 at DC and max frequency): odd number of
    # taps and anti-symmetric
    return scipy.signal.firwin2(numtaps * 2 + 1,
                                freq,
                                gain,
                                window=('kaiser', beta),
                                antisymmetric=True,
                                nyq=fn)


def consecutive_nans(x, y):
    """
    """
    nan = NP.isnan(x) | NP.isnan(y)
    lengths = []
    for key, group in groupby(nan):
        if key == True:
            lengths.append(len(list(group)))
    if lengths:
        return max(lengths)
    return 0


def process(hdf_fname,
            iaga_fnames,
            mode='valid',
            remove_mean=True,
            save_raw=False,
            strict=False):
    """
    ???
    """
    # parse IAGA2002 files
    headers = []
    data = []
    intervals = []
    for iaga_fname in iaga_fnames:
        logger.info('reading {}'.format(iaga_fname))
        header_i, data_i = parse(iaga_fname, strict=strict)
        intervals_i = set(NP.diff(data_i.keys()))
        if len(intervals_i) > 1:
            raise ValueError('multiple sample intervals detected in {}'.format(iaga_fanme))
        intervals.append(intervals_i.pop())
        headers.append(header_i)
        data.append(data_i)
    intervals = set(intervals)
    if len(intervals) > 1:
        raise ValueError('no common sample interval found: {}'.format(', '.join(iaga_fnames)))
    interval = intervals.pop().total_seconds()
    # gather data
    dt = []
    x = []
    y = []
    for data_i in data:
        for dt_i, col_i in data_i.iteritems():
            dt.append(dt_i)
            # 1e-9 to convert nT to T
            x.append(col_i.x * 1e-9)
            y.append(col_i.y * 1e-9)
    n = consecutive_nans(x, y)
    if n > 0:
        logger.warning('longest contiguous gap = {:.2f} minutes'.format(n * interval / 60))
    # fill data gaps via linear interpolation
    x = nan_interp(x)
    y = nan_interp(y)
    # apply 1 - 100 mHz bandpass filter
    if interval == 1.0:
        h = second_interval_filter()
    elif interval == 60.0:
        h = minute_interval_filter()
    else:
        raise ValueError('1 to 100 mHz filter not yet synthesized for {} s interval data'.format(interval))
    x_filtered = lp_fir_filter(h, x, mode='valid')
    y_filtered = lp_fir_filter(h, y, mode='valid')
    D_min = len(h)
    D_max = len(x)
    dt_filtered = dt[(D_min - 1):D_max]
    # remove mean
    if remove_mean:
        x -= NP.mean(x)
        y -= NP.mean(y)
    # build DataFrame and store to disk
    df_filtered = PD.DataFrame(index=dt_filtered,
                               data={'Bx': x_filtered,
                                     'By': y_filtered})
    df_filtered.to_hdf(hdf_fname,
                       key='filtered',
                       mode='w')
    if save_raw:
        df = PD.DataFrame(index=dt,
                          data={'Bx': x,
                                'By': y})
        df.to_hdf(hdf_fname,
                  key='raw',
                  mode='a')
    return hdf_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Process magnetometer time-series spanning IAGA input and save to HDF.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('hdf_fname',
                        type=str,
                        help='ouput HDF file name')
    parser.add_argument('iaga_fnames',
                        type=str,
                        nargs='*',
                        metavar='iaga_fname',
                        help='input IAGA2002 file')
    parser.add_argument('--mode',
                        '-m',
                        choices=['same', 'valid'],
                        default='valid')
    parser.add_argument('--raw',
                        '-r',
                        action='store_true',
                        help='include unprocessed time-series')
    args = parser.parse_args(argv[1:])

    process(args.hdf_fname,
            args.iaga_fnames,
            save_raw=args.raw,
            mode=args.mode)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
