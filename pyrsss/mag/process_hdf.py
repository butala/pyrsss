from __future__ import division

import sys
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from itertools import groupby

import numpy as NP
import pandas as PD
import scipy.signal

from iaga2hdf import read_hdf, write_hdf
from ..util.nan import nan_interp
from ..util.signal import lp_fir_filter

logger = logging.getLogger('pyrsss.mag.process_hdf')


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
    """
    DOCUMENT 4001 CHOICE!!!
    """
    N_remez = 4001
    return scipy.signal.remez(N_remez,
                              [0, 0.25e-3, 1e-3, 100e-3, 101e-3, 0.5],
                              [0, 1, 0],
                              [1, 1, 1],
                              Hz=1)



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



def process_timeseries(dt,
                       Bx,
                       By,
                       c1='Bx',
                       c2='By',
                       remove_mean=True):
    """
    """
    n = consecutive_nans(Bx, By)
    interval = (dt[1] - dt[0]).total_seconds()
    if n > 0:
        logger.warning('longest contiguous gap = {:.2f} minutes'.format(n * interval / 60))
    # fill data gaps via linear interpolation
    Bx = nan_interp(Bx)
    By = nan_interp(By)
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
            remove_mean=True):
    """
    ???
    """
    logger.info('processing {}'.format(hdf_fname))
    df_raw, header = read_hdf(hdf_fname, source_key)
    dt = [PD.to_datetime(x).to_pydatetime() for x in df_raw.index.values]
    Bx_raw = df_raw['Bx'].values * 1e-9
    By_raw = df_raw['By'].values * 1e-9
    df_filtered = process_timeseries(dt,
                                     Bx_raw,
                                     By_raw,
                                     remove_mean=remove_mean)
    if he:
        Bh_raw = df_raw['Bh'].values * 1e-9
        Be_raw = df_raw['Be'].values * 1e-9
        df_he_filtered = process_timeseries(dt,
                                            Bh_raw,
                                            Be_raw,
                                            c1='Bh',
                                            c2='Be',
                                            remove_mean=remove_mean)
        df_filtered = df_filtered.join(df_he_filtered)
    write_hdf(hdf_fname, df_filtered, key, header)
    return hdf_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('',
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
