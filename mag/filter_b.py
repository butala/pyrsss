from __future__ import division

import sys
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import timedelta
from contextlib import closing

import numpy as NP
import scipy.signal
from tables import open_file, IsDescription, Time64Col, Float64Col

from iaga2002 import parse, fname2date
from ..util.signal import nextpow2
from ..util.date import UNIX_EPOCH


def load_data(iaga2002_fnames):
    """ ??? """
    assert len(iaga2002_fnames) == 3
    headers, records = zip(*map(parse, iaga2002_fnames))
    site = set([x.IAGA_CODE for x in headers])
    if len(site) != 1:
        raise ValueError('data must all originate from the same site')
    dates = map(fname2date, iaga2002_fnames)
    if (dates[2] - dates[0] + timedelta(days=1)).days != 3:
        raise ValueError('date must span 3 consecutive days')
    return ([dt for record in records for dt in record.iterkeys()],
            [v.x for record in records for v in record.itervalues()],
            [v.y for record in records for v in record.itervalues()],
            [v.z for record in records for v in record.itervalues()])


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
    raise NotImplementedError('VERIFY')
    Ts = 1
    fs = 1./ Ts
    fn = fs / 2
    width = 0.1e-3

    freq = [0, 1e-3 - width, 1e-3, 100e-3, 100e-3 + width, fn]
    gain = [0, 0, 1, 1, 0, 0]

    # THE VALUES BELOW ARBITRARY --- VERIFY
    width_1min = 0.25e-3
    fn_1min = 1. / 60 / 2
    numtaps, beta = scipy.signal.kaiserord(ripple, width_1min / fn_1min)

    # want type 3 filter (0 at DC and max frequency): odd number of
    # taps and anti-symmetric
    return scipy.signal.firwin2(numtaps * 2 + 1,
                                freq,
                                gain,
                                window=('kaiser', beta),
                                antisymmetric=True,
                                nyq=fn)


def filter_b(Bx, By, Bz):
    """ ??? """
    assert len(Bx) == len(By) == len(Bz)
    if len(Bx) == 3 * 60 * 24:
        h = minute_interval_filter()
    elif len(Bx) == 3 * 60 * 60 * 24:
        h = second_interval_filter()
    else:
        raise ValueError('expected 1 minute or 1 second interval data')
    K = len(h) + len(Bx) - 1
    N = nextpow2(K)
    H = NP.fft.fft(h, n=N)
    F_Bx = NP.fft.fft(Bx, n=N)
    F_By = NP.fft.fft(By, n=N)
    F_Bz = NP.fft.fft(Bz, n=N)
    M = len(Bx)
    return (NP.real(NP.fft.ifft(H * F_Bx))[:M],
            NP.real(NP.fft.ifft(H * F_By))[:M],
            NP.real(NP.fft.ifft(H * F_Bz))[:M])


""" ??? """
class Table(IsDescription):
    dt = Time64Col()
    Bx = Float64Col()
    By = Float64Col()
    Bz = Float64Col()


def output_table(dt, B, h5file, where, name, title, center_day_only=True):
    """ ??? """
    table = h5file.create_table(where,
                                name,
                                Table,
                                title)
    Bx, By, Bz = B
    if center_day_only:
        N = len(Bx)
        M = int(N / 3)
        Bx = Bx[M:2*M]
        By = By[M:2*M]
        Bz = Bz[M:2*M]
    row = table.row
    for dt_i, Bx_i, By_i, Bz_i in zip(dt, Bx, By, Bz):
        row['dt'] = (dt_i - UNIX_EPOCH).total_seconds()
        row['Bx'] = Bx_i
        row['By'] = By_i
        row['Bz'] = Bz_i
        row.append()
        table.flush()
    return table


def save_result(h5_fname,
                dt,
                filtered_B,
                raw_B=None,
                center_day_only=True):
    """ ??? """
    with closing(open_file(h5_fname,
                           mode='w',
                           title='pyrsss.mag.filter_b output')) as h5file:
        mag_group = h5file.create_group('/',
                                        'mag_data',
                                        'Magnetometer data')
        output_table(dt,
                     filtered_B,
                     h5file,
                     mag_group,
                     'filtered',
                     'Band-pass filtered magnetometer data',
                     center_day_only=center_day_only)
        if raw_B:
            output_table(dt,
                         raw_B,
                         h5file,
                         mag_group,
                         'raw',
                         'Magnetometer data',
                         center_day_only=center_day_only)
    return h5_fname


def process(h5_fname, iaga2002_fnames):
    """ ??? """
    dt, Bx, By, Bz = load_data(iaga2002_fnames)
    filtered_B = filter_b(Bx, By, Bz)
    save_result(h5_fname,
                dt,
                filtered_B,
                raw_B=(Bx, By, Bz))
    return h5_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Filter magnetic field time series, retaining frequencies between 1 and 100 [mHz].',
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            epilog='Three days of data are processed to reduce edge effects for the data retained on the central day.')
    parser.add_argument('h5_fname',
                        type=str,
                        help='output H5 file')
    parser.add_argument('iaga2002_fnames',
                        nargs=3,
                        metavar='iaga2002_fname',
                        type=str,
                        help='input daily IAGA2002')
    args = parser.parse_args(argv[1:])

    process(args.h5_fname,
            args.iaga2002_fnames)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
