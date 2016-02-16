from __future__ import division

import logging
import math
import itertools
from datetime import datetime, timedelta

import numpy as NP
import scipy.signal
import pandas as PD
import pynfftls
from intervals import DateTimeInterval

from repository import get_data_frame, get_root
from constants import LUNAR_SYNODIC_MONTH, SOLAR_TROPICAL_YEAR
from ..util.date import toJ2000
from ..util.search import find_le, find_ge, find_gt

logger = logging.getLogger('pyrsss.mag.love_gannon')


"""
J. J. Love and J. L. Gannon, "Revised Dst and the epicycles of
magnetic disturbance: 1958--2009," Ann. Geophys., vol 27,
pp. 3101--3131, 2009.
"""


SAMPLES = {'min': 60 * 24}
"""
Mapping between the number of samples per day for a given sampling
period.
"""


def quiet_days(series,
               n_days=5,
               period='min'):
    """
    Given a pandas time series *series*, return a mapping between
    months and the *n_days* quietest days for that month. See (4) of
    Love and Gannon. The parameter *period* gives the sampling period
    (`'min'` = 1 measurement per minute).
    """
    quiet_day_map = {}
    delta_H_i = series.diff().abs().\
                groupby(PD.TimeGrouper(freq='D')).\
                filter(lambda x: x.isnull().sum() <= int(0.5 * SAMPLES[period])).\
                groupby(PD.TimeGrouper(freq='D')).mean()
    for month, delta_H_i_month in delta_H_i.groupby(PD.TimeGrouper(freq='M')):
        quiet_day_map[month] = delta_H_i_month.nsmallest(n_days).sort_index().index
    return quiet_day_map


def chebyshev_fit(series,
                  quiet_day_map,
                  deg=10):
    """
    ???
    """
    grouped_by_day = series.groupby(PD.TimeGrouper(freq='D'))
    x = []
    y = []
    for month_end in sorted(quiet_day_map):
        for quiet_day in quiet_day_map[month_end]:
            quiet_day_H = grouped_by_day.get_group(quiet_day)
            quiet_day_H = quiet_day_H[quiet_day_H.notnull()]
            x.extend([toJ2000(time_stamp.to_datetime()) for time_stamp in quiet_day_H.index])
            y.extend(quiet_day_H.values)
    return NP.polynomial.chebyshev.chebfit(x, y, deg, full=True)


def external_field_variation(series, coeffs):
    """
    ???
    """
    x = map(toJ2000, series.index.to_datetime())
    return PD.Series(index=series.index,
                     data=series.values - NP.polynomial.chebyshev.chebval(x, coeffs))


def slice_series_by_dt(series,
                       min_dt,
                       max_dt):
    """
    ???
    """
    i1, _ = find_le(series.index, PD.Timestamp(min_dt))
    i2, _ = find_ge(series.index, PD.Timestamp(max_dt))
    return series[i1:i2]


def remove_active_days(E_series,
                       min_threshold,
                       quiet_min,
                       quiet_max,
                       delta=timedelta(hours=4)):
    """
    ???
    """
    Q_series = E_series.copy()
    active_intervals = []
    while Q_series.min() < min_threshold:
        dt_min = Q_series.argmin().to_datetime()
        Q_series_min = Q_series.index[0].to_datetime()
        Q_series_max = Q_series.index[-1].to_datetime()
        active_interval = DateTimeInterval([max(dt_min - timedelta(hours=12),
                                               Q_series_min),
                                           min(dt_min + timedelta(hours=12),
                                               Q_series_max)])
        # expand left
        while True:
            Q_series_left = slice_series_by_dt(Q_series,
                                               active_interval.lower,
                                               active_interval.lower + delta)
            if Q_series_left.min() > quiet_min and Q_series_left.max() < quiet_max:
                break
            active_interval = DateTimeInterval([max(active_interval.lower - delta,
                                                   Q_series_min),
                                               active_interval.upper])
            if active_interval.lower == Q_series_min:
                break
        # expand right
        while True:
            Q_series_right = slice_series_by_dt(Q_series,
                                                active_interval.upper - delta,
                                                active_interval.upper)
            if Q_series_right.min() > quiet_min and Q_series_right.max() < quiet_max:
                break
            active_interval = DateTimeInterval([active_interval.lower,
                                               min(active_interval.upper + delta,
                                                   Q_series_max)])
            if active_interval.upper == Q_series_max:
                break
        # remove active period from
        Q_series = PD.concat([slice_series_by_dt(Q_series,
                                                 Q_series_min,
                                                 active_interval.lower - timedelta(seconds=1)),
                              slice_series_by_dt(Q_series,
                                                 active_interval.upper + timedelta(seconds=1),
                                                 Q_series_max)])
        active_intervals.append(active_interval)
    return Q_series, active_intervals


def find_harmonic_peaks(f_range,
                        q_series_dft,
                        threshold=10,
                        center_names=['1', '2', '4', '1/2', '0.8', '4/3', '2/3', '4/7'],
                        I_year_harmonics=[-4, -2, -1, 0, 1, 2, 4],
                        I_month_harmonics=[-6, -4, -2, -1, 0, 1, 2, 4, 6],
                        BETA=1,
                        DELTA=1000):
    """
    """
    # convert center name to center period
    def eval_fraction(s):
        tokens = s.split('/')
        if len(tokens) == 1:
            return float(tokens[0])
        elif len(tokens) == 2:
            return float(tokens[0]) / float(tokens[1])
        else:
            raise ValueError('could not convert {} to float'.format(s))
    # compute relative difference between inputs
    def relative_pct_difference(x, y):
        return (x - y) / y * 100
    # only consider (strictly) positive frequencies
    I = f_range > 0
    f_range = f_range[I]
    q_series_dft = q_series_dft[I]
    # expand combinations of harmonics
    year_harmonics = [i * SOLAR_TROPICAL_YEAR for i in I_year_harmonics]
    month_harmonics = [i * LUNAR_SYNODIC_MONTH for i in I_month_harmonics]
    harmonics = list(itertools.product(year_harmonics, month_harmonics))
    year_names = [str(i) + 'Y' if i != 0 else '' for i in I_year_harmonics]
    month_names = [str(i) + 'M' if i != 0 else '' for i in I_month_harmonics]
    harmonic_names = list(itertools.product(year_names, month_names))
    # compute magnitudes and background level at harmonic frequency
    # locations
    center_periods = map(eval_fraction, center_names)
    magnitudes = {}
    for center_period, center_name in zip(center_periods, center_names):
        for (y_harmonic, m_harmonic), (y_name, m_name) in zip(harmonics, harmonic_names):
            if y_harmonic == m_harmonic == 0:
                continue
            # build name
            name = center_name + 'D'
            if y_name:
                name += ' ' + y_name
            if m_name:
                name += ' ' + m_name
            # compute frequency from period
            freq = 1/center_period
            if y_harmonic != 0:
                freq += 1/y_harmonic
            if m_harmonic != 0:
                freq += 1/m_harmonic
            # compute magnitude and background at frequency
            i, _ = find_gt(f_range, freq)
            i1 = max(0, i - BETA)
            i2 = min(len(q_series_dft), i + BETA)
            magnitude = math.log10(NP.max(NP.abs(q_series_dft[i1:i2+1])))
            j1 = max(0, i - DELTA)
            j2 = min(len(q_series_dft), i + DELTA)
            background = math.log10(NP.median(NP.abs(q_series_dft[j1:j2+1])))
            # print(relative_difference(magnitude, background))
            if relative_pct_difference(magnitude, background) >= threshold:
                magnitudes[name] = (freq, magnitude, background)
    return magnitudes


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    # root = get_root()

    # logger.info('loading')
    # dataframe, _ = get_data_frame(root,
    #                               'frd',
    #                               0,
    #                               columns=['time', 'H'])
    #                               # years=[2002])

    # series = PD.Series(index=dataframe.time,
    #                    data=dataframe.H)

    # quiet_day_map = quiet_days(series)

    # logger.info('fitting')
    # c, info = chebyshev_fit(series, quiet_day_map)

    # logger.info('evaluating')
    # E_series = external_field_variation(series, c)

    # E_series.to_hdf('/tmp/E.hdf', 'E_frd', mode='w', complevel=9, complib='zlib')

    # import pylab as PL

    # fig = PL.figure()

    # PL.plot_date(E_series.index.to_datetime(),
    #              E_series.values,
    #              c='r',
    #              mec='None',
    #              ms=0.5)

    # PL.show()

    E_series = PD.read_hdf('/tmp/E.hdf')
    E_series = E_series[E_series.notnull()]

    summary = E_series.describe()

    # print(summary)
    # assert False

    (Q_series, active_intervals) = remove_active_days(E_series,
                                                      -100,
                                                      -20,
                                                       20)

    from cPickle import dump

    with open('Q_series.pkl', 'w') as fid:
        dump((Q_series, active_intervals), fid, -1)

    import ipdb; ipdb.set_trace()

    # print(dir(summary))

    # print(E_series.argmin())


    hist, bin_edges = NP.histogram(E_series.values,
                                   bins=1024,
                                   normed=True)

    delta = bin_edges[1] - bin_edges[0]
    bin_centers = [x + delta / 2 for x in bin_edges[:-1]]

    cdf = NP.cumsum(hist * delta)

    print(cdf[-1])

    # assert False

    # print(E_series.min())
    # print(E_series.max())

    import pylab as PL

    PL.plot(bin_centers,
            cdf)

    PL.autoscale(axis='x', tight=True)

    PL.show()

    # print(len(E_series))
    # print(E_series.min())
    # print(E_series.max())

    # import pylab as PL
    # PL.hist(E_series.)


    # N = 1024
    # days_to_seconds = 60 * 60 * 24

    # periods = NP.linspace(0.9 * days_to_seconds, 1.1 * days_to_seconds, N)
    # freqs = 1 / periods
    # angular_freqs = 2 * math.pi * freqs

    # # time = [toJ2000(x.to_datetime()) for x in E_series.index]
    # time = E_series.index.astype(NP.int64) / 1e9
    # time -= time[0]

    # # print(time[:10])

    # # print(E_series.index[0])
    # # print(toJ2000(E_series.index[0].to_datetime()))
    # # print(toJ2000(E_series.index[1].to_datetime()))
    # # # print(dir(E_series.index[0]))
    # # # print(type(E_series.index[0]))
    # # assert False

    # # pgram = scipy.signal.lombscargle(time, E_series.values, angular_freqs)
    # # normalized_pgram = NP.sqrt(4 * (pgram / len(E_series)))

    # ofac = 2
    # hifac = 1

    # # print('begin nfftls')
    # # (f_nfftls, p_nfftls) = pynfftls.period(time, E_series.values, ofac, hifac)
    # # print('end nfftls')

    # # from cPickle import dump
    # # with open('/tmp/love_gannon.pkl', 'w') as fid:
    # #     dump((f_nfftls, p_nfftls), fid, -1)

    # from cPickle import load
    # with open('/tmp/love_gannon.pkl') as fid:
    #     (f_nfftls, p_nfftls) = load(fid)

    # import pylab as PL

    # fig = PL.figure()
    # # PL.plot(1 / f_nfftls / days_to_seconds,
    # #         p_nfftls)
    # PL.semilogy(1 / f_nfftls / days_to_seconds,
    #          p_nfftls)
    # # PL.autoscale(axis='x', tight=True)
    # # PL.xlim(0.9, 1.1)

    # PL.show()
