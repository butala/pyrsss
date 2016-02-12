from __future__ import division

import logging
import math

import numpy as NP
import scipy.signal
import pandas as PD
import pynfftls

from repository import get_data_frame, get_root
from ..util.date import toJ2000

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


def remove_active_days(E_series):
    """
    ???
    """
    pass


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

    N = 1024
    days_to_seconds = 60 * 60 * 24

    periods = NP.linspace(0.9 * days_to_seconds, 1.1 * days_to_seconds, N)
    freqs = 1 / periods
    angular_freqs = 2 * math.pi * freqs

    # time = [toJ2000(x.to_datetime()) for x in E_series.index]
    time = E_series.index.astype(NP.int64) / 1e9
    time -= time[0]

    # print(time[:10])

    # print(E_series.index[0])
    # print(toJ2000(E_series.index[0].to_datetime()))
    # print(toJ2000(E_series.index[1].to_datetime()))
    # # print(dir(E_series.index[0]))
    # # print(type(E_series.index[0]))
    # assert False

    # pgram = scipy.signal.lombscargle(time, E_series.values, angular_freqs)
    # normalized_pgram = NP.sqrt(4 * (pgram / len(E_series)))

    ofac = 2
    hifac = 1

    # print('begin nfftls')
    # (f_nfftls, p_nfftls) = pynfftls.period(time, E_series.values, ofac, hifac)
    # print('end nfftls')

    # from cPickle import dump
    # with open('/tmp/love_gannon.pkl', 'w') as fid:
    #     dump((f_nfftls, p_nfftls), fid, -1)

    from cPickle import load
    with open('/tmp/love_gannon.pkl') as fid:
        (f_nfftls, p_nfftls) = load(fid)

    import pylab as PL

    fig = PL.figure()
    # PL.plot(1 / f_nfftls / days_to_seconds,
    #         p_nfftls)
    PL.semilogy(1 / f_nfftls / days_to_seconds,
             p_nfftls)
    # PL.autoscale(axis='x', tight=True)
    # PL.xlim(0.9, 1.1)

    PL.show()
