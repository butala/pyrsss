import logging

import numpy as NP
import pandas as PD

from repository import get_data_frame, get_root
from ..util.date import toJ2000


"""
J. J. Love and J. L. Gannon, "Revised Dst and the epicycles of
magnetic disturbance: 1958--2009", Ann. Geophys., vol 27,
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
        quiet_day_map[month] = delta_H_i_month.nsmallest(n_days).index
    return quiet_day_map


def chebyshev_fit(series,
                  quiet_day_map,
                  deg=20):
    """
    ???
    """
    fit_map = {}
    grouped_by_day = series.groupby(PD.TimeGrouper(freq='D'))
    for month_end in sorted(quiet_day_map):
        x = []
        y = []
        for quiet_day in quiet_day_map[month_end]:
            quiet_day_H = grouped_by_day.get_group(quiet_day)
            quiet_day_H = quiet_day_H[quiet_day_H.notnull()]
            x.extend([toJ2000(time_stamp.to_datetime()) for time_stamp in quiet_day_H.index])
            y.extend(quiet_day_H.values)
        fit_map[month_end] = NP.polynomial.chebyshev.chebfit(x, y, deg, full=True)
    return fit_map


def external_field_variation(series,
                             fit_map):
    """
    ???
    """
    E_series_list = []
    for month, group in series.groupby(PD.TimeGrouper(freq='M')):
        # NEED TO DO FIT QUALITY CONTROL!!!
        coeffs, _ = fit_map[month]
        j2000 = [toJ2000(time_stamp.to_datetime()) for time_stamp in group.index]
        poly_val = NP.polynomial.chebyshev.chebval(j2000, coeffs)
        # E_group = group.values - poly_val
        E_group = poly_val
        E_series_list.append(PD.Series(index=group.index, data=E_group))
    return PD.concat(E_series_list)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    root = get_root()

    dataframe, _ = get_data_frame(root,
                                  'frd',
                                  0,
                                  columns=['time', 'H'],
                                  years=[2002])
    series = PD.Series(index=dataframe.time,
                       data=dataframe.H)

    quiet_day_map = quiet_days(series)

    fit_map = chebyshev_fit(series, quiet_day_map)

    E_series = external_field_variation(series,
                                        fit_map)

    # print(len(series))
    # print(len(E_series))

    import pylab as PL

    fig = PL.figure()

    PL.plot_date(series.index.to_datetime(),
                 series.values,
                 mec='None',
                 ms=1)

    PL.plot_date(E_series.index.to_datetime(),
                 E_series.values,
                 c='r',
                 mec='None',
                 ms=0.5)

    PL.show()
