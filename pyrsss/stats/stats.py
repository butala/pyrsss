import logging
import sys
import math
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections.abc import Iterable

import scipy.stats
import numpy as np

logger = logging.getLogger('pyrsss.stats.stats')


def mad(l):
    """
    Compute the median absolute deviation (a robust measure of spread)
    of the list of values *l*.
    """
    median = np.median(l)
    return np.median(np.abs(l - median)), median


def robust_std(l, alpha=1/scipy.stats.norm.ppf(0.75)):
    """
    Compute a robust estimate of the standard deviation for the list
    of values *l* (by default, for normally distributed samples ---
    see https://en.wikipedia.org/wiki/Median_absolute_deviation).
    """
    return alpha * mad(l)[0]


def max_consectuive_true(df):
    """
    Return the tuple of the number of `True` elements in each column
    of *df* and the largest consecutive run of `True` elements.
    """
    # see https://stackoverflow.com/questions/52717996/how-can-i-count-the-number-of-consecutive-trues-in-a-dataframe
    cumsum = df.cumsum()
    c = cumsum.sub(cumsum.mask(df).ffill().fillna(0)).astype(int)
    mask = df.any()
    return cumsum.iloc[-1].tolist(), np.where(mask, c.max(), -1).tolist()


def despike(df, window=31, l=6):
    """
    Remove outliers from the columns of :class:`DataFrame` by
    comparing the absolute deviation from the windowed median to the
    windowed robust standard deviation (see :func:`robust_std`). Use a
    centered window of length *window* (must be odd). Replace values
    that are *l* robust standard deviations from the absolute
    difference from the median with the median.

    Reference: Hampel F. R., "The influence curve and its role in
    robust estimation," Journal of the American Statistical
    Association, 69, 382-393, 1974.
    """
    if window % 2 == 0:
        raise ValueError('window length must be odd')
    df_rolling = df.rolling(window, center=True)
    df_rolling_median = df_rolling.median()
    df_robust_std = df_rolling.apply(robust_std, raw=True)
    # it may be possible to speed this up:
    # https://stackoverflow.com/questions/42865103/numpy-version-of-rolling-mad-mean-absolute-deviation
    I = (df - df_rolling_median).abs() > l * df_robust_std
    cumsum, consecutive = max_consectuive_true(I)
    for column, cumsum_i, consecutive_i in zip(df.columns, cumsum, consecutive):
        if cumsum_i > 0:
            logger.info('column {}: {} points ({:.2f}%) removed --- max run removed = {} points'.format(column, cumsum_i, cumsum_i / len(df) * 100, consecutive_i))
    df_despike = df.copy()
    df_despike[I] = df_rolling_median
    return df_despike.iloc[(window-1):-(window-1)]


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.

    References:
    - http://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
    - https://en.wikipedia.org/wiki/Mean_square_weighted_deviation

    Note: The method is biased (division by N and not (N-1)).
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))


class Stats(object):
    def __init__(self, *args):
        self.reset()
        self(*args)

    def reset(self):
        self.mu = 0
        self.s2 = 0
        self.N = 0
        self._max = -float('inf')
        self._min = float('inf')

    @property
    def mean(self):
        return self.mu

    @property
    def var(self):
        return self.s2 / (self.N-1) if self.N >= 2 else float('nan')

    @property
    def sigma(self):
        return math.sqrt(self.var)

    @property
    def rms(self):
        return math.sqrt(self.mean**2 + self.var**2)

    @property
    def max(self):
        return self._max

    @property
    def min(self):
        return self._min

    def __call__(self, *x):
        """
        """
        if len(x) == 1 and isinstance(x[0], Iterable):
            x = x[0]
        for x_i in x:
            self.N += 1
            mu_previous = self.mu
            self.mu += (x_i - self.mu) / self.N
            self.s2 += (x_i - mu_previous) * (x_i - self.mu)
            if x_i < self._min:
                self._min = x_i
            if x_i > self._max:
                self._max = x_i
        return self

    def __repr__(self):
        return 'N={} mean={:f} sigma={:f} min={:f} max={:f}'.format(self.N,
                                                                    self.mean,
                                                                    self.sigma,
                                                                    self.min,
                                                                    self.max)


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Report mean and standard deviation for the stream of numbers read from stdin',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    args = parser.parse_args(argv[1:])

    stats = Stats()
    for line in sys.stdin:
        stats(float(line))
    print(stats)


if __name__ == '__main__':
    sys.exit(main())
