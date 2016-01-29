from __future__ import division

import sys
import math
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import scipy.stats
import numpy as NP


def mad(l):
    """
    Compute the median absolute deviation (a robust measure of spread)
    of the list of values *l*.
    """
    median = NP.median(l)
    return NP.median([abs(x - median) for x in l]), median


def robust_std(l, alpha=1/scipy.stats.norm.ppf(0.75)):
    """
    Compute a robust estimate of the standard deviation for the list
    of values *l* (by default, for normally distributed samples ---
    see https://en.wikipedia.org/wiki/Median_absolute_deviation).
    """
    return alpha * mad(l)[0]


class Stats(object):
    def __init__(self):
        self.reset()

    def reset(self):
        self.mu = 0
        self.s2 = 0
        self.N = 0

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

    def __call__(self, *x):
        """
        """
        for x_i in x:
            self.N += 1
            mu_previous = self.mu
            self.mu += (x_i - self.mu) / self.N
            self.s2 += (x_i - mu_previous) * (x_i - self.mu)
        return self

    def __str__(self):
        return str(self.mean) + ' ' + str(self.sigma)


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
