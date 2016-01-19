from __future__ import division

import sys
import math
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


class Stats(object):
    def __init__(self):
        self.mu = 0
        self.s2 = 0
        self.N = 0

    @property
    def mean(self):
        return self.mu

    @property
    def var(self):
        return self.s2 / (self.N-1) if self.N-1 >= 2 else float('nan')

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


def main(args):
    parser = ArgumentParser('Report mean and standard deviation for the stream of numbers read from stdin',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    args = parser.parse_args(args)

    stats = Stats()
    for line in sys.stdin:
        stats(float(line))
    print(stats)


if __name__ == '__main__':
    main(sys.argv[1:])
