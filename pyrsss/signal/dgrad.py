from __future__ import division

import numpy as NP
import scipy.signal


def differentiator(n, Hz=1):
    """
    Return linear phase impulse response for a length *n* filter that
    approximates the differential operator. The sampling frequency is
    *Hz*.

    The filter length *n* must be even. The remez function returns
    type 3 (for *n* odd) and 4 (for *n* even) linear phase
    filters. However, type 3 linear phase filters are 0 at $\omega =
    0$ and $\omega = \pi$.
    """
    if n % 2 == 1:
        raise ValueError('the filter length n must be even')
    return scipy.signal.remez(n,
                              [0, Hz / 2],
                              [1],
                              Hz=Hz,
                              type='differentiator') * Hz * 2 * NP.pi


if __name__ == '__main__':
    from collections import OrderedDict

    import numpy as NP
    import pylab as PL

    from pyrsss.signal.spectrum import spectrum

    l = [2, 4, 10]
    h_map = OrderedDict()
    for l_i in l:
        h_map[l_i] = differentiator(l_i)

    oversample = 8
    H_map = OrderedDict()
    for l_i, h_i in h_map.iteritems():
        H_map[l_i] = spectrum(h_i, oversample=oversample)

    f_ideal = NP.linspace(0, 0.5, 129)
    H_ideal = NP.abs(f_ideal) * 2 * NP.pi

    PL.figure(figsize=(8, 4))
    PL.plot(f_ideal,
            H_ideal,
            c='C0',
            zorder=10,
            label='l=$\infty$')
    for i, (l_i, (H_i, f_i)) in enumerate(H_map.iteritems(), 1):
        PL.plot(f_i,
                NP.abs(H_i),
                c='C{}'.format(i),
                label='l={}'.format(l_i))
    PL.legend()
    PL.xlabel('Frequency (Hz)')
    PL.ylabel('Amplitude')
    PL.title('Comparison of Differentiator Magnitude Responses')

    PL.show()
