from __future__ import division

import logging
from itertools import repeat

import numpy as NP
import scipy as SP
import scipy.signal

from spectrum import nextpow2

logger = logging.getLogger('pyrsss.signal.lfilter')


def lp_fir_type(h):
    """
    Determine if FIR filter impulse response *h* is symmetric or
    antisymmetric. Return {1, 2, 3, 4} depending on FIR filter type or
    None if the FIR filter is not linear phase.
    """
    M = len(h) - 1
    n_range = range(M + 1)
    if M % 2 == 0:
        if all([NP.isclose(h[n], h[M - n]) for n in n_range]):
            return 1
        elif all([NP.isclose(h[n], -h[M - n]) for n in n_range]):
            return 3
        else:
            return None
    else:
        if all([NP.isclose(h[n], h[M - n]) for n in n_range]):
            return 2
        elif all([NP.isclose(h[n], -h[M - n]) for n in n_range]):
            return 4
        else:
            return None
    assert False


def lp_fir_filter(h, x, real=True, mode='same', index=None):
    """
    Apply a linear phase FIR filter with impulse response *h* to the
    signal *x* and return the output (with same length as *x*) after
    compensating for the constant group delay. If *real*, return only
    the real part of the filter output.

    The *mode* parameter specifies the portion of the convolution
    returned. If `same` the output will be the same shape as *x*. If
    `full` the entire convolution is returned (`len(h) + len(x) -
    1`). Finally if mode is 'valid', return only that portion for
    which *h* and *x* completely overlap (i.e., the portion where no 0
    boundary values are included).

    If *index* is provided, return the slice of *index* corresponding
    to *mode*. The purpose is to associate the correct time indices
    with the filter output.
    """
    lp_type = lp_fir_type(h)
    if lp_type is None:
        raise ValueError('FIR filter is not linear phase')
    if lp_type in [2, 4]:
        logger.warning('linear phase FIR filter is type {} --- cannot compensate for half sample delay (compensating for integer portion only)')
    N = len(x)
    K = len(h)
    L = N + K - 1
    L_fft = nextpow2(L)
    H = NP.fft.fft(h, n=L_fft)
    X = NP.fft.fft(x, n=L_fft)
    y = NP.fft.ifft(H * X)
    if real:
        y = NP.real(y)
    M = len(h) - 1
    if mode == 'same':
        y_out = y[int(M/2):int(M/2) + N]
        if index is not None:
            index_out = index
    elif mode == 'valid':
        D_min = min(N, K)
        D_max = max(N, K)
        y_out = y[(D_min - 1):D_max]
        if index is not None:
            index_out = index[(D_min - int(M/2) - 1):(D_max - int(M/2))]
    elif mode == 'full':
        y_out = y[:L]
        if index is not None:
            delta = index[1] - index[0]
            J = int(M/2)
            index_out = [index[0] - (J-i) * delta for i in range(1, J + 1)] + \
                        index + \
                        [index[-1] + i * delta for i in range(1, J + 1)]
    else:
        raise ValueError('unknown convolution mode {} (choices are same, valid, or full)')
    if index is not None:
        return y_out, index_out
    else:
        return y_out


def fir_response(h, bands, desired, Hz=1, names=None, verbose=True):
    """
    Report on the frequency response magnitude characteristics of the
    filter with impulse response *h*. The list-like *bands* is twice
    the length of *desired*, i.e., one filter desired magnitude per
    two value frequency band specification, see
    :func:`scipy.signal.remez`). The parameter *Hz* is the sample rate
    and *names*, if given, associates a string name with each
    band. Return a mapping from band identifiers to absolute desired
    to filter magnitude per band :class:`Stats` and the median
    value. If *vervose*, dump a report to stdout.
    """
    # check for argument consistency
    bands = NP.array(bands)
    if len(bands) % 2 != 0:
        raise ValueError('# of band boundaries must be even')
    if len(bands) / 2 != len(desired):
        raise ValueError('# of band boundaries must equal to twice the length of desired (i.e., 2 boundaries per band and 1 desired amplitude per band)')
    if any(bands < 0):
        raise ValueError('no value of bands may be negative')
    if any(bands > Hz / 2):
        raise ValueError('no value of bands may be larger than the Nyquist rate (Hz / 2)')
    if names:
        if len(names) != len(desired):
            raise ValueError('there should be an equal number of bands as names')
    # compute filter frequency response
    L = nextpow2(10 * len(h))
    H = NP.fft.rfft(h, n=L)
    f = NP.fft.rfftfreq(L, 1/Hz)
    # gather stats per band
    band_stats = OrderedDict()
    medians = OrderedDict()
    band_tuples = zip(bands[::2], bands[1::2])
    for index, ((b1, b2), d) in enumerate(zip(band_tuples, desired)):
        I = [i for  i, f_i in enumerate(f) if b1 <= f_i < b2]
        diff = d - NP.abs([H[i] for i in I])
        band_stats[index] = Stats(*NP.abs(diff)), NP.median(diff)
        if names:
            band_stats[names[index]] = band_stats[index]
    if verbose:
        # report on linear phase determined from symmetry of h
        fir_type = lp_fir_type(h)
        if fir_type:
            print('FIR filter is linear phase type {}'.format(fir_type))
        else:
            print('FIR filter is NOT linear phase')
        print('')
        # output report
        for index, d in enumerate(desired):
            b1, b2 = band_tuples[index]
            b1_str = sistr(b1, 'Hz')
            b2_str = sistr(b2, 'Hz')
            stats, medians = band_stats[index]
            if names:
                print('{} ({}): {} -- {}'.format(names[index],
                                                 index,
                                                 b1_str,
                                                 b2_str))
            else:
                print('band {}: {} -- {}'.format(index,
                                                 b1_str,
                                                 b2_str))
            print('abs deviations from {} statistics:'.format(d))
            print('min = {:.3e}  (db={:f})'.format(stats.min,
                                                   20 * math.log10(stats.min)))
            print('med = {:.3e}'.format(medians))
            if d == 1:
                print('std = {:.3e}'.format(stats.sigma))
            print('max = {:.3e} (db={:f})'.format(stats.max,
                                                  20 * math.log10(stats.max)))
    return band_stats


def miso_lfilter(b, a, x, zi=None):
    """
    Multiple input, single output (MISO) IIR filter. The transfer
    function numerator terms are stored in the list *b* and the
    denominator terms in the list *a* (the length of *a* and *b* must
    be equal). The multiple inputs are stored in the list *x* (same
    length as *a* and *b* and each element of *x* must be the same
    length). The initial filter delay initial conditions are stored in
    *zi* (use 0 initial condition if not provided). If *zi* is given,
    then return the tuple containing the filter output and the final
    delay values. Otherwise, return just the filter output.

    This routine produces the same output as the Matlab routine *sim*.
    """
    assert len(b) == len(a) == len(x)

    N_set = set(map(len, x))
    assert len(N_set) == 1

    if zi is None:
        zi = repeat(zi)
        return_zf = False
    else:
        assert len(zi) == len(b)
        return_zf = True

    y_i = [SP.signal.lfilter(b_i, a_i, x_i, zi=zi_i) for b_i, a_i, x_i, zi_i in zip(b, a, x, zi)]

    if return_zf:
        y_i, zf = zip(*y_i)
        return NP.sum(y_i, axis=0), zf
    else:
        return NP.sum(y_i, axis=0)


def difference_eq(b, a, x, prior=None):
    """
    A simple implementation of a linear, constant-coefficient,
    difference equation calculator and a crude reimplementation of
    `scipy.signal.lfilter` for testing purposes. The values *b* and
    *a* are the MA and AR filter parameters and *x* is the input
    sequence. If *prior* is given, it is the tuple containing the
    prior input and output values (otherwise a initial conditions of
    0s are used).
    """
    N = max(len(b), len(a)) - 1
    b = NP.concatenate((b, NP.zeros(N + 1 - len(b))))
    a = NP.concatenate((a, NP.zeros(N + 1 - len(a))))
    if prior:
        xn, yn = prior
        assert len(xn) == len(yn) == N
    else:
        xn = NP.zeros(N)
        yn = NP.zeros(N)
    Nx = len(x)
    x = NP.concatenate((xn, x))
    y = NP.concatenate((yn, NP.empty(Nx)))
    a_rev = a[1:][::-1]
    b_rev = b[::-1]
    for i in range(N, N + Nx):
        y[i] = -NP.dot(a_rev, y[i-N:i]) + NP.dot(b_rev, x[i-N:i+1])
    return y[N:]


def directform2Tzi(b, a, x, y):
    """
    Compute transposed direct form 2 (the structure used by the
    :func:`scipy.signal.filter`) initial conditions. The coefficients
    *a* and *b* are the AR and MA filter parameters. The filter prior
    input and output are given in *x* and *y* indexed using the
    convention that x[-1] (in the Python indexing convention) is the
    input at time n=-1 (mathematical convention), x[-2] (in the Python
    indexing convention) is the input at time n=-2 (mathematical
    convention), and so on.

    In particular, if

    zi = directform2Tzi(b, a, xn, yn)
    y_lfilter = SP.signal.lfilter(b, a, x, zi=zi)[0]

    then

    difference_eq(b, a, x, prior=(xn, yn)) == y_lfilter
    """
    assert NP.allclose(a[0], 1)
    assert len(x) == len(y) == max(len(b), len(a)) - 1
    N = len(x)
    if N == 0:
        return []
    z = NP.empty(N)
    b = NP.concatenate((b, NP.zeros(N + 1 - len(b))))
    a = NP.concatenate((a, NP.zeros(N + 1 - len(a))))
    x = x[::-1]
    y = y[::-1]
    z[0] = NP.dot(b[1:], x) - NP.dot(a[1:], y)
    for i in range(1, N):
        z[i] = NP.dot(b[i+1:], x[:-i]) - NP.dot(a[i+1:], y[:-i])
    return z
