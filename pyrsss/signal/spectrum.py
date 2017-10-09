from __future__ import division

import math

import numpy as NP
import scipy.signal


def rect(t, a):
    """
    Return a vector of the same length as $t$ that is equal to 1 for
    absolute values of $t$ less than $a$/2 and 0 otherwise.
    """
    x = NP.zeros_like(t)
    x[NP.abs(t) < a/2] = 1
    x[NP.abs(t) == a/2] = 1/2
    return x


def nextpow2(N):
    """
    Return the power of 2 greater than or equal to *N*.
    """
    return 2**int(math.ceil(math.log(N, 2)))


def spectrum(x,
             n0=0,
             T_s=1,
             oversample=1,
             only_positive=True):
    """
    Return the spectrum for the signal *x* calculated via FFT and the
    associated frequencies as a tuple. The *n0* parameter gives the
    index in *x* for time index 0 (*n0* = 0 means that `x[0]` is at
    time 0). The number of spectral samples returned is the next power
    of 2 greater than the length of *x* multiplied by *oversample*. If
    *only_positive*, return the spectrum only for positive frequencies
    (assuming *x* is real).
    """
    assert oversample >= 1 and isinstance(oversample, int)
    N = nextpow2(len(x)) * 2**(oversample - 1)
    if only_positive:
        X = NP.fft.rfft(x, n=N) * T_s
        f = NP.fft.rfftfreq(N, d=T_s)
    else:
        X = NP.fft.fft(x, n=N) * T_s
        f = NP.fft.fftfreq(N, d=T_s)

        X = NP.fft.fftshift(X)
        f = NP.fft.fftshift(f)
    if n0 != 0:
        X *= NP.exp(-1j * 2 * math.pi * NP.arange(N) * n0 / N)
    return f, X


def blackman_tukey(x,
                   M,
                   L,
                   y=None,
                   window='boxcar',
                   window_args=[],
                   d=1,
                   full=False):
    """
    Compute the Blackman-Tukey cross power spectral density (PSD)
    estimate between the time-domain signals *x* and *y* (must be the
    same length as *x*). If *y* is not given, compute the power
    spectral density estimate of *x*.  Use the spectral window with
    identifier *window* (see the options in
    :func:scipy.`signal.get_window`, e.g., a tuple can be used to pass
    arguments to the window function) and length *M* (i.e., the
    maximum auto-correlation lag to include in the estimate). Compute
    the estimate at *L* uniformly spaced frequency samples where *d*
    is the time domain sample interval. If not *full*, return the
    tuple containing the length *L* PSD estimate and length *L*
    corresponding frequencies. If *full*, also return the estimated
    cross correlation and window function (i.e., a tuple with four
    elements).
    """
    N = len(x)
    assert M <= N
    if y is None:
        y = x
    else:
        assert len(y) == N
    Rxy = scipy.signal.correlate(x, y) / N
    Rxy_window = Rxy[(N - 1) - M:(N - 1) + M + 1]
    window = scipy.signal.get_window(window, 2*M + 1, fftbins=False)
    k_range = NP.arange(0, L)
    shift = NP.exp(2j * NP.pi * k_range * M / L)
    Sxy = NP.fft.fft(window * Rxy_window, n=L) * shift * d
    f = NP.fft.fftfreq(L, d=d)
    if full:
        return (Sxy, f, Rxy, window)
    else:
        return (Sxy, f)


def periodogram(x,
                L,
                y=None,
                d=1,
                full=False):
    """
    Compute the periodogram of the cross power spectral density of *x*
    and *y*. The implementation is based on :func:`blackman-tukey`,
    following the same input and output conventions.
    """
    return blackman_tukey(x, len(x) - 1, L, y=y, d=d, full=full)


if __name__ == '__main__':
    import pylab as PL

    # Reproduction of Oppenheim and Schafer (O&S), 3rd edition,
    # Example 10.4. The example considers the effects of windowing and
    # frequency sampling of a two sinusoid superposition example.

    N = 64  # number of samples (window length)

    fs = 10e3   # sampling frequency (Hz)
    T = 1 / fs  # sample period (s)

    def W_r(w, L):
        """
        DTFT of rectangular window of length *L* evaluated an angular
        frequencies *w*. See O&S (10.11).
        """
        return NP.exp(-1j * w * (L - 1) / 2) * NP.sin(w * L / 2) / NP.sin(w / 2)

    K = 2048    # number of frequency samples
    w = NP.linspace(-math.pi, math.pi, K)
    W = W_r(w, N)

    n = NP.arange(N)
    v = NP.cos(2*math.pi/14 * n) + 0.75 * NP.cos(4*math.pi/15 * n)

    def V(w, A0, w0, A1, w1, L):
        """
        DTFT of the superposition of two sinusoids with amplitudes *A0*
        and *A1*, angular frequencies *w0* and *w1*, and a recangular
        window length *L* and angular frequencies *w*.
        """
        V1 = A0 / 2 * W_r(w - w0, L)
        V2 = A0 / 2 * W_r(w + w0, L)
        V3 = A1 / 2 * W_r(w - w1, L)
        V4 = A1 / 2 * W_r(w + w1, L)
        return V1 + V2 + V3 + V4

    A0 = 1
    A1 = 0.75

    f2 = NP.linspace(-fs/2, fs/2, K)
    w2 = f2 / fs * 2 * math.pi
    V_dtft = V(w2, A0, 2*math.pi/14, A1, 4*math.pi/15, N)

    V_spectrum, f_spectrum = spectrum(v,
                                      T_s=T,
                                      only_positive=False)

    V_spectrum2, f_spectrum2 = spectrum(v,
                                        T_s=T,
                                        oversample=2,
                                        only_positive=False)

    V_spectrum3, f_spectrum3 = spectrum(v,
                                        T_s=T)

    PL.figure(figsize=(10, 4))
    PL.plot(f2,
            NP.abs(V_dtft) * T,
            c='C0',
            label='DTFT (scaled)')
    PL.scatter(f_spectrum,
               NP.abs(V_spectrum),
               c='C1',
               s=10,
               label='spectrum')
    PL.scatter(f_spectrum2,
               NP.abs(V_spectrum2),
               c='C2',
               zorder=-1,
               s=20,
               label='spectrum (oversample=2)')
    PL.scatter(f_spectrum3,
               NP.abs(V_spectrum3),
               c='C3',
               zorder=-2,
               s=40,
               label='spectrum (positive-only)')
    PL.legend()
    PL.xlabel('Frequency (Hz)')
    PL.ylabel('Amplitude')
    PL.title('Comparison of pyrsss spectrum and scaled DTFT')

    PL.show()
