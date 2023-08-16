import math

import numpy as np
import scipy as sp


def rect(t, a):
    """
    Return a vector of the same length as $t$ that is equal to 1 for
    absolute values of $t$ less than $a$/2 and 0 otherwise.
    """
    x = np.zeros_like(t)
    x[np.abs(t) < a/2] = 1
    x[np.abs(t) == a/2] = 1/2
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
        X = sp.fft.rfft(x, n=N) * T_s
        f = sp.fft.rfftfreq(N, d=T_s)
    else:
        X = sp.fft.fft(x, n=N) * T_s
        f = sp.fft.fftfreq(N, d=T_s)

        X = sp.fft.fftshift(X)
        f = sp.fft.fftshift(f)
    if n0 != 0:
        X *= np.exp(-1j * 2 * math.pi * np.arange(N) * n0 / N)
    return f, X


def blackman_tukey(x,
                   M,
                   L,
                   y=None,
                   window='boxcar',
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
    # For an explanation of the circular shift below:
    # https://dsp.stackexchange.com/questions/34642/practical-cross-spectrum-estimation-using-blackman-tukey-approach
    #
    # Also, consider adding support for a lab parameter. A
    # cross-correlation does not necessarily peak at k=0 and a lag,
    # k_0 in Stoica and Moses section 2.8.4, can be used to avoid
    # having the window attenuate the cross-correlation where it is
    # most significant.
    N = len(x)
    assert M <= N
    if y is None:
        y = x
        autocorrelation = True
    else:
        assert len(y) == N
        autocorrelation = False

    Rxy = sp.signal.correlate(x, y) / N
    Rxy_window = Rxy[(N - 1) - M:(N - 1) + M + 1]
    window = sp.signal.get_window(window, 2*M + 1, fftbins=False)
    #k_range = np.arange(0, L)
    #shift = np.exp(2j * np.pi * k_range * M / L)
    #Sxy = np.fft.fft(window * Rxy_window, n=L) * shift

    # I assume a circular shift (as below) would be more efficient
    # than multiplying by a complex exponential (as above) to achieve
    # the same result.
    zp = np.pad(window * Rxy_window, (0, L - (2*M+1)), constant_values=0)
    Sxy = sp.fft.fft(np.roll(zp, -M))
    if autocorrelation:
        Sxy = np.real(Sxy)
    f = sp.fft.fftfreq(L, d=d)
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


def etfe(x,
         y,
         M,
         L,
         d=1,
         window='parzen'):
    """
    Compute the empirical transfer function estimate (ETFE) relating
    the input time series *x* to the output time series *y*. Compute
    the response at *L* equally spaced frequency samples (where the
    sampling period is *d*). Limit the correlations to a lag of *M*
    (and *M* <= len(*x*) - 1) and use the window function *window*
    (see :func:`scipy.signal.get_window`). Return the tuple containing
    the ETFE and the frequency sample points.

    See Section 6.3 of Ljung, System Identification Theory for the
    User, 2nd Edition.
    """
    Phi_yu, f = blackman_tukey(y, M, L, y=x, d=d, window=window)
    Phi_u, _ = blackman_tukey(x, M, L, d=d, window=window)
    return Phi_yu / Phi_u, f


def etfe_welch(x,
               y,
               fs=1,
               **kwds):
    """
    """
    assert len(x) == len(y)
    fxx, Pxx = sp.signal.welch(x, fs=fs, **kwds)
    fxy, Pxy = sp.signal.csd(x, y, fs=fs, **kwds)
    return fxx, Pxy / Pxx


def bartlett(x, fs=1.0, window='boxcar', nfft=None):
    """Estimate the power spectrum of the sequence *x* using Bartlett's
    method (the method of averaged periodograms). The sequence is
    split into consecutive, disjoint, length *nfft* segments. Each
    segment is windowed and real FFT'ed. The periodogram is estimated
    as the mean of the squared magnitude over the set of processed
    segments. See `scipy.signal.get_window` for available options for
    *window*. Return a tuple with the vector of frequency sample
    points (taking into consideration the sampling frequency *fs* and
    the power spectral density estimate.
    """
    if nfft is None:
        nfft = len(x)
    nseg = len(x) // nfft
    Pxx = np.zeros(nfft // 2 + 1)
    window = sp.signal.get_window(window, nfft)
    while True:
        x_i, x = x[:nfft], x[nfft:]
        if len(x_i) < nfft:
            break
        Pxx += np.abs(sp.fft.rfft(x_i * window, n=nfft))**2 / nseg / nfft
    return sp.fft.rfftfreq(nfft, d=1/fs), Pxx


def bartlett_csd(y, x, fs=1.0, window='boxcar', nfft=None):
    """Estimate the cross power spectrum of the sequences *x* and *y*
    using Bartlett's method (the method of averaged periodograms). The
    sequence is split into consecutive, disjoint, length *nfft*
    segments. Each segment is windowed and real FFT'ed. The
    periodogram is estimated as the mean of the squared magnitude over
    the set of processed segments. See `scipy.signal.get_window` for
    available options for *window*. Return a tuple with the vector of
    frequency sample points (taking into consideration the sampling
    frequency *fs* and the power spectral density estimate.
    """
    # NOTE: this code could be used to generalize the bartlett
    # function above (with x=y). However, the code would then
    # needlessly calculate rfft twice (once for x and y when
    # x=y). Instead of writing complicating logic to avoid the
    # needless calculation, I have left the two functions separate.
    # --- MDB
    assert len(y) == len(x)
    if nfft is None:
        nfft = len(x)
    nseg = len(x) // nfft
    Pyx = np.zeros(nfft // 2 + 1, dtype=complex)
    window = sp.signal.get_window(window, nfft)
    while True:
        y_i, y = y[:nfft], y[nfft:]
        x_i, x = x[:nfft], x[nfft:]
        if len(x_i) < nfft:
            break
        Pyx += sp.fft.rfft(y_i * window, n=nfft) * np.conj(sp.fft.rfft(x_i * window, n=nfft)) / nseg / nfft
    return sp.fft.rfftfreq(nfft, d=1/fs), Pyx


def etfe_bartlett(x,
                  y,
                  fs=1,
                  **kwds):
    """
    """
    assert len(x) == len(y)
    fxx, Pxx = bartlett(x, fs=fs, **kwds)
    fyx, Pyx = bartlett_csd(y, x, fs=fs, **kwds)
    return fxx, np.conj(Pyx) / Pxx


if __name__ == '__main__':
    import matplotlib.pyplot as plt

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
        return np.exp(-1j * w * (L - 1) / 2) * np.sin(w * L / 2) / np.sin(w / 2)

    K = 2048    # number of frequency samples
    w = np.linspace(-math.pi, math.pi, K)
    W = W_r(w, N)

    n = np.arange(N)
    v = np.cos(2*math.pi/14 * n) + 0.75 * np.cos(4*math.pi/15 * n)

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

    f2 = np.linspace(-fs/2, fs/2, K)
    w2 = f2 / fs * 2 * math.pi
    V_dtft = V(w2, A0, 2*math.pi/14, A1, 4*math.pi/15, N)

    f_spectrum, V_spectrum = spectrum(v,
                                      T_s=T,
                                      only_positive=False)

    f_spectrum2, V_spectrum2 = spectrum(v,
                                        T_s=T,
                                        oversample=2,
                                        only_positive=False)

    f_spectrum3, V_spectrum3 = spectrum(v,
                                        T_s=T)

    plt.figure(figsize=(10, 4))
    plt.plot(f2,
            np.abs(V_dtft) * T,
            c='C0',
            label='DTFT (scaled)')
    plt.scatter(f_spectrum,
               np.abs(V_spectrum),
               c='C1',
               s=10,
               label='spectrum')
    plt.scatter(f_spectrum2,
               np.abs(V_spectrum2),
               c='C2',
               zorder=-1,
               s=20,
               label='spectrum (oversample=2)')
    plt.scatter(f_spectrum3,
               np.abs(V_spectrum3),
               c='C3',
               zorder=-2,
               s=40,
               label='spectrum (positive-only)')
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.title('Comparison of pyrsss spectrum and scaled DTFT')

    plt.show()
