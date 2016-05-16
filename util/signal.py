from __future__ import division

import math

import numpy as NP


def rect(t, a):
    """
    ???
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
             oversample=1):
    """
    ???
    """
    N = nextpow2(len(x)) * oversample
    X = NP.fft.fft(x, n=N) * T_s
    f = NP.fft.fftfreq(N, d=T_s)
    if n0 != 0:
        X *= NP.exp(-1j * 2 * math.pi * NP.arange(N) * (n0 + (N - len(x))) / N)
        #print(NP.exp(-1j * 2 * math.pi * NP.arange(N) * (n0 + (N - len(x))) / N))
        # HACK
        X[f < 0] *= -1
    return (NP.fft.fftshift(X),
            NP.fft.fftshift(f))


# def convolve(x, y, T_s=1):
#     """
#     ???
#     """
#     N_x = len(x)
#     N_y = len(y)
#     # compute the power next largest power of 2
#     N2_x = int(math.ceil(math.log(N_x, 2)))
#     N2_y = int(math.ceil(math.log(N_y, 2)))
#     N2 = max(N2_x, N2_y)
#     # compute FFTs
#     X = NP.fft.fft(x, n=N2)
#     Y = NP.fft.fft(y, n=N2)



if __name__ == '__main__':
    import pylab as PL

    # Test 1: Plot spectrum of sinusoid
    phi = 0  # phase offset [radians]
    f0  = 5  # frequency [Hz]
    A   = 1  # amplitude

    T_max = 1   # time max (time range from -T_max to T_max) [s]
    T_s = 1e-2  # sampling frequency [s]

    t = NP.arange(-T_max,
                  T_max + T_s / 2,
                  T_s)

    assert len(t) % 2 == 1 # even case not implemented
    n0 = (len(t) - 1) / 2

    x = NP.sin(2 * math.pi * f0 * t + phi)

    X, f = spectrum(x, n0=n0, T_s=T_s, oversample=4)

    fig = PL.figure()
    PL.subplot(211)
    PL.scatter(t, x, edgecolors='None')
    PL.xlim(-T_max, T_max)
    PL.xlabel('Time [s]')
    PL.subplot(212)
    PL.plot(f, NP.abs(X), label='DFT')
    PL.axvline(f0, c='r', label='FT')
    PL.axvline(-f0, c='r')
    PL.xlim(-2 * f0, 2 * f0)
    PL.legend()
    PL.xlabel('Frequency [Hz]')
    PL.suptitle('Sinusoid Example')

    # Test 3: Plot spectrum of rect
    a = 1  # width of rect function

    # T_max = 1   # time max (time range from -T_max to T_max) [s]
    T_max = 10   # time max (time range from -T_max to T_max) [s]
    #T_s = 1e-2  # sampling frequency (instead, calculate given a)
    T_s = 1e-2  # sampling frequency (instead, calculate given a)

    t = NP.arange(-T_max,
                  T_max + T_s / 2,
                  T_s)
    assert len(t) % 2 == 1 # even case not implemented
    n0 = (len(t) - 1) / 2
    # n0 = len(t) / 2
    # BUT: n0 = len(t) / 2 results in 0 imaginary part

    x = rect(t, a)

    X, f = spectrum(x, n0=n0, T_s=T_s, oversample=4)

    print(NP.linalg.norm(NP.imag(X)))

    N_FT = 4 * len(f)
    f_FT = NP.linspace(f[0], f[-1], N_FT)
    X_FT = abs(a) * NP.sinc(f_FT * a)

    fig = PL.figure()
    PL.subplot(211)
    PL.stem(t, x)
    PL.xlim(-T_max, T_max)
    PL.xlabel('Time [s]')
    PL.subplot(212)
    PL.plot(f, NP.real(X), label='DFT')
    PL.plot(f_FT, X_FT, label='FT', color='r')

    PL.suptitle('Rect Example')

    PL.show()
