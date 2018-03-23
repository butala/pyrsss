from collections import namedtuple

import numpy as NP
import scipy as SP
import scipy.signal

from ..signal.lfilter import miso_lfilter


def arma_predictor_model(x, y, Na, Nb, Nk=1):
    """
    Return matrix and vector relating (*Na*, *Nb*) ARMA model to the
    input *x* and output *y*.

    In other words, construct the matrix with A_M = N - max(*Na*, *Nb*
    + *Nk* - 1) rows and *Na* + *Nb* columns such that

    y[n] + a_1 y[n - 1] + ... + a_Na y[n - Na] = b_1 x[n - k] + ... + b_Nb x[n - Nb - k + 1]

    for N - A_M <= n < N (where N is the length of *x* and *y*).
    """
    assert len(x) == len(y)
    N = len(x)

    A_M = N - max(Na, Nb + Nk - 1)
    if A_M <= 0:
        raise ValueError('problem specification results in matrix with {} rows'.format(A_M))

    A1r = -y[(N - A_M - Na):(N - A_M)][::-1]
    A1c = -y[(N - A_M - 1):-1]
    A1 = SP.linalg.toeplitz(A1c, r=A1r)

    A2r = x[(N - A_M - Nk - Nb + 1):(N - A_M - Nk + 1)][::-1]
    A2c = x[(N - A_M - Nk):(N - Nk)]
    A2 = SP.linalg.toeplitz(A2c, r=A2r)

    A = NP.hstack((A1, A2))
    b = y[(N - A_M):]
    return A, b


def arma_fit_linear(x, y, Na, Nb, Nk=1):
    """
    Return the (*Na*, *Nb*) ARMA predictor trained on the input
    *x* with a delay of *Nk* and output *y*. The output is the
    tuple of the *Na* AR and *Nb* MA coefficients.

    The implemented functionality compares with the single-input,
    single-output Matlab routine `arx` (see
    http://www.mathworks.com/help/ident/ref/arx.html).
    """
    assert len(x) == len(y)
    A, b = arma_predictor_model(x, y, Na, Nb, Nk=Nk)
    x_hat = NP.linalg.lstsq(A, b)[0]
    a_hat = x_hat[:Na]
    b_hat = x_hat[Na:]
    a_hat = NP.insert(a_hat, 0, 1)
    b_hat = NP.insert(x_hat[Na:], 0, NP.zeros(Nk))
    return a_hat, b_hat


def arma_residual(x_hat, Na, Nb, Nk, x, y):
    """
    Return the residual between the ARMA filter *x_hat* (the tuple of
    *Na* AR and *Nb* MA coefficients with input delay of *Nk* time
    steps) for the input *x* and output *y*.
    """
    assert(len(x) == len(y))
    assert len(x_hat) == Na + Nb
    a_hat = x_hat[:Na]
    b_hat = x_hat[Na:]
    a_hat = NP.insert(a_hat, 0, 1)
    b_hat = NP.insert(x_hat[Na:], 0, NP.zeros(Nk))
    y_hat = scipy.signal.lfilter(b_hat, a_hat, x)
    return y_hat - y


def arma_fit_nonlinear(x, y, Na, Nb, Nk=1, x_hat0=None, **kwds):
    """
    Calculate the nonlinear fit between the (*Na*, *Nb*, with *Nk*
    time step delay on the input) ARMA model and the input *x* and
    output *y*. The optimization starts at *x_hat0* (a vector with all
    0s when `None`). The output is the tuple of the AR and MA
    coefficients.

    The implemented functionality compares with the single-input,
    single-output Matlab routine `oe` (see
    http://www.mathworks.com/help/ident/ref/oe.html).
    """
    if x_hat0 is None:
        x_hat0 = NP.zeros(Na + Nb)
    (x_hat,
     cov_x,
     info,
     mesg,
     ier) = SP.optimize.leastsq(arma_residual,
                                x_hat0,
                                args=(Na, Nb, Nk, x, y),
                                full_output=True,
                                **kwds)
    if ier not in [1, 2, 3, 4]:
        raise RuntimeError('optimization failed (ier={}) --- {}'.format(ier,
                                                                        mesg))
    a_hat = x_hat[:Na]
    b_hat = x_hat[Na:]
    a_hat = NP.insert(a_hat, 0, 1)
    b_hat = NP.insert(x_hat[Na:], 0, NP.zeros(Nk))
    return a_hat, b_hat


def miso_pack(b, a):
    """
    Pack MISO transfer function components (numerator coefficients *b*
    and denominator coefficients *a*, both lists of vectors) into a
    vector representation.
    """
    assert len(b) == len(a)
    x = []
    for b_i, a_i in zip(b, a):
        if a_i[0] != 1:
            a_i /= a[0]
            b_i /= a[0]
        x.extend(b_i)
        x.extend(a_i[1:])
    return NP.array(x)


def miso_unpack(x, Na, Nb, Nk):
    """
    Unpack the MISO vector transfer function representation *x* to
    lists of numerator and denominator coefficients. *Na* is the list
    of the number of denominator terms per system, *Nb* is the list of
    the number of numerator terms per system, and *Nk* is the list of
    input delays per system.
    """
    assert len(x) == sum(Na) + sum(Nb)
    i = 0
    b = []
    a = []
    for Na_i, Nb_i, Nk_i in zip(Na, Nb, Nk):
        b_i = x[i:i+Nb_i]
        b.append(NP.insert(b_i, 0, NP.zeros(Nk_i)))
        i += Nb_i
        a_i = x[i:i+Na_i]
        a.append(NP.insert(a_i, 0, 1))
        i += Na_i
    return b, a


def miso_arma_residual(x_hat, Na, Nb, Nk, x, y, I=slice(None)):
    """
    Return the residual between the MISO ARMA filter *x_hat* (see
    func:`miso_pack`). the for the inputs *x* and output *y*. The
    lists *Na*, *Nb*, and *Nk* store the number of denominator,
    numerator, and input delay time steps per system. *I* indicates
    the slice of the residual to return (useful for leaving off the
    portion dominated by the transient response).
    """
    assert len(Na) == len(Nb) == len(Nk) == len(x)
    for x_i in x:
        assert len(x_i) == len(y)
    assert len(x_hat) == sum(Na) + sum(Nb)
    b_hat, a_hat = miso_unpack(x_hat, Na, Nb, Nk)
    y_hat = miso_lfilter(b_hat, a_hat, x)
    return (y_hat - y)[I]


def miso_arma_fit_nonlinear(x, y, Na, Nb, Nk=None, x_hat0=None, I=slice(None), **kwds):
    """
    Fit MISO ARMA model to list of inputs *x* and corresponding output
    *y*. Use a model with the list *Na* and *Nb* of the number of
    denominator and numerator parameters, respectively, per
    system. Use the list of *Nk* input delays per system (assume 1 for
    all systems if not provided). Start the optimization with the
    vector parameter set *x_hat0* (see :func:`pack` for the ordering
    and start from all 0s if not provided). *I* indicates the slice of
    the residuals to consider (see :func:`miso_arma_residual`).

    The implemented functionality compares with the multiple-input,
    single-output Matlab routine `oe` (see
    http://www.mathworks.com/help/ident/ref/oe.html).
    """
    assert len(x) == len(Na) == len(Nb)
    if Nk is None:
        Nk = [1] * len(Na)
    else:
        assert len(Nk) == len(x)
    if x_hat0 is None:
        x_hat0 = NP.zeros(sum(Na) + sum(Nb))
    N_set = set(map(len, x))
    assert len(N_set) == 1
    assert N_set.pop() == len(y)
    J = lambda *args: miso_arma_residual(*args, **{'I': I})
    (x_hat,
     cov_x,
     info,
     mesg,
     ier) = SP.optimize.leastsq(J,
                                x_hat0,
                                args=(Na, Nb, Nk, x, y),
                                full_output=True,
                                **kwds)
    if ier not in [1, 2, 3, 4]:
        raise RuntimeError('optimization failed (ier={}) --- {}'.format(ier,
                                                                        mesg))
    return miso_unpack(x_hat, Na, Nb, Nk)


"""ARMA fit assessment metrics""",
Assessment = namedtuple('Assessment',
                        'RSE '\
                        'R2 '\
                        'F_stat '\
                        'Cp '\
                        'AIC '\
                        'BIC '\
                        'adjusted_R2')


def assessment(x, y, m, n, x_hat):
    """
    WARNING --- UNTESTED!
    """
    # assert n + 1 == len(x_hat) - m
    # a_hat = x_hat[:m+1]
    # b_hat = x_hat[m+1:]
    a_hat, b_hat = x_hat
    m = len(a_hat)
    n = len(b_hat)
    k = max(m, n)
    y_hat = scipy.signal.lfilter(b_hat, a_hat, x[k:-1])
    N = len(y_hat)
    P = len(x_hat)
    RSS = NP.sum((y[k+1:] - y_hat)**2)
    RSE = NP.sqrt(RSS / (N - P - 1))
    TSS = NP.sum((y - NP.mean(y))**2)
    R2 = 1 - RSS / TSS
    F_stat = (TSS - RSS) / P / (RSS / (N - P - 1))
    sigma_hat = RSE
    d = P
    Cp = (RSS + 2 * d * sigma_hat**2) / N
    AIC = (RSS + 2 * d * sigma_hat**2) / (N * sigma_hat**2)
    BIC = (RSS + NP.log(N) * d * sigma_hat**2) / N
    adjusted_R2 = 1 - (RSS / (N - d - 1)) / (TSS / (N - 1))
    return Assessment(RSE,
                      R2,
                      F_stat,
                      Cp,
                      AIC,
                      BIC,
                      adjusted_R2)


if __name__ == '__main__':
    poles = [0.2, 0.6, 0.9]

    a_true = NP.poly(poles)
    b_true = [0, 0, 1, 0.5, 0, 0]

    m = len(a_true) - 1
    n = len(b_true) - 1
    k = max(m, n)

    N = 100
    x = NP.zeros(N)
    x[10] = 1

    y = scipy.signal.lfilter(b_true, a_true, x)

    a_hat, b_hat = arma_predictor_linear(x, y, m, n)
    a_hat_nl, b_hat_nl = arma_predictor_nonlinear(x, y, m, n)

    print(a_hat, a_hat_nl, a_true)
    print(b_hat, b_hat_nl, b_true)

    # print(a_hat)

    # print(m, n)
    # print(len(a_hat), len(b_hat))

    print(assessment(x, y, m, n, NP.concatenate((a_hat, b_hat))))
