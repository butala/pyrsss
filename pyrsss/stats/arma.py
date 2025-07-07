from collections import namedtuple

import numpy as np
import scipy as sp
import scipy.signal

from ..signal.lfilter import miso_lfilter


def arma_predictor_model(x, y, Na, Nb, Nk=1):
    """
    Return matrix and vector relating (*Na*, *Nb*) ARMA model to the
    input *x* and output *y*.

    In other words, construct the matrix with A_M = N - max(*Na*, *Nb*
    + *Nk* - 1) rows and *Na* + *Nb* columns such that

    y[n] + a_1 y[n - 1] + ... + a_Na y[n - Na] = b_1 x[n - Nk] + ... + b_Nb x[n - Nb - Nk + 1]

    for N - A_M <= n < N (where N is the length of *x* and *y*).
    """
    assert len(x) == len(y)
    N = len(x)

    A_M = N - max(Na, Nb + Nk - 1)
    if A_M <= 0:
        raise ValueError('problem specification results in matrix with {} rows'.format(A_M))

    A1r = -y[(N - A_M - Na):(N - A_M)][::-1]
    A1c = -y[(N - A_M - 1):-1]
    A1 = sp.linalg.toeplitz(A1c, r=A1r)

    A2r = x[(N - A_M - Nk - Nb + 1):(N - A_M - Nk + 1)][::-1]
    A2c = x[(N - A_M - Nk):(N - Nk)]
    A2 = sp.linalg.toeplitz(A2c, r=A2r)

    A = np.hstack((A1, A2))
    b = y[(N - A_M):]
    return A, b


def get_a_b(theta, Na, Nb, Nk=1):
    """
    Convert the vector of ARMA model parameters *theta* =
    [theta_a, theta_b] to the tuple (theta_a, theta_b) of polynomial
    transfer function parameters (a, b) where a = (1, theta_a) are the
    denominator and b = (zeros(Nk), theta_b) are the numerator
    coefficients, respectively.
    """
    a = theta[:Na]
    b = theta[Na:]
    a = np.insert(a, 0, 1)
    b = np.insert(b, 0, np.zeros(Nk))
    return a, b


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
    x_hat = np.linalg.lstsq(A, b)[0]
    a_hat, b_hat = get_a_b(x_hat, Na, Nb, Nk=Nk)
    return a_hat, b_hat


def arma_residual(x_hat, Na, Nb, Nk, x, y):
    """
    Return the residual between the ARMA filter *x_hat* (the tuple of
    *Na* AR and *Nb* MA coefficients with input delay of *Nk* time
    steps) for the input *x* and output *y*.
    """
    assert(len(x) == len(y))
    assert len(x_hat) == Na + Nb
    a_hat, b_hat = get_a_b(x_hat, Na, Nb, Nk=Nk)
    y_hat = scipy.signal.lfilter(b_hat, a_hat, x)
    return y_hat - y


def arma_sensitivity(b, a, x, Nk, zi=None):
    """
    Return the matrix \\partial y(n_i, theta) / \\partial theta_j
    where theta = [a, b] and

    y(n_i, theta) = lfilter(b, a, x[:n_i], zi=zi).

    where n_i < x.shape[0] (number of ARMA inputs/outputs) and theta_j
    < len(a) + len(b) (the number of ARMA filter parameters).
    """
    assert np.isclose(a[0], 1)
    assert Nk < len(b)
    if Nk > 0:
        assert np.isclose(b[:Nk], 0).all()

    # SWITCH TO SOS FILTER!
    w = sp.signal.lfilter(1, a, x)
    z = sp.signal.lfilter(-1, a, w)
    v = sp.signal.lfilter(b, 1, z)

    if zi is not None:
        w2 = sp.signal.lfilter(1, a, np.pad(zi, (0, len(x)-len(zi))))
        v += sp.signal.lfilter(-1, a, w2)

    Da_r = np.zeros(len(a) - 1)
    Da_c = np.pad(v[:-1], (1, 0))
    Da = sp.linalg.toeplitz(Da_c, r=Da_r)

    Db_r = np.zeros(len(b) - Nk)
    Db_c = np.pad(w[:(len(w) - Nk)], (Nk, 0))
    if Nk == 0:
        Db_r[0] = Db_c[0]
    Db = sp.linalg.toeplitz(Db_c, r=Db_r)

    return np.hstack([Da, Db])


def arma_zi_sensitivity(Nb, a, N, Nk=0):
    """
    Return the matrix \\partial y(n_i, theta) / \\partial v_j
    where theta = [a, b], 1 \\leq i \\leq N, and v[j] for 1 \\leq j
    \\leq min(Nb, len(a)-1) are the initial filter states zi for a
    type 2 transposed ARMA filter as is implemented by, e.g.,
    `:func:scipy.signal.lfilter`.
    """
    assert N >= 1
    # SWITCH TO SOS FILTER!
    c = sp.signal.lfilter(1, a, np.pad([1], (0, N-1)))
    r = np.zeros(max(Nb+Nk, len(a))-1)
    return sp.linalg.toeplitz(c, r=r)


def arma_l2_norm_sensitivity(b, a, x, y_target, Nk):
    """
    Return the gradient of the L2 norm term function

    || lfilter(b, a, x) - y_target ||_2^2

    as a len(a)-1 + len(b)-Nk vector. Note that a[0] = 1 and b[:Nk] = 0.
    """
    y = sp.signal.lfilter(b, a, x)
    return 2 * arma_sensitivity(b, a, x) @ (y - y_target)


def arma_jacobian(x_hat, Na, Nb, Nk, x, y):
    """
    Return the Jacobian of the function

    lfilter(b, a, x)

    as a len(a)-1 + len(b)-Nk vector = len(x_hat). See :func:`get_a_b`
    for the mapping from *x_hat* to the polynomial transfer function
    parameters and b. Note that a[0] = 1 and b[:Nk] = 0.
    """
    a_hat, b_hat = get_a_b(x_hat, Na, Nb, Nk=Nk)
    return arma_sensitivity(b_hat, a_hat, x, Nk)


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
        x_hat0 = np.zeros(Na + Nb)

    result = sp.optimize.least_squares(arma_residual,
                                       x_hat0,
                                       jac=arma_jacobian,
                                       args=(Na, Nb, Nk, x, y),
                                       **kwds)

    if not result.success:
        raise RuntimeError(result.message)
    x_hat = result.x
    return get_a_b(x_hat, Na, Nb, Nk=Nk)


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
    return np.array(x)


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
        b.append(np.insert(b_i, 0, np.zeros(Nk_i)))
        i += Nb_i
        a_i = x[i:i+Na_i]
        a.append(np.insert(a_i, 0, 1))
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
        x_hat0 = np.zeros(sum(Na) + sum(Nb))
    N_set = set(map(len, x))
    assert len(N_set) == 1
    assert N_set.pop() == len(y)
    J = lambda *args: miso_arma_residual(*args, **{'I': I})
    (x_hat,
     cov_x,
     info,
     mesg,
     ier) = sp.optimize.leastsq(J,
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
    RSS = np.sum((y[k+1:] - y_hat)**2)
    RSE = np.sqrt(RSS / (N - P - 1))
    TSS = np.sum((y - np.mean(y))**2)
    R2 = 1 - RSS / TSS
    F_stat = (TSS - RSS) / P / (RSS / (N - P - 1))
    sigma_hat = RSE
    d = P
    Cp = (RSS + 2 * d * sigma_hat**2) / N
    AIC = (RSS + 2 * d * sigma_hat**2) / (N * sigma_hat**2)
    BIC = (RSS + np.log(N) * d * sigma_hat**2) / N
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

    a_true = np.poly(poles)
    b_true = [0, 0, 1, 0.5, 0, 0]

    m = len(a_true) - 1
    n = len(b_true) - 1
    k = max(m, n)

    N = 10000
    #x = np.zeros(N)
    #x[10] = 1
    x = np.random.randn(N)

    y = scipy.signal.lfilter(b_true, a_true, x)

    a_hat, b_hat = arma_fit_linear(x, y, m, n)
    a_hat_nl, b_hat_nl = arma_fit_nonlinear(x, y, m, n)

    print(f'truth a:     {a_true}')
    print(f'linear a:    {a_hat}')
    print(f'nonlinear a: {a_hat_nl}')

    print()

    print(f'truth b:     {b_true}')
    print(f'linear b:    {b_hat}')
    print(f'nonlinear b: {b_hat_nl}')

    # print(assessment(x, y, m, n, np.concatenate((a_hat, b_hat))))
