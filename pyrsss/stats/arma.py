import numpy as NP
import scipy as SP
import scipy.signal


def arma_predictor_model(x, y, m, n):
    """
    Return matrix and vector relating (*m*, *n*) ARMA model to the
    input *x* and output *y*. In other words, construct the max(*m*,
    *n*) - 1 by (*m* - 1 + *n* - 1) matrix A such that

    y[k] + a_1 y[k - 1] + ... + a_m y[k - m] = b_1 x[k - 1] + ... + b_n x[k - n]

    and the vector b corresponds to y[k] for k >= max(*m*, *n*).
    """
    assert len(x) == len(y)
    k = max(m, n)
    A1 = SP.linalg.toeplitz(-y[k:-1], r=-y[(k - m):k][::-1])
    A2 = SP.linalg.toeplitz(x[k:-1], r=x[(k - n):k][::-1])
    A = NP.hstack((A1, A2))
    b = y[k+1:]
    return A, b


def arma_predictor_linear(x, y, m, n):
    """
    Return the (*m*, *n*) one-step ARMA predictor trained on the input
    *x* and output *y*. The output is the tuple of the *m* AR and *n*
    MA coefficients.
    """
    assert len(x) == len(y)
    N = len(x)
    k = max(m, n)
    A, b = arma_predictor_model(x, y, m, n)
    x_hat = NP.linalg.lstsq(A, b)[0]
    a_hat = x_hat[:m]
    b_hat = x_hat[m:]
    a_hat = NP.insert(a_hat, 0, 1)
    return a_hat, b_hat


def residual(x_hat, m, x, y):
    """
    Return the residual between the ARMA filter *x_hat* (the tuple of
    *m* AR and MA coefficients) for the input *x* and output *y*.
    """
    n = len(x_hat) - m
    k = max(m, n)
    x = x[k:-1]
    y = y[k+1:]
    assert(len(x) == len(y))
    a_hat = x_hat[:m]
    b_hat = x_hat[m:]
    a_hat = NP.insert(a_hat, 0, 1)
    y_hat = scipy.signal.lfilter(b_hat, a_hat, x)
    return NP.abs(y_hat - y)


def Dfun(x_hat, m, x, y):
    """
    Return the Jacobian matrix for the ARMA model, i.e., the ith row
    corresponds to the ith data equation and the jth column is the jth
    parameter partial derivative.
    """
    n = len(x_hat) - m
    A, _ = arma_predictor_model(x, y, m, n)
    return -A


def arma_predictor_nonlinear(x, y, m, n, x_hat0=None):
    """
    Calculate the nonlinear fit between the (*m*, *n*) ARMA model and
    the input *x* and output *y*. The optimization starts at *x_hat*
    (a vector with all 0s when `None`). The output is the tuple of the
    *m* AR and *n* MA coefficients.
    """
    if x_hat0 is None:
        x_hat0 = NP.zeros(m + n)
    (x_hat,
     cov_x,
     info,
     mesg,
     ier) = SP.optimize.leastsq(residual,
                                x_hat0,
                                args=(m, x, y),
                                Dfun=Dfun,
                                full_output=True)
    if ier not in [1, 2, 3, 4]:
        raise RuntimeError('optimization failed (ier={}) --- {}'.fomat(ier,
                                                                       mesg))
    a_hat = x_hat[:m]
    b_hat = x_hat[m:]
    a_hat = NP.insert(a_hat, 0, 1)
    return a_hat, b_hat


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
