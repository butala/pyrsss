from __future__ import division

from itertools import izip, repeat
from collections import namedtuple

import numpy as NP
import scipy as SP
from scipy.linalg import cho_factor, cho_solve


"""
Class to store Kalman filter results.
"""
class FilterResult(namedtuple('FilterResult', 'x_hat P')):
    pass


def kalman_filter(y, H, R, F, Q, mu, PI, z=None):
    """
    Given the following sequences (one item of given dimension for
    each time step):
    - *y*: measurements (M)
    - *H*: measurement operator (MxN)
    - *R*: measurement noise covariance (MxM)
    - *F*: time update operator (NxN)
    - *Q*: time update noise covariance (NxN)
    - *mu*: initial state (N)
    - *PI*: initial state covariance (NxN)
    - *z*: (optional) systematic time update input (N)

    Return the :class:`FilterResult` containing lists of posterior
    state estimates and error covariances.
    """
    x_hat = []
    P = []
    x_hat_prior = mu
    P_prior = PI
    if z is None:
        z = repeat(None)
    for i, (y_i, H_i, R_i, F_i, Q_i, z_i) in enumerate(izip(y, H, R, F, Q, z)):
        # measurement update
        A = cho_factor(NP.matmul(H_i, NP.matmul(P_prior, H_i.T)) + R_i)
        B = cho_solve(A, NP.matmul(H_i, P_prior))
        b = cho_solve(A, y_i - NP.matmul(H_i, x_hat_prior))
        C = NP.matmul(P_prior, H_i.T)
        x_hat.append(x_hat_prior + NP.matmul(C, b))
        P.append(P_prior - NP.matmul(C, B))
        # time update
        x_hat_prior = NP.matmul(F_i, x_hat[-1])
        if z_i is not None:
            x_hat_prior += z_i
        P_prior = NP.matmul(F_i, NP.matmul(P[-1], F_i.T)) + Q_i
    return FilterResult(x_hat, P)


"""
Class to store square root  Kalman filter results.
"""
class SqrtFilterResult(namedtuple('SqrtFilterResult', 'x_hat P_sqrt')):
    pass


def sqrt_kf_mu(x_hat_prior,
               P_sqrt_prior,
               y_i,
               H_i,
               R_sqrt_i):
    """
    Square root Kalman filter measurement update: Given the following:
    - *x_hat_prior*: prior state estimate (N)
    - *P_sqrt_prior: prior error covariance square root (NxN)
    - *y_i*: measurements (M)
    - *H_i*: measurement operator (MxN)
    - *R_sqrt_i*: measurement noise covariance square root (MxM)

    Return the tuple containing the posterior state estimate and
    posterior error covariance square root.
    """
    M, N = H_i.shape
    A_T = NP.block([[R_sqrt_i,         NP.matmul(H_i, P_sqrt_prior)],
                    [NP.zeros((N, M)), P_sqrt_prior]])
    R_T = NP.linalg.qr(A_T.T, mode='r')
    P_sqrt_posterior = R_T.T[-N:, -N:]
    R_e_i_sqrt = R_T.T[:M,:M]
    b = y_i - NP.matmul(H_i, x_hat_prior)
    b2 = SP.linalg.solve_triangular(R_e_i_sqrt, b, lower=True, overwrite_b=True)
    A2 = R_T.T[-N:, :M]
    x_hat_posterior = x_hat_prior + NP.matmul(A2, b2)
    return x_hat_posterior, P_sqrt_posterior


def sqrt_kf_tu(x_hat_posterior,
               P_sqrt_posterior,
               F_i,
               Q_sqrt_i,
               z_i=None):
    """
    Square root Kalman filter time update. Given the following:
    - *x_hat_posterior*: posterior state estimate (N)
    - *P_sqrt_posterior*: posterior error covariance square root (NxN)
    - *F_i*: time update operator (NxN)
    - *Q_sqrt_i*: time update noise covariance square root (NxN)
    - *z_i*: optional) systematic time update input (N)

    Return the tuple containing the one time step prediction of the
    state and the square root of the error covariance.
    """
    N, _ = F_i.shape
    x_hat_prior = NP.matmul(F_i, x_hat_posterior)
    if z_i is not None:
        x_hat_prior += z_i
    A_T = NP.block([NP.matmul(F_i, P_sqrt_posterior), Q_sqrt_i])
    R_T = NP.linalg.qr(A_T.T, mode='r')
    P_sqrt_prior = R_T.T[:, :N]
    return x_hat_prior, P_sqrt_prior


def sqrt_kalman_filter(y,
                       H,
                       R_sqrt,
                       F,
                       Q_sqrt,
                       mu,
                       PI_sqrt,
                       z=None,
                       callback=None):
    """
    Given the following sequences (one item of given dimension for
    each time step):
    - *y*: measurements (M)
    - *H*: measurement operator (MxN)
    - *R_sqrt*: measurement noise covariance square root (MxM)
    - *F*: time update operator (NxN)
    - *Q_sqrt*: time update noise covariance square root (NxN)
    - *mu*: initial state (N)
    - *PI_sqrt*: initial state covariance square root (NxN)
    - *z*: (optional) systematic time update input (N)

    Return the :class:`SqrtFilterResult` containing lists of posterior
    state estimates and error covariance square roots.

    With infinite numerical precision, the state estimates computed
    here will be identical to :func:`kalman_filter`. With finite
    precision, the posterior estimates computed by this routine are
    more robust: covariances are operated upon in square root form
    (ensuring covariances are always positive definite) and numerical
    manipulations are through unitary operations.

    The function *callback*, if provided, is called at the end of each
    measurement / time update cycle. One argument is passed: the time
    index (starting from 0).

    Reference: Kailath, Sayed, and Hassibi, Linear Estimation, Chapter
    12
    """
    x_hat = []
    P_sqrt = []
    x_hat_prior = mu
    P_sqrt_prior = PI_sqrt
    if z is None:
        z = repeat(None)
    for i, (y_i, H_i, R_sqrt_i, F_i, Q_sqrt_i, z_i) in enumerate(izip(y,
                                                                      H,
                                                                      R_sqrt,
                                                                      F,
                                                                      Q_sqrt,
                                                                      z)):
        # measurement update
        (x_hat_posterior,
         P_sqrt_posterior) = sqrt_kf_mu(x_hat_prior,
                                        P_sqrt_prior,
                                        y_i,
                                        H_i,
                                        R_sqrt_i)
        x_hat.append(x_hat_posterior)
        P_sqrt.append(P_sqrt_posterior)
        # time update
        (x_hat_prior,
         P_sqrt_prior) = sqrt_kf_tu(x_hat_posterior,
                                    P_sqrt_posterior,
                                    F_i,
                                    Q_sqrt_i,
                                    z_i=z_i)
        if callback:
            callback(i)
    return SqrtFilterResult(x_hat, P_sqrt)
