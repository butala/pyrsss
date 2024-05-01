from itertools import repeat
from dataclasses import dataclass
import math

import numpy as np
import numpy.typing as npt
from scipy.linalg import solve_triangular, cho_factor, cho_solve


"""
Class to store Kalman filter results.
"""
try:
    @dataclass
    class FilterResult:
        x_hat: list[npt.ArrayLike]
        P: list[npt.ArrayLike]
except Exception:
    from collections import namedtuple

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
    for i, (y_i, H_i, R_i, F_i, Q_i, z_i) in enumerate(zip(y, H, R, F, Q, z)):
        # measurement update
        A = cho_factor(H_i @ (P_prior @ H_i.T) + R_i)
        B = cho_solve(A, H_i @ P_prior)
        b = cho_solve(A, y_i - H_i @ x_hat_prior)
        C = P_prior @ H_i.T
        x_hat.append(x_hat_prior + C @ b)
        P.append(P_prior - C @ B)
        # time update
        x_hat_prior = F_i @ x_hat[-1]
        if z_i is not None:
            x_hat_prior += z_i
        P_prior = F_i @ P[-1] @ F_i.T + Q_i
    return FilterResult(x_hat, P)


def extended_kalman_filter(y, H_fun, H_jac, R, F_fun, F_jac, Q, mu, PI, z=None):
    """
    Given the following sequences (one item of given dimension for
    each time step):
    - *y*: measurements (M)
    - *H_fun*: function that calculates nonlinear observation (N) -> (M)
    - *H_jac*: function that calculates the observation Jacobian (N) -> (MxN)
    - *R*: measurement noise covariance (MxM)
    - *F_fun*: function that calculates the nonlinear state update (N) -> (N)
    - *F_jac*: function that calculates the state update Jacobian (N) -> (NxN)
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
    for i, (y_i, R_i, Q_i, z_i) in enumerate(zip(y, R, Q, z)):
        # measurement update
        H_i = H_jac(x_hat_prior)
        A = cho_factor(H_i @ (P_prior @ H_i.T) + R_i)
        B = cho_solve(A, H_i @ P_prior)
        b = cho_solve(A, y_i - H_fun(x_hat_prior))
        C = P_prior @ H_i.T
        x_hat.append(x_hat_prior + C @ b)
        P.append(P_prior - C @ B)
        # time update
        F_i = F_jac(x_hat[-1])
        x_hat_prior = F_fun(x_hat[-1])
        if z_i is not None:
            x_hat_prior += z_i
        P_prior = F_i @ P[-1] @ F_i.T + Q_i
    return FilterResult(x_hat, P)


"""
Class to store square root Kalman filter results.
"""
try:
    @dataclass
    class SqrtFilterResult:
        x_hat: list[npt.ArrayLike]
        P_sqrt: list[npt.ArrayLike]
except Exception:
    from collections import namedtuple

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
    A_T = np.block([[R_sqrt_i,         H_i @ P_sqrt_prior],
                    [np.zeros((N, M)), P_sqrt_prior]])
    R_T = np.linalg.qr(A_T.T, mode='r')
    P_sqrt_posterior = R_T.T[-N:, -N:]
    R_e_i_sqrt = R_T.T[:M,:M]
    b = y_i - H_i @ x_hat_prior
    b2 = solve_triangular(R_e_i_sqrt, b, lower=True, overwrite_b=True)
    A2 = R_T.T[-N:, :M]
    x_hat_posterior = x_hat_prior + A2 @ b2
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
    x_hat_prior = F_i @ x_hat_posterior
    if z_i is not None:
        x_hat_prior += z_i
    A_T = np.block([F_i @ P_sqrt_posterior, Q_sqrt_i])
    R_T = np.linalg.qr(A_T.T, mode='r')
    P_sqrt_prior = R_T.T[:, :N]
    return x_hat_prior, P_sqrt_prior


def sqrt_kalman_filter_yield(y,
                             H,
                             R_sqrt,
                             F,
                             Q_sqrt,
                             mu,
                             PI_sqrt,
                             z=None,
                             callback=None):
    """Given the following sequences (one item of given dimension for
    each time step):
    - *y*: measurements (M)
    - *H*: measurement operator (MxN)
    - *R_sqrt*: measurement noise covariance square root (MxM)
    - *F*: time update operator (NxN)
    - *Q_sqrt*: time update noise covariance square root (NxN)
    - *mu*: initial state (N)
    - *PI_sqrt*: initial state covariance square root (NxN)
    - *z*: (optional) systematic time update input (N)

    Yields the tuple (x_hat_posterior, P_sqrt_posterior) after each
    measurement update.

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
    x_hat_prior = mu
    P_sqrt_prior = PI_sqrt
    if z is None:
        z = repeat(None)
    for i, (y_i, H_i, R_sqrt_i, F_i, Q_sqrt_i, z_i) in enumerate(zip(y,
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
        yield x_hat_posterior, P_sqrt_posterior
        # time update
        (x_hat_prior,
         P_sqrt_prior) = sqrt_kf_tu(x_hat_posterior,
                                    P_sqrt_posterior,
                                    F_i,
                                    Q_sqrt_i,
                                    z_i=z_i)
        if callback:
            callback(i)


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
    for x_hat_prior, P_sqrt_prior in sqrt_kalman_filter_yield(y, H, R_sqrt, F, Q_sqrt, mu, PI_sqrt,
                                                              z=z, callback=callback):
        x_hat.append(x_hat_prior)
        P_sqrt.append(P_sqrt_prior)
    return SqrtFilterResult(x_hat, P_sqrt)


if __name__ == '__main__':
    I = 10
    N = 4
    M = 3

    sigma_R = 1e-2
    sigma_Q = 1e-2

    y = []
    H = []
    R_sqrt = []
    F = []
    Q_sqrt = []

    R = []
    Q = []

    mu = np.zeros(N)
    PI_sqrt = np.random.rand(N, N)

    PI = PI_sqrt @ PI_sqrt.T

    x = [np.random.randn(N)]

    for i in range(I):
        # measurements
        H.append(np.random.randn(M, N))
        R_sqrt.append(sigma_R * np.random.randn(M, M))
        R.append(R_sqrt[-1] @ R_sqrt[-1].T)
        v = R_sqrt[-1] @ np.random.randn(M)
        y.append(H[-1] @ x[-1] + v)
        # dynamics
        F.append(np.random.randn(N, N))
        F[-1] /= np.linalg.norm(F[-1])
        Q_sqrt.append(sigma_Q * np.random.randn(N, N))
        Q.append(Q_sqrt[-1] @ Q_sqrt[-1].T)
        u = Q_sqrt[-1] @ np.random.randn(N)
        x.append(F[-1] @ x[-1] + u)


    kf_result = kalman_filter(y, H, R, F, Q, mu, PI)
    sqrt_kf_result = sqrt_kalman_filter(y, H, R_sqrt, F, Q_sqrt, mu, PI_sqrt)

    for i in range(I):
        assert np.allclose(kf_result.x_hat[i], sqrt_kf_result.x_hat[i])
        assert np.allclose(kf_result.P[i], sqrt_kf_result.P_sqrt[i] @ sqrt_kf_result.P_sqrt[i].T)

    from ..util.h5pymat import savemat

    m = {'I': I,
         'N': N,
         'D': math.ceil(math.log(I + 1)),
         'mu': mu,
         'PI': PI,
         'PI_sqrt': PI_sqrt}
    suffix_template = f"_{{:0{m['D']}d}}"
    for i in range(I):
        suffix = suffix_template.format(i)
        m[f'M{suffix}'] = len(y[i])
        m[f'y{suffix}'] = y[i]
        m[f'H{suffix}'] = H[i]
        m[f'R{suffix}'] = R[i]
        m[f'R_sqrt{suffix}'] = R_sqrt[i]
        m[f'F{suffix}'] = F[i]
        m[f'Q{suffix}'] = Q[i]
        m[f'Q_sqrt{suffix}'] = Q_sqrt[i]

    savemat('/tmp/test.mat', m)
