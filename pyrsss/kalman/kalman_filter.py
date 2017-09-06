from __future__ import division

from itertools import izip, repeat

import numpy as NP
from scipy.linalg import cho_factor, cho_solve


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

    Return the tuple of lists of posterior state estimates and error
    covariances.
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
    return x_hat, P
