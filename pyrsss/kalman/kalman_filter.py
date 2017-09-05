from __future__ import division

from itertools import izip

from scipy.linalg import cho_factor, cho_solve


def is_col_vector(x):
    """
    ???
    """
    return x.ndim == 2 and x.shape[1] == 1


def kalman_filter(y, H, R, F, Q, mu, PI, z=None):
    """
    ???
    """
    x_hat = []
    P = []
    assert is_col_vector(mu)
    x_hat_prior = mu
    P_prior = PI
    for i, (y_i, H_i, R_i, F_i, Q_i) in enumerate(izip(y, H, R, F, Q)):
        assert is_col_vector(y_i)
        # measurement update
        A = cho_factor(NP.dot(H_i, NP.dot(P_prior, H_i.T)) + R_i)
        B = cho_solve(A, NP.dot(H_i, P_prior))
        b = cho_solve(A, y_i - NP.dot(H_i, x_hat_prior))
        C = NP.dot(P_prior, H_i.T)
        x_hat.append(x_hat_prior + NP.dot(C, b))
        P.append(P_prior - NP.dot(C, B))
        # time update
        x_hat_prior = NP.dot(F_i, x_hat[-1])
        if z is not None:
            assert is_col_vector(z[i])
            x_hat_prior += z[i]
        P_prior = NP.dot(F_i, NP.dot(P[-1], F_i.T)) + Q_i
    return x_hat, P
