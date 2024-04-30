import sys
import math
import logging
from dataclasses import astuple
from collections import defaultdict
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import numpy as np

from pyrsss.kalman.kalman_filter import kalman_filter, sqrt_kalman_filter


def randn_pd(N):
    """
    Generate a random NxN positive definite matrix with norm 1. Return
    the tuple containing this matrix and its square root.
    """
    A_sqrt = np.random.randn(N, N)
    A = np.dot(A_sqrt, A_sqrt.T)
    norm = np.linalg.norm(A)
    return A / norm, A_sqrt / math.sqrt(norm)


def setup_random_sim(N, M, I):
    """
    Generate a random, standard state-space form simulation for the
    purpose of testing. The state dimension is *N*, number of
    measurements is *M*, and *I* is the number of time
    steps. State-space model parameters are generated randomly.
    """
    sim = defaultdict(list)
    sim['N'], sim['M'], sim['I'] = N, M, I
    sim['mu'] = np.random.randn(N)
    sim['PI'], sim['PI_sqrt'] = randn_pd(N)
    for _ in range(I):
        # measurement
        sim['H'].append(np.random.randn(M, N))
        R_i, R_i_sqrt = randn_pd(M)
        sim['R'].append(R_i)
        sim['R_sqrt'].append(R_i_sqrt)
        # time update
        sim['F'].append(np.random.randn(N, N))
        sim['F'][-1] /= np.linalg.norm(sim['F'][-1])
        Q_i, Q_i_sqrt = randn_pd(N)
        sim['Q'].append(Q_i)
        sim['Q_sqrt'].append(Q_i_sqrt)
    return sim


def run_random_sim(sim, L):
    """
    Run *L* simulations of the state space model specified by *sim*
    (see :func:`setup_random_sim`). Each simulation is added to *sim*
    index by an integer identifier.
    """
    sim['L'] = L
    for l in range(L):
        sim[l] = defaultdict(list)
        x_i = sim['mu'] + np.matmul(sim['PI_sqrt'], np.random.randn(sim['N']))
        for i in range(sim['I']):
            sim[l]['x'].append(x_i)
            # measurement
            v_i = np.matmul(sim['R_sqrt'][i], np.random.randn(sim['M']))
            sim[l]['y'].append(np.matmul(sim['H'][i], sim[l]['x'][i]) + v_i)
            # time update
            u_i = np.matmul(sim['Q_sqrt'][i], np.random.randn(sim['N']))
            x_i = np.matmul(sim['F'][i], x_i) + u_i
    return sim


def setup_random_test(N, M, I, L):
    """
    Generate a random, standard state-space form simulation for the
    purpose of testing. The state dimension is *N*, number of
    measurements is *M*, number of time steps is *I*, and *L* is the
    number of simulations. State-space model parameters are generated
    randomly. The random simulations are returned as a mapping between
    parameter names and integer simulation indices to lists of vectors
    or matrices.
    """
    sim = setup_random_sim(N, M, I)
    return run_random_sim(sim, L)


def kf_sim(sim):
    """
    Process each simulation trial generated with
    :func:`setup_random_test` with a Kalman filter and return the
    posterior state estimates and error covariances.
    """
    post = defaultdict(dict)
    for l in range(sim['L']):
        x_hat_l, P_l = astuple(kalman_filter(sim[l]['y'],
                                             sim['H'],
                                             sim['R'],
                                             sim['F'],
                                             sim['Q'],
                                             sim['mu'],
                                             sim['PI']))
        post[l]['x_hat'] = x_hat_l
        if l == 0:
            post['P'] = P_l
        post[l]['error'] = []
        for x_i, x_hat_i in zip(sim[l]['x'], post[l]['x_hat']):
            post[l]['error'].append(x_hat_i - x_i)
    return post


def sqrt_kf_sim(sim):
    """
    Process each simulation trial generated with
    :func:`setup_random_test` with a Kalman filter and return the
    posterior state estimates and error covariances.
    """
    post = defaultdict(dict)
    for l in range(sim['L']):
        x_hat_l, P_sqrt_l = astuple(sqrt_kalman_filter(sim[l]['y'],
                                                       sim['H'],
                                                       sim['R_sqrt'],
                                                       sim['F'],
                                                       sim['Q_sqrt'],
                                                       sim['mu'],
                                                       sim['PI_sqrt']))
        post[l]['x_hat'] = x_hat_l
        if l == 0:
            post['P'] = [np.matmul(x, x.T) for x in P_sqrt_l]
        post[l]['error'] = []
        for x_i, x_hat_i in zip(sim[l]['x'], post[l]['x_hat']):
            post[l]['error'].append(x_hat_i - x_i)
    return post


def analysis(N, M, I, L, sqrt=False):
    """
    Conduct a Kalman filter validation experiment. Output results
    (concerning the error, i.e., x_hat - x) for the last time step
    only. If *sqrt* use the square root form Kalman filter.
    """
    sim = setup_random_test(N, M, I, L)
    if sqrt:
        post = sqrt_kf_sim(sim)
    else:
        post = kf_sim(sim)
    # output statistics of \hat{x}_{I|I}
    error_I = []
    for l in range(sim['L']):
        error_I.append(post[l]['error'][-1])
    E_I = np.stack(error_I, 1)
    E_I_mean = np.mean(E_I, 1)
    P_I = np.cov(E_I)
    print('Mean of error at time step I={}'.format(I))
    for E_I_mean_n in E_I_mean:
        print('{:9.2e}'.format(E_I_mean_n))
    print('')
    print(f'Norm of mean of error = {np.linalg.norm(E_I_mean):.3e}')
    print('')
    print('True posterior covariance at time step I')
    print(np.array_str(post['P'][-1], precision=2))
    print('')
    print('Empirical posterior covariance at time step I')
    print(np.array_str(P_I, precision=2))
    print('')
    print(f"||P_{{I|I}}(true) - P_{{I|I}}(empirical)||_F = {np.linalg.norm(post['P'][-1] - P_I):.3e}")
    return sim, post


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Validate Kalman filter implementation with a randomly generated simulation.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('N',
                        type=int,
                        help='state dimension')
    parser.add_argument('M',
                        type=int,
                        help='number of measurements')
    parser.add_argument('I',
                        type=int,
                        help='number of time steps')
    parser.add_argument('L',
                        type=int,
                        help='number of simulation trials')
    parser.add_argument('--sqrt',
                        '-s',
                        action='store_true',
                        help='use square root form Kalman filter')
    parser.add_argument('--seed',
                        type=int,
                        default=None,
                        help='random number generator seed')
    args = parser.parse_args(argv[1:])

    np.random.seed(args.seed)

    analysis(args.N,
             args.M,
             args.I,
             args.L,
             sqrt=args.sqrt)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
