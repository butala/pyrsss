import math
from itertools import repeat

import numpy as NP
import pylab as PL

from kalman_filter import kalman_filter


if __name__ == '__main__':
    # true initial elevation angle and velocity of projectile
    el_deg_true = 45  # deg
    v0_true = 100     # m/s

    # true projectile launch position (m)
    px0_true = 0
    py0_true = 0

    el_rad_true = math.radians(el_deg_true)
    vx0_true = v0_true * math.cos(el_rad_true)
    vy0_true = v0_true * math.sin(el_rad_true)

    # gravitational constant (m/s^2)
    g_true = 9.81

    # time step (s)
    delta = 0.1

    # number of time steps
    N_STEPS = 144
    n = NP.arange(N_STEPS)

    # calculate projection positions and velocities
    px_true = px0_true + vx0_true * delta * n
    py_true = py0_true + vy0_true * delta * n - 1/2*g_true*(delta * n)**2

    vx_true = vx0_true * NP.ones(N_STEPS)
    vy_true = vy0_true - g_true * delta * n


    # compute true state (position x, velocity x, position y, velocity y)
    x = [NP.array([px, vx, py, vy]) for px, vx, py, vy in zip(px_true,
                                                              vx_true,
                                                              py_true,
                                                              vy_true)]

    # measurement operator (position x and y)
    H_i = NP.array([[1., 0, 0,  0],
                    [0,  0, 1., 0]])

    H = repeat(H_i)

    # measurement noise covariance
    var_v = 20
    R_i = var_v * NP.eye(2)
    R = repeat(R_i)

    # measurements without noise
    y_true = [NP.dot(H_i, x_i) for x_i in x]

    # measurements with noise
    y = [y_true_i + NP.matmul(R_i, NP.random.randn(2)) for y_true_i in y_true]

    # time update operator (based on ballistic motion)
    F_i = NP.array([[1, delta, 0,      0],
                    [0,     1, 0,      0],
                    [0,     0, 1,  delta],
                    [0,     0, 0,      1]])
    F = repeat(F_i)

    # time update covariance (no uncertainty)
    Q_i = NP.zeros((4, 4))
    Q = repeat(Q_i)

    # time update systematics (gravitation)
    z_i = NP.array([0,
                    0,
                    -1/2*g_true*delta**2,
                    -g_true*delta])
    z = repeat(z_i)

    # provided initial elevation angle and velocity of projectile
    el_deg_0 = 45  # deg
    v0 = 100       # m/s

    el_rad_0 = math.radians(el_deg_0)

    vx0 = v0 * math.cos(el_rad_0)
    vy0 = v0 * math.sin(el_rad_0)

    # provided initial position of projectile (m)
    px0 = 0
    py0 = 500

    mu = NP.array([px0, vx0, py0, vy0])

    # initial state variance parameters
    var_px0 = 10.
    var_py0 = 100.

    var_vx0 = 1.
    var_vy0 = 1.

    # initial state covariance
    PI = NP.array([[var_px0,       0,       0,       0],
                   [      0, var_vx0,       0,       0],
                   [      0,       0, var_py0,       0],
                   [      0,       0,       0, var_vy0]])

    # execute filter
    x_hat, P = kalman_filter(y, H, R, F, Q, mu, PI, z=z)

    # gather results and plot
    px_hat = [x_hat_i[0] for x_hat_i in x_hat]
    py_hat = [x_hat_i[2] for x_hat_i in x_hat]

    fig = PL.figure(figsize=(11, 4))
    ax = fig.add_subplot(111)
    PL.plot([x_i[0] for x_i in x],
            [x_i[2] for x_i in x],
            label='Truth')
    PL.plot(px_hat,
            py_hat,
            color='r',
            ls='-',
            mec='None',
            label='Posterior estimate')
    PL.scatter([y_i[0] for y_i in y],
               [y_i[1] for y_i in y],
               color='g',
               edgecolors='None',
               label='Measurement')
    PL.xlim(0, px_true[-1])
    PL.ylim(0, max(py_true) + 10)
    ax.set_aspect('equal')
    PL.title('Ballistic Tracking Example')
    PL.xlabel('Horizontal distance from cannon [m]')
    PL.ylabel('Vertical distance from cannon [m]')
    PL.legend()

    PL.show()
