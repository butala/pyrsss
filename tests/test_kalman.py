"""Tests for pyrsss.kalman.kalman_filter."""
import numpy as np
import pytest

from pyrsss.kalman.kalman_filter import kalman_filter


def test_scalar_two_step_hand_computed():
    """
    Scalar state/measurement with zero process noise: the posterior
    moments are available in closed form (precision-weighted mean).
    """
    y = [np.array([1.0]), np.array([2.0])]
    H = [np.eye(1), np.eye(1)]
    R = [np.eye(1), np.eye(1)]
    F = [np.eye(1), np.eye(1)]
    Q = [np.zeros((1, 1)), np.zeros((1, 1))]
    mu = np.zeros(1)
    PI = np.eye(1)
    res = kalman_filter(y, H, R, F, Q, mu, PI)
    # step 1: prior N(0, 1), measurement N(1, 1) -> N(1/2, 1/2)
    assert res.x_hat[0][0] == pytest.approx(0.5)
    assert res.P[0][0, 0] == pytest.approx(0.5)
    # step 2: prior N(1/2, 1/2), measurement N(2, 1) -> N(1, 1/3)
    assert res.x_hat[1][0] == pytest.approx(1.0)
    assert res.P[1][0, 0] == pytest.approx(1 / 3)


def test_constant_velocity_perfect_measurements():
    """
    Two-state constant-velocity model with tiny measurement noise: the
    filter tracks the deterministic trajectory.
    """
    n = 20
    dt = 1.0
    F1 = np.array([[1.0, dt], [0.0, 1.0]])
    H1 = np.array([[1.0, 0.0]])
    truth = [np.array([0.5 * i, 0.5]) for i in range(n)]
    y = [H1 @ x for x in truth]
    res = kalman_filter(y,
                        [H1] * n,
                        [np.eye(1) * 1e-12] * n,
                        [F1] * n,
                        [np.zeros((2, 2))] * n,
                        np.zeros(2),
                        np.eye(2) * 10)
    for i in range(5, n):
        assert res.x_hat[i][0] == pytest.approx(truth[i][0], abs=0.05)
        assert res.x_hat[i][1] == pytest.approx(truth[i][1], abs=0.05)


def test_covariance_shrinks_with_measurements():
    y = [np.array([0.0])] * 5
    res = kalman_filter(y,
                        [np.eye(1)] * 5,
                        [np.eye(1)] * 5,
                        [np.eye(1)] * 5,
                        [np.zeros((1, 1))] * 5,
                        np.zeros(1),
                        np.eye(1) * 10)
    variances = [P[0, 0] for P in res.P]
    assert all(np.diff(variances) < 0)
    assert variances[-1] < variances[0]
