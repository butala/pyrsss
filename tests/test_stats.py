"""Tests for pyrsss.stats (weighted statistics, ARMA, SOS)."""
import numpy as np
import pytest
import scipy as sp

from pyrsss.stats.arma import ARMA, arma_sensitivity
from pyrsss.stats.sos import Parallel, Series
from pyrsss.stats.stats import weighted_avg_and_std


def test_weighted_avg_and_std_uniform():
    avg, std = weighted_avg_and_std(np.array([1.0, 2.0, 3.0]), np.ones(3))
    assert avg == pytest.approx(2.0)
    # biased estimator (division by N)
    assert std == pytest.approx(np.sqrt(2 / 3))


def test_weighted_avg_and_std_weighted():
    values = np.array([0.0, 10.0])
    weights = np.array([3.0, 1.0])
    avg, std = weighted_avg_and_std(values, weights)
    assert avg == pytest.approx(2.5)
    assert std == pytest.approx(np.sqrt(18.75))


def _finite_difference(y_of, theta, eps=1e-6):
    dtype = complex if np.iscomplexobj(theta) else float
    J = np.empty((len(y_of(theta)), len(theta)), dtype=dtype)
    for j in range(len(theta)):
        e = np.zeros_like(theta)
        e[j] = eps
        J[:, j] = (y_of(theta + e) - y_of(theta - e)) / (2 * eps)
    return J


def test_arma_sensitivity_finite_difference():
    rng = np.random.default_rng(0)
    x = rng.standard_normal(200)
    a = np.array([1.0, -0.7])
    b = np.array([0.5, 0.2])
    J = arma_sensitivity(b, a, x, 0)
    assert J.shape == (200, (len(a) - 1) + len(b))

    def y_of(theta):
        a_t = np.r_[1.0, theta[:len(a) - 1]]
        b_t = theta[len(a) - 1:]
        return sp.signal.lfilter(b_t, a_t, x)

    theta = np.r_[a[1:], b]
    J_fd = _finite_difference(y_of, theta)
    np.testing.assert_allclose(J, J_fd, rtol=1e-5, atol=1e-6)


def test_arma_sensitivity_with_delay():
    rng = np.random.default_rng(1)
    x = rng.standard_normal(150)
    a = np.array([1.0, -0.7, 0.1])
    b = np.array([0.0, 0.5, 0.2])  # Nk=1: leading zero
    Nk = 1
    J = arma_sensitivity(b, a, x, Nk)
    assert J.shape == (150, (len(a) - 1) + (len(b) - Nk))

    def y_of(theta):
        a_t = np.r_[1.0, theta[:len(a) - 1]]
        b_t = np.r_[np.zeros(Nk), theta[len(a) - 1:]]
        return sp.signal.lfilter(b_t, a_t, x)

    theta = np.r_[a[1:], b[Nk:]]
    J_fd = _finite_difference(y_of, theta)
    np.testing.assert_allclose(J, J_fd, rtol=1e-5, atol=1e-6)


def test_arma_sensitivity_complex():
    rng = np.random.default_rng(2)
    x = rng.standard_normal(100) + 1j * rng.standard_normal(100)
    a = np.array([1.0, 0.2j])
    b = np.array([0.5 - 0.1j, 0.2])
    J = arma_sensitivity(b, a, x, 0)
    assert np.iscomplexobj(J)

    def y_of(theta):
        a_t = np.r_[1.0, theta[:1]]
        b_t = theta[1:]
        return sp.signal.lfilter(b_t, a_t, x)

    theta = np.r_[a[1:], b]
    J_fd = _finite_difference(y_of, theta, eps=1e-7)
    np.testing.assert_allclose(J, J_fd, rtol=1e-4, atol=1e-5)


def test_arma_class_roundtrip_and_apply():
    b = np.array([0.5, 0.2])
    a = np.array([1.0, -0.7])
    model = ARMA(b, a)
    assert model.Na == 1
    assert model.Nb == 2
    theta = model.theta
    model2 = ARMA.from_theta(theta, Nb=2, Na=1)
    np.testing.assert_allclose(model2.theta, theta)
    rng = np.random.default_rng(3)
    x = rng.standard_normal(50)
    np.testing.assert_allclose(model(x), sp.signal.lfilter(b, a, x))
    np.testing.assert_allclose(model.residual(x, np.zeros(50)),
                               -sp.signal.lfilter(b, a, x))


def test_arma_fit_nonlinear_at_truth():
    b = np.array([0.5, 0.2])
    a = np.array([1.0, -0.7])
    rng = np.random.default_rng(4)
    x = rng.standard_normal(500)
    y = sp.signal.lfilter(b, a, x)
    theta_true = ARMA(b, a).theta
    fitted = ARMA.fit_nonlinear(x, y, Nb=2, Na=1, theta0=theta_true)
    np.testing.assert_allclose(fitted.theta, theta_true, rtol=1e-6, atol=1e-8)


def test_sos_from_theta_roundtrip():
    theta = np.array([0.5, -0.3, 0.1, 0.2, 0.05, -0.4, 0.3])
    Nb, Na = 4, 3
    for cls in (Series, Parallel):
        obj = cls.from_theta(theta, Nb, Na)
        np.testing.assert_allclose(obj.theta, theta)


def test_sos_series_polynomial_identity():
    theta = np.array([0.5, -0.3, 0.1, 0.2, 0.05, -0.4, 0.3])
    s = Series.from_theta(theta, 4, 3)
    # cascade of biquads == convolution of the individual polynomials
    b = np.r_[1.0]
    a = np.r_[1.0]
    for i in range(s.sos.shape[0]):
        b = np.convolve(b, s.sos[i, :3])
        a = np.convolve(a, s.sos[i, 3:])
    np.testing.assert_allclose(s.b, b)
    np.testing.assert_allclose(s.a, a)


def test_sos_series_call_equals_sosfilt():
    theta = np.array([0.5, -0.3, 0.1, 0.2, 0.05, -0.4, 0.3])
    s = Series.from_theta(theta, 4, 3)
    rng = np.random.default_rng(5)
    x = rng.standard_normal(80)
    y = s(x)
    expected = sp.signal.sosfilt(s.sos, x)
    np.testing.assert_allclose(y, expected)


def test_sos_parallel_call_equals_sum_of_sections():
    theta = np.array([0.5, -0.3, 0.1, 0.2, 0.05, -0.4, 0.3])
    p = Parallel.from_theta(theta, 4, 3)
    rng = np.random.default_rng(6)
    x = rng.standard_normal(80)
    expected = np.zeros_like(x)
    for i in range(p.sos.shape[0]):
        expected += sp.signal.lfilter(p.sos[i, :3], p.sos[i, 3:], x)
    np.testing.assert_allclose(p(x), expected)


def test_sos_jacobians_finite_difference():
    rng = np.random.default_rng(7)
    u = rng.standard_normal(120)
    theta = np.array([0.4, -0.2, 0.1, 0.1, 0.05, -0.2, 0.15])
    Nb, Na = 4, 3
    for cls in (Series, Parallel):
        model = cls.from_theta(theta, Nb, Na)
        J = model.jacobian(u)
        assert J.shape == (120, Nb + Na)

        def y_of(th):
            return cls.from_theta(th, Nb, Na)(u)

        J_fd = _finite_difference(y_of, theta, eps=1e-6)
        np.testing.assert_allclose(J, J_fd, rtol=1e-4, atol=1e-5)
