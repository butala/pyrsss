"""Tests for pyrsss.emission rate/coefficient formulas (Makela dissertation)."""
import math

import pytest

from pyrsss.emission.v6300 import (A_1D, A_6300, BETA_1D, Oplus,
                                   Oplus_simple, alpha1, alpha2, emission_v6300,
                                   k1, k2, k3, k4, k5)
from pyrsss.emission.v6300 import OplusType
from pyrsss.emission.v7774 import (ALPHA1, BETA_7774, K1, K2, K3,
                                   emission_v7774, emission_v7774_ii,
                                   emission_v7774_rr)


def test_alpha_rate_coefficients():
    assert alpha1(300.0) == pytest.approx(1.95e-7)
    assert alpha1(600.0) == pytest.approx(1.95e-7 * 2.0**(-0.7))
    assert alpha2(300.0) == pytest.approx(4.00e-7)
    assert alpha2(150.0) == pytest.approx(4.00e-7 * 0.5**(-0.9))


def test_k_rate_coefficients():
    assert k3(300.0) == pytest.approx(2.0e-11 * math.exp(111.8 / 300.0))
    assert k4(300.0) == pytest.approx(2.9e-11 * math.exp(67.5 / 300.0))
    assert k5(300.0) == pytest.approx(1.6e-12 * 300.0**0.91)
    assert k1(300.0) == pytest.approx(3.23e-12 * math.exp(3.72 - 1.87))
    assert k2(300.0) == pytest.approx(2.78e-13 * math.exp(2.07 - 0.61))


def test_oplus():
    assert Oplus_simple(1e5) == 1e5
    # zero neutrals -> unclustered O+
    assert Oplus(1e5, 300.0, 300.0, 0.0, 0.0) == pytest.approx(1e5)
    # neutral clustering reduces the O+ density
    assert Oplus(1e5, 300.0, 300.0, 1e8, 1e8) < 1e5


def test_emission_v6300():
    kwargs = dict(Te=300.0, Ti=300.0, Tn=300.0, O2=1e8, N2=1e8)
    rate_simple = emission_v6300(1e5, oplus_type=OplusType.ne, **kwargs)
    rate_full = emission_v6300(1e5, oplus_type=OplusType.charge_neutrality,
                               **kwargs)
    assert rate_simple > 0
    assert rate_full > 0
    assert rate_full < rate_simple  # clustering reduces emission
    with pytest.raises(NotImplementedError):
        emission_v6300(1e5, oplus_type='bogus', **kwargs)


def test_emission_constants():
    assert A_1D == pytest.approx(7.45e-3)
    assert A_6300 == pytest.approx(5.63e-3)
    assert BETA_1D == pytest.approx(1.1)


def test_emission_v7774():
    assert emission_v7774_rr(2.0, 0.0) == pytest.approx(ALPHA1 * 4.0)
    # no O+ -> no intercombination emission
    assert emission_v7774_ii(1e5, 1e8, 0.0) == pytest.approx(0.0)
    ii = emission_v7774_ii(1e5, 1e8, 1e5)
    assert ii == pytest.approx((BETA_7774 * K1 * K2 * 1e8 * 1e5 * 1e5) /
                               (K2 * 1e5 + K3 * 1e8))
    assert emission_v7774(1e5, 1e8, 1e5) == pytest.approx(
        emission_v7774_rr(1e5, 1e8) + ii)
