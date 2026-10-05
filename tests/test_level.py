"""Leveling tests for pyrsss.gnss.level."""
import numpy as np
import pytest

from pyrsss.gnss.constants import LAMBDA_1, LAMBDA_2, M_TO_TECU
from pyrsss.gnss.level import LeveledArc, convert_phase_m, level

from helpers import L_TRUE, synthetic_arc


def test_level_recovers_constant():
    arcs = level(synthetic_arc(noise=1e-4))
    assert len(arcs) == 1
    arc = arcs[0]
    assert isinstance(arc, LeveledArc)
    assert arc.L == pytest.approx(L_TRUE * M_TO_TECU, rel=1e-3)
    assert arc.L_scatter == pytest.approx(0.0, abs=1e-2)


def test_level_rejects_short_arc_time():
    # 2 min of data < 18 min minimum_arc_time
    assert level(synthetic_arc(n=5, dt_s=30)) == []


def test_level_rejects_few_points():
    # under minimum_arc_points with a long time span
    assert level(synthetic_arc(n=2, dt_s=1000)) == []


def test_level_rejects_low_elevation():
    assert level(synthetic_arc(el=5.0)) == []


def test_level_rejects_scatter():
    # noise far above the modeled RMS -> scatter rejection
    assert level(synthetic_arc(noise=50.0)) == []


def test_convert_phase_m_gps():
    df = synthetic_arc(n=3)
    L1m, L2m = convert_phase_m(df, 'G01')
    np.testing.assert_allclose(L1m, df.L1 * LAMBDA_1)
    np.testing.assert_allclose(L2m, df.L2 * LAMBDA_2)


def test_convert_phase_m_glonass(monkeypatch):
    monkeypatch.setattr('pyrsss.gnss.level.glonass_lambda',
                        lambda slot, dt: (0.187, 0.242))
    df = synthetic_arc(n=3, sat='R04')
    L1m, L2m = convert_phase_m(df, 'R04')
    np.testing.assert_allclose(L1m, df.L1 * 0.187)
    np.testing.assert_allclose(L2m, df.L2 * 0.242)
