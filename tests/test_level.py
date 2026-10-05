"""Leveling tests for pyrsss.gnss.level."""
from datetime import timedelta

import numpy as np
import pytest

from pyrsss.gnss.constants import GPS_EPOCH, LAMBDA_1, LAMBDA_2, M_TO_TECU
from pyrsss.gnss.level import LeveledArc, convert_phase_m, level
from pyrsss.gnss.rinex import RinexDump

L_TRUE = 5.5  # leveling constant [m]


def synthetic_arc(n=40,
                  dt_s=30,
                  el=45.0,
                  noise=0.0,
                  l_true=L_TRUE,
                  sat='G01',
                  arc=0,
                  rng=None):
    """
    Build a synthetic geometry-free arc with phase level *l_true* [m]:
    P_I = P2 - P1 and L_Im = L1m - L2m = P_I - l_true (plus noise).
    """
    rng = rng or np.random.default_rng(42)
    t0 = GPS_EPOCH + timedelta(days=1800 * 7, seconds=259200)
    times = [t0 + timedelta(seconds=dt_s * i) for i in range(n)]
    P_I = 100.0 + 0.01 * np.arange(n) + noise * rng.standard_normal(n)
    L_Im = P_I - l_true + noise * rng.standard_normal(n)
    P1 = 21e6 * np.ones(n)
    P2 = P1 + P_I
    L2m = np.zeros(n)
    L1m = L_Im + L2m
    df = RinexDump({'gps_time': times,
                    'sat': [sat] * n,
                    'C1': P1,
                    'P1': P1,
                    'P2': P2,
                    'L1': L1m / LAMBDA_1,
                    'L2': L2m / LAMBDA_2,
                    'el': el * np.ones(n),
                    'az': 180.0 * np.ones(n),
                    'satx': 1000.0 * np.ones(n),
                    'saty': 2000.0 * np.ones(n),
                    'satz': 3000.0 * np.ones(n),
                    'arc': arc * np.ones(n, dtype=int)})
    df.xyz = [1.0, 2.0, 3.0]
    df.llh = [40.0, -88.0, 200.0]
    df.stn = 'TEST'
    df.recv_type = 'TEST'
    df.recv_p1c1 = 1
    df.p1c1_table = {}
    return df


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
