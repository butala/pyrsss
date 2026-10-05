"""Shared synthetic-data builders for the test suite."""
from datetime import timedelta

import numpy as np

from pyrsss.gnss.constants import GPS_EPOCH, LAMBDA_1, LAMBDA_2
from pyrsss.gnss.level import LeveledArc
from pyrsss.gnss.rinex import RinexDump

L_TRUE = 5.5  # leveling constant [m]


def synthetic_rinex_dump(n=40,
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


# backwards-compatible alias
synthetic_arc = synthetic_rinex_dump


def synthetic_leveled_arc(sat='G01', n=10, el=45.0, L_I=None, P_I=None):
    """Build a synthetic :class:`LeveledArc` (post-leveling schema)."""
    t0 = GPS_EPOCH + timedelta(days=1800 * 7, seconds=259200)
    arc = LeveledArc({'gps_time': [t0 + timedelta(seconds=30 * i)
                                   for i in range(n)],
                      'az': np.full(n, 180.0),
                      'el': np.full(n, el),
                      'satx': np.full(n, 1000.0),
                      'saty': np.full(n, 2000.0),
                      'satz': np.full(n, 3000.0),
                      'P_I': P_I if P_I is not None else np.full(n, 10.0),
                      'L_I': L_I if L_I is not None else np.full(n, 5.0)})
    arc.xyz = [1.0, 2.0, 3.0]
    arc.llh = [40.0, -88.0, 200.0]
    arc.stn = 'TEST'
    arc.recv_type = 'TEST'
    arc.sat = sat
    arc.L = 5.0
    arc.L_scatter = 0.1
    return arc
