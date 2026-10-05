"""Tests for pyrsss.gnss.phase_edit (DiscFix log parsing and edits)."""
from collections import namedtuple
from datetime import timedelta

import pandas as pd
import pytest

from pyrsss.gnss.constants import GPS_EPOCH
from pyrsss.gnss.phase_edit import (ArcInfo,
                                    apply_phase_adjustments,
                                    apply_rejections,
                                    label_phase_arcs,
                                    parse_discfix_log)

Interval = namedtuple('Interval', 'lower upper')

DISCFIX_LOG_FIXTURE = """\
Some header line
Fine Arc 3 10 G01 9 0 1800 259200.000 1800 262800.000 3600.000 GC1W GL1C
Fine Arc 1 5 R04 4 0 1800 259200.000 1800 261000.000 1800.000 RC1P RL1C
"""


def test_parse_discfix_log(tmp_path):
    p = tmp_path / 'test.df.log'
    p.write_text(DISCFIX_LOG_FIXTURE)
    arcs = parse_discfix_log(str(p))
    assert len(arcs) == 2
    a = arcs[0]
    assert (a.gap, a.tot, a.sat, a.ok, a.s) == (3, 10, 'G01', 9, 0)
    assert a.start == GPS_EPOCH + timedelta(days=1800 * 7, seconds=259200)
    assert a.end == GPS_EPOCH + timedelta(days=1800 * 7, seconds=262800)
    assert a.dt == pytest.approx(3600.0)
    assert a.obs_types == 'GC1W GL1C'
    assert arcs[1].sat == 'R04'


def _dump_with_times():
    t0 = GPS_EPOCH + timedelta(days=1800 * 7, seconds=259200)
    return pd.DataFrame({'gps_time': [t0,
                                      t0 + timedelta(seconds=30),
                                      t0 + timedelta(seconds=60),
                                      t0 + timedelta(seconds=90)],
                         'sat': ['G01'] * 3 + ['G02'],
                         'L1': [1.0, 2.0, 3.0, 4.0]})


def test_label_phase_arcs():
    t0 = GPS_EPOCH + timedelta(days=1800 * 7, seconds=259200)
    phase_breaks = [
        ArcInfo(0, 3, 'G01', 3, 0,
                t0,
                t0 + timedelta(seconds=60),
                60.0, 'GC1W GL1C'),
    ]
    dump = _dump_with_times()
    # G02 row and G01 tail row (t0+90) fall outside all arcs and are dropped
    out = label_phase_arcs(dump, phase_breaks)
    assert set(out.arc) == {0}
    assert out.shape[0] == 3
    assert (out.sat == 'G01').all()


def test_apply_rejections():
    t0 = GPS_EPOCH + timedelta(days=1800 * 7, seconds=259200)
    dump = _dump_with_times()
    n0 = dump.shape[0]
    out = apply_rejections(dump, {'G01': [Interval(t0,
                                                   t0 + timedelta(seconds=30))]})
    assert out.shape[0] == n0 - 2


def test_apply_phase_adjustments():
    t0 = GPS_EPOCH + timedelta(days=1800 * 7, seconds=259200)
    dump = _dump_with_times()
    out = apply_phase_adjustments(dump, {'G01': [(t0 + timedelta(seconds=30),
                                                  'L1',
                                                  1.0)]})
    # adjustment applies at dt and after
    assert out.L1[0] == pytest.approx(1.0)
    assert out.L1[1] == pytest.approx(1.0)
    assert out.L1[2] == pytest.approx(2.0)
    # other satellites untouched
    assert out.L1[3] == pytest.approx(4.0)
