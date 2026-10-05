"""Tests for pyrsss.gnss.constants, pyrsss.gnss.rinex_new."""
from datetime import timedelta

import pandas as pd
import pytest

from pyrsss.gnss.constants import GPS_EPOCH, week_sec2dt
from pyrsss.gnss.rinex_new import RinexDump, correct_p1c1


def test_week_sec2dt():
    assert week_sec2dt(0, 0) == GPS_EPOCH
    assert week_sec2dt(1, 0) == GPS_EPOCH + timedelta(days=7)
    assert week_sec2dt(1800, 30.5) == GPS_EPOCH + timedelta(days=1800 * 7,
                                                            seconds=30.5)


RINDUMP_FIXTURE = """\
# Data GC1C GC1W GL1C GC2W GL2W ELE AZI SVX SVY SVZ
# Refpos XYZ(m): 1000.0 2000.0 3000.0  LLH(ddm): 40.0N 272.0E 200.0
# Station ID: JPLM
# Receiver type: SEPT POLARX2
# Receiver p1c1 type: 2
# P1-C1 [m]: G01: 0.100
1800 259200.0 G01 21000000.1 21000000.2 110000000.0 210000100.2 88000000.0 45.0 180.0 100.0 200.0 300.0
1800 259230.0 G01 21000000.3 21000000.4 110000001.0 210000100.4 88000001.0 45.1 181.0 101.0 201.0 301.0
"""


def test_load(tmp_path):
    p = tmp_path / 'test.dump'
    p.write_text(RINDUMP_FIXTURE)
    dump = RinexDump.load(str(p), p1c1=False)
    assert list(dump.columns) == ['gps_time', 'sat', 'C1', 'P1', 'L1',
                                  'P2', 'L2', 'el', 'az',
                                  'satx', 'saty', 'satz']
    assert dump.shape[0] == 2
    assert dump.stn == 'JPLM'
    assert dump.xyz == [1000.0, 2000.0, 3000.0]
    assert dump.llh == [40.0, 272.0, 200.0]
    assert dump.recv_type == 'SEPT POLARX2'
    assert dump.recv_p1c1 == 2
    assert dump.p1c1_table == {'G01': 0.100}
    assert dump.iloc[0].gps_time == week_sec2dt(1800, 259200.0)
    assert dump.iloc[0].C1 == pytest.approx(21000000.1)


def test_correct_p1c1_receiver_type_2():
    dump = pd.DataFrame({'sat': ['G01', 'G01'],
                         'C1': [21.0, 21.5],
                         'P1': [22.0, float('nan')],
                         'P2': [23.0, 23.5]},
                        index=[0, 1])
    dump.recv_p1c1 = 2
    dump.p1c1_table = {'G01': 0.5}
    correct_p1c1(dump)
    # type 2: C1 corrected only
    assert dump.C1[0] == pytest.approx(21.5)
    assert dump.C1[1] == pytest.approx(22.0)
    assert dump.P2[0] == pytest.approx(23.0)
    # missing P1 filled from corrected C1
    assert dump.P1[0] == pytest.approx(22.0)
    assert dump.P1[1] == pytest.approx(22.0)


def test_correct_p1c1_receiver_type_1():
    dump = pd.DataFrame({'sat': ['G01'],
                         'C1': [21.0],
                         'P1': [22.0],
                         'P2': [23.0]})
    dump.recv_p1c1 = 1
    dump.p1c1_table = {'G01': 0.5}
    correct_p1c1(dump, replace_p1_with_c1=False)
    # type 1: C1 and P2 corrected
    assert dump.C1[0] == pytest.approx(21.5)
    assert dump.P2[0] == pytest.approx(23.5)
    assert dump.P1[0] == pytest.approx(22.0)


def test_correct_p1c1_unknown_receiver_type():
    dump = pd.DataFrame({'sat': ['G01'], 'C1': [21.0], 'P1': [22.0], 'P2': [23.0]})
    dump.recv_p1c1 = 9
    dump.p1c1_table = {'G01': 0.5}
    with pytest.raises(ValueError):
        correct_p1c1(dump)
