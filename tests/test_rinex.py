"""Tests for pyrsss.gnss.rinex (RinexDump) and pyrsss.gnss.constants."""
from datetime import timedelta

import numpy as np
import pytest

from pyrsss.gnss.constants import GPS_EPOCH, week_sec2dt
from pyrsss.gnss.rinex import RinexDump, correct_p1c1
from pyrsss.util.position import Position


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

# same observables as RINDUMP_FIXTURE as a RINEX 2 OBS record
_OBS_VALUES = [(21000000.1, 21000000.2, 110000000.0, 210000100.2, 88000000.0),
               (21000000.3, 21000000.4, 110000001.0, 210000100.4, 88000001.0)]


def _hdr(content, label):
    return f'{content:<60}{label}'


def _obs_line(dt, sats, flag=0):
    return (f' {dt.year % 100:2d} {dt.month:2d} {dt.day:2d} {dt.hour:2d}'
            f' {dt.minute:2d}{dt.second + dt.microsecond / 1e6:11.7f}'
            f'  {flag:d}{len(sats):3d}' + ''.join(sats))


def _obs_vals(vals):
    return ''.join(f'{v:14.3f}  ' for v in vals)


def _obs_fixture_text():
    t0 = week_sec2dt(1800, 259200.0)
    t1 = week_sec2dt(1800, 259230.0)
    return '\n'.join([
        _hdr('     2.11           OBSERVATION DATA    G (GPS)',
             'RINEX VERSION / TYPE'),
        _hdr('JPLM', 'MARKER NAME'),
        _hdr('  1000.0000        2000.0000        3000.0000',
             'APPROX POSITION XYZ'),
        _hdr('                    SEPT POLARX2', 'REC # / TYPE'),
        _hdr('     5    C1    P1    L1    P2    L2', '# / TYPES OF OBSERV'),
        _hdr('', 'END OF HEADER'),
        _obs_line(t0, ['G01']),
        _obs_vals(_OBS_VALUES[0]),
        _obs_line(t1, ['G01']),
        _obs_vals(_OBS_VALUES[1]),
    ]) + '\n'


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


def test_load_with_p1c1(tmp_path):
    p = tmp_path / 'test.dump'
    p.write_text(RINDUMP_FIXTURE)
    dump = RinexDump.load(str(p), p1c1=True)
    # receiver type 2: C1 corrected by +0.1 and missing P1 filled from C1
    assert dump.iloc[0].C1 == pytest.approx(21000000.2)
    assert dump.iloc[0].P1 == pytest.approx(21000000.2)
    # P2 untouched by type 2
    assert dump.iloc[0].P2 == pytest.approx(210000100.2)


def test_from_rinex_parity_with_load(tmp_path):
    """
    RinexDump.from_rinex and RinexDump.load must agree on the shared
    schema given equivalent input (observables + station metadata).
    """
    p = tmp_path / 'jplm0010.14o'
    p.write_text(_obs_fixture_text())
    d = tmp_path / 'test.dump'
    d.write_text(RINDUMP_FIXTURE)
    from_rinex = RinexDump.from_rinex(str(p), p1c1=False)
    from_load = RinexDump.load(str(d), p1c1=False)
    shared = ['gps_time', 'sat', 'C1', 'P1', 'L1', 'P2', 'L2']
    for column in shared:
        if column == 'sat':
            assert from_rinex[column].tolist() == from_load[column].tolist()
        elif column == 'gps_time':
            assert from_rinex[column].tolist() == from_load[column].tolist()
        else:
            np.testing.assert_allclose(from_rinex[column],
                                       from_load[column], rtol=1e-12)
    assert from_rinex.stn == from_load.stn == 'JPLM'
    assert from_rinex.recv_type == from_load.recv_type == 'SEPT POLARX2'
    assert from_rinex.xyz == from_load.xyz == [1000.0, 2000.0, 3000.0]


def test_from_rinex_no_nav(tmp_path):
    p = tmp_path / 'jplm0010.14o'
    p.write_text(_obs_fixture_text())
    dump = RinexDump.from_rinex(str(p), p1c1=False)
    assert dump.shape[0] == 2
    # geometry columns require the NAV file
    assert np.isnan(dump['el']).all()
    assert np.isnan(dump['satx']).all()
    # llh derived from the header position (WGS84)
    np.testing.assert_allclose(dump.llh,
                               list(Position(1000.0, 2000.0, 3000.0).llh))
    assert dump.recv_p1c1 is None


def test_from_rinex_with_nav(tmp_path):
    obs = tmp_path / 'jplm0010.14o'
    obs.write_text(_obs_fixture_text())
    nav = tmp_path / 'jplm0010.14n'

    def d19(v):
        return f'{v:19.12E}'.replace('E', 'D')

    def nav_line1(prn, dt, v1, v2, v3):
        return (f'{prn:2d} {dt.year % 100:2d} {dt.month:2d} {dt.day:2d}'
                f' {dt.hour:2d} {dt.minute:2d}{dt.second:5.1f}'
                + d19(v1) + d19(v2) + d19(v3))

    def nav_cont(vals):
        return '   ' + ''.join(d19(v) for v in vals)

    clock = [1e-4, 2e-11, 0.0]
    # circular equatorial-ish elements with Omega0=0, M0=0 -> pos ~ (A, 0, 0)
    A = 26560e3
    import math
    orb25 = [13.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, math.sqrt(A),
             259200.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 0.0, 0.0,
             0.0, 0.0, 1800.0, 0.0,
             2.0, 0.0, 1e-3, 1.0,
             0.0]
    toc = week_sec2dt(1800, 259200.0)
    nav.write_text('\n'.join([
        _hdr('     2.11           N: GPS NAV DATA', 'RINEX VERSION / TYPE'),
        _hdr('', 'END OF HEADER'),
        nav_line1(1, toc, *clock),
        *[nav_cont(orb25[i:i + 4]) for i in range(0, 24, 4)],
    ]) + '\n')
    dump = RinexDump.from_rinex(str(obs), nav_fname=str(nav), p1c1=False)
    assert dump.shape[0] == 2
    for column in ('el', 'az', 'satx', 'saty', 'satz'):
        assert not dump[column].isnull().any()
    # at tk=0 the satellite sits at argument-of-latitude 0 in the plane
    # rotated by the ICD -omega_e*Toe term
    from pyrsss.gnss.orbit import OMEGA_E
    OmegaK = -OMEGA_E * 259200.0
    np.testing.assert_allclose([dump.iloc[0].satx, dump.iloc[0].saty,
                                dump.iloc[0].satz],
                               [A * math.cos(OmegaK),
                                A * math.sin(OmegaK),
                                0.0], rtol=1e-6, atol=1e3)
    # cross-check az/el against the geometry module
    from pyrsss.gnss.orbit import azel
    az, el = azel(dump.llh, [dump.iloc[0].satx, dump.iloc[0].saty,
                             dump.iloc[0].satz])
    assert dump.iloc[0].az == pytest.approx(az, abs=1e-9)
    assert dump.iloc[0].el == pytest.approx(el, abs=1e-9)


def test_correct_p1c1_receiver_type_2():
    dump = RinexDump({'sat': ['G01', 'G01'],
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
    dump = RinexDump({'sat': ['G01'],
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
    dump = RinexDump({'sat': ['G01'], 'C1': [21.0], 'P1': [22.0], 'P2': [23.0]})
    dump.recv_p1c1 = 9
    dump.p1c1_table = {'G01': 0.5}
    with pytest.raises(ValueError):
        correct_p1c1(dump)
