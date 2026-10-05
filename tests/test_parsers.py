"""Fixture-based parser tests (IAGA-2002, Bartels, GLONASS, P1C1, PRN/SVN)."""
from datetime import date, datetime
from pathlib import Path

import pytest

from pyrsss.gnss.glonass import GLONASS_Status
from pyrsss.gnss.p1c1 import P1C1Table
from pyrsss.gnss.prn_gps import Table as PrnGpsTable
from pyrsss.l1.bartels import Bartels, efloat, parse as bartels_parse
from pyrsss.mag.iaga2002 import convert_float, fname2date, parse as iaga_parse

FIXTURES = Path(__file__).parent / 'fixtures'


def test_efloat():
    assert efloat('1.5') == (1.5, False)
    assert efloat('e2.5') == (2.5, True)
    assert efloat('-') == (None, False)


def test_bartels_parse(tmp_path):
    p = tmp_path / 'bartels.txt'
    p.write_text("""\
some header line
-----
1800 12/26/1993 360 1234.5 12345.6 e123.4
1801 01/15/1994 15 100.0 20000.0 200.5
1802 02/03/1994 34 200.0 30000.0 -
""")
    data_map = bartels_parse(str(p))
    assert set(data_map) == {1800, 1801, 1802}
    field = data_map[1800]
    assert field.start == datetime(1993, 12, 26)
    assert field.doy == 360
    assert field.ut_sec == pytest.approx(1234.5)
    assert field.spacecraft_clock == pytest.approx(123.4)
    assert field.spacecraft_clock_estimated is True
    assert data_map[1802].spacecraft_clock is None

    b = Bartels(str(p))
    # rotation lookup by the interval midpoints
    assert b(datetime(1994, 1, 1)) == 1800
    assert b(datetime(1994, 1, 20)) == 1801
    assert b(datetime(2020, 1, 1)) is None


def test_glonass_status(tmp_path):
    p = tmp_path / 'GLO_STATUS'
    p.write_text("""\
# comment
2000-01-01 00:00 2001-06-01 12:00 2005-01-01 00:00 1 23 4 5 6
2000-01-01 00:00 2003-01-01 00:00 0000-00-00 00:00 2 -3 4 5 6
""")
    status = GLONASS_Status(str(p))
    info = status(1, datetime(2002, 1, 1))
    assert info.slot == 1
    assert info.freq == 23
    assert info.plane == 4
    # open end date (0000-00-00) is unbounded
    assert status(2, datetime(2020, 1, 1)).freq == -3
    # outside the valid interval
    with pytest.raises(KeyError):
        status(1, datetime(2010, 1, 1))


def test_p1c1_table(tmp_path):
    p = tmp_path / 'p1c1'
    p.write_text("""\
# comment
2014-01-01 1 2 0.123
2014-01-01 2 3 -0.5
2014-06-01 1 2 0.2
""")
    table = P1C1Table(str(p))
    assert table[datetime(2014, 1, 1)]['prn'][1] == pytest.approx(0.123)
    assert table[datetime(2014, 1, 1)]['svn'][3] == pytest.approx(-0.5)
    # nearest entry within the default 32 day window
    assert table(datetime(2014, 1, 15))['prn'][1] == pytest.approx(0.123)
    # too far from any entry
    with pytest.raises(AssertionError):
        table(datetime(2015, 1, 1))


def test_prn_gps_table(tmp_path):
    p = tmp_path / 'PRN_GPS'
    p.write_text("""\
header line
2010-01-01 2015-06-01 63 1 IIR-M
2012-03-01 0000 64 2 IIR-M 2
""")
    table = PrnGpsTable(str(p))
    assert table.prn(63, date(2012, 1, 1)) == 1
    assert table.svn(1, date(2012, 1, 1)) == 63
    # active satellite (0000 deactivation date)
    assert table.prn(64, date(2020, 1, 1)) == 2
    with pytest.raises(RuntimeError):
        table.prn(63, date(2020, 1, 1))


def test_iaga2002_parse():
    fname = str(FIXTURES / 'hon20000807d.min')
    assert fname2date(fname) == datetime(2000, 8, 7)
    header, data_map = iaga_parse(fname)
    assert header.IAGA_CODE == 'HON'
    assert header.Station_Name == 'Honolulu'
    assert header.Reported == 'HDZF'
    assert header.Geodetic_Latitude == pytest.approx(21.301)
    assert header.Geodetic_Longitude == pytest.approx(202.0)
    assert len(data_map) > 100
    assert min(data_map).date() == date(2000, 8, 7)
    rec = data_map[min(data_map)]
    for value in (rec.H, rec.D, rec.z, rec.f):
        assert isinstance(value, float)


def test_iaga2002_convert_float():
    assert convert_float('1.5') == pytest.approx(1.5)
    assert convert_float('-2.25') == pytest.approx(-2.25)
