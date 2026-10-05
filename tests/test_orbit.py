"""Tests for pyrsss.gnss.orbit (broadcast ephemeris propagation, az/el)."""
import math
from datetime import datetime, timedelta

import numpy as np
import pytest

from pyrsss.gnss.orbit import (GPS_EPOCH, MU_GLO, OMEGA_E, azel,
                               glonass_position, gps_position,
                               satellite_position)
from pyrsss.util.position import Position


def _circular_gps_eph(A=26560e3):
    """Minimal circular, equator-crossing GPS ephemeris (Omega0=0, M0=0)."""
    return dict(sqrtA=math.sqrt(A), DeltaN=0.0, M0=0.0, Cuc=0.0,
                Eccentricity=0.0, Cus=0.0, Toe=0.0, Cic=0.0, Omega0=0.0,
                Cis=0.0, Io=math.radians(55.0), Crc=0.0, Crs=0.0,
                omega=0.0, OmegaDot=0.0, IDOT=0.0, GPSWeek=1800)


def test_gps_position_closed_form():
    # at tk=0 with these elements the position is exactly (A, 0, 0)
    A = 26560e3
    eph = _circular_gps_eph(A)
    p = gps_position(eph, GPS_EPOCH + timedelta(weeks=1800))
    np.testing.assert_allclose(p, [A, 0.0, 0.0], atol=1e-6)


def test_gps_position_radius_stable():
    eph = _circular_gps_eph()
    for mins in (0, 30, 60, 120):
        p = gps_position(eph, GPS_EPOCH + timedelta(weeks=1800, minutes=mins))
        assert np.linalg.norm(p) == pytest.approx(26560e3, rel=1e-5)


def test_gps_position_eccentric():
    eph = _circular_gps_eph()
    eph['Eccentricity'] = 0.01
    eph['M0'] = 0.5
    t0 = GPS_EPOCH + timedelta(weeks=1800)
    p0 = gps_position(eph, t0)
    p1 = gps_position(eph, t0 + timedelta(hours=1))
    # well inside the MEO shell at both epochs
    for p in (p0, p1):
        assert 25000e3 < np.linalg.norm(p) < 28000e3
    assert np.linalg.norm(p1 - p0) > 1000e3  # moved along the orbit


def _circular_glo_geph():
    r = 25510e3
    v_rot = math.sqrt(MU_GLO / r) - OMEGA_E * r
    geph = dict(X=25510.0, Y=0.0, Z=0.0, VX=0.0, VY=v_rot / 1e3, VZ=0.0,
                AX=0.0, AY=0.0, AZ=0.0)
    geph['_epoch'] = datetime(2014, 6, 7)
    return geph


def test_glonass_position_exact_at_epoch():
    geph = _circular_glo_geph()
    p = glonass_position(geph, geph['_epoch'])
    np.testing.assert_allclose(p, [25510e3, 0.0, 0.0])


def test_glonass_position_bounded():
    geph = _circular_glo_geph()
    p = glonass_position(geph, geph['_epoch'] + timedelta(minutes=15))
    assert np.linalg.norm(p) == pytest.approx(25510e3, rel=2e-5)


def test_glonass_position_velocity():
    geph = _circular_glo_geph()
    p0 = glonass_position(geph, geph['_epoch'])
    p1 = glonass_position(geph, geph['_epoch'] + timedelta(seconds=1))
    v = math.hypot(0.0, geph['VY'] * 1e3)
    assert np.linalg.norm(p1 - p0) == pytest.approx(v, rel=5e-3)


def test_satellite_position_dispatch():
    eph = _circular_gps_eph()
    t = GPS_EPOCH + timedelta(weeks=1800)
    p = satellite_position({'G01': eph}, 'G01', t)
    np.testing.assert_allclose(p, gps_position(eph, t))
    with pytest.raises(ValueError):
        satellite_position({}, 'E01', t)


def test_azel_zenith_and_east():
    stn = Position(0.0, 0.0, 0.0, Position.CoordinateSystem['geodetic'])
    # straight up
    up = np.array(stn.xyz) * (1 + 1e6 / np.linalg.norm(stn.xyz))
    az, el = azel(stn.llh, up)
    assert el == pytest.approx(90.0, abs=1e-6)
    # due east (+Y at the equator/prime meridian): az=90, el=0
    east = np.array(stn.xyz) + np.array([0.0, 2e7, 0.0])
    az, el = azel(stn.llh, east)
    assert az == pytest.approx(90.0, abs=1e-9)
    assert el == pytest.approx(0.0, abs=1e-9)
