"""Tests for pyrsss.util.position (pure-Python WGS84 geometry)."""
import math

import pytest

from pyrsss.gnss.ipp import ipp_from_azel
from pyrsss.util.position import Position, point

# Millstone Hill IS Radar Madrigal validation case (see point.__doc__):
# az=100 [deg], range=1000 [km], expected (lat [deg], lon [deg], alt [km])
MADRIGAL = [(45, 41.361338, -64.012651, 742.038442),
            (55, 41.645771, -65.488196, 841.794521),
            (65, 41.932738, -67.098375, 918.473216)]


def _wrap180(lon):
    return lon - 360 if lon > 180 else lon


def test_point_madrigal_validation():
    stn = Position(42.619, 288.51, 0.146,
                   Position.CoordinateSystem['geodetic'])
    for el, elat, elon, ealt_km in MADRIGAL:
        lat, lon, alt = point(stn, 100, el, 1000).llh
        # agreement with the published Madrigal web service values
        assert abs(lat - elat) < 1e-3
        assert abs(_wrap180(lon) - elon) < 1e-3
        assert abs(alt / 1e3 - ealt_km) < 1.0


def test_roundtrip_geodetic():
    for lat, lon, alt in [(42.619, 288.51, 146.0),
                          (-60.0, 180.0, 0.0),
                          (0.0, 0.0, 450e3),
                          (89.0, 359.0, -100.0)]:
        pos = Position(lat, lon, alt,
                       Position.CoordinateSystem['geodetic'])
        lat2, lon2, alt2 = pos.llh
        assert lat2 == pytest.approx(lat, abs=1e-9)
        assert _wrap180(lon2) == pytest.approx(_wrap180(lon), abs=1e-9)
        assert alt2 == pytest.approx(alt, abs=1e-6)


def test_roundtrip_cartesian():
    xyz = (4696986.004, 723992.717, 4239681.595)
    pos = Position(*xyz)
    assert pos.xyz == pytest.approx(xyz, abs=1e-6)
    pos2 = Position(*pos.xyz,
                    s=Position.CoordinateSystem['cartesian'])
    assert pos2.xyz == pytest.approx(xyz, abs=1e-9)


def test_geocentric_constructor():
    lat, lon, r = 45.0, 100.0, 6500e3
    pos = Position(lat, lon, r, Position.CoordinateSystem['geocentric'])
    assert pos.geocentricLatitude == pytest.approx(lat, abs=1e-9)
    assert pos.longitude == pytest.approx(lon, abs=1e-9)
    assert pos.radius == pytest.approx(r, abs=1e-6)


def test_spherical_constructor():
    # point on the +Z axis
    pos = Position(0.0, 0.0, 6400e3, Position.CoordinateSystem['spherical'])
    assert pos.xyz == pytest.approx((0.0, 0.0, 6400e3), abs=1e-6)
    # zenith angle 90 at phi=90 lands on +Y
    pos = Position(90.0, 90.0, 6400e3, Position.CoordinateSystem['spherical'])
    assert pos.xyz == pytest.approx((0.0, 6400e3, 0.0), abs=1e-6)


def test_llh_matches_geo_xyz2geodetic():
    from pyrsss.gnss.geo import xyz2geodetic
    xyz = (4696986.004, 723992.717, 4239681.595)
    lat, lon, alt = Position(*xyz).llh
    lat2, lon2, alt2 = xyz2geodetic(*xyz)
    assert lat == pytest.approx(lat2, abs=1e-9)
    assert _wrap180(lon) == pytest.approx(_wrap180(lon2), abs=1e-9)
    assert alt == pytest.approx(alt2, abs=1e-6)


def test_geodetic_vs_geocentric_latitude():
    pos = Position(45.0, 20.0, 0.0, Position.CoordinateSystem['geodetic'])
    # geocentric latitude is smaller in magnitude away from the equator
    assert 0 < pos.geocentricLatitude < 45.0


def test_point_zenith_geometry():
    stn = Position(42.619, 288.51, 0.146,
                   Position.CoordinateSystem['geodetic'])
    p = point(stn, 0.0, 90.0, 1.0)  # straight up, 1 km
    assert p.radius == pytest.approx(stn.radius + 1e3, abs=1e-3)


def test_ipp_zenith():
    stn = Position(42.619, 288.51, 146.0,
                   Position.CoordinateSystem['geodetic'])
    ipp = ipp_from_azel(stn, 180.0, 90.0, ht=450)
    lat, lon, alt = ipp.llh
    # point() at el=90 follows the geocentric radial (Madrigal POINT
    # semantics), which differs from the geodetic normal by ~0.013 deg
    # at this latitude
    assert lat == pytest.approx(42.619, abs=0.05)
    assert _wrap180(lon) == pytest.approx(_wrap180(288.51), abs=0.05)
    assert alt == pytest.approx(450e3, rel=1e-2)


def test_unknown_coordinate_system():
    with pytest.raises(ValueError):
        Position(1.0, 2.0, 3.0, s=99)
