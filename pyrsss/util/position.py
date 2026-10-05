"""
Pure-Python WGS84 position and geometry utilities. Drop-in replacement
for the retired Cython ``pyrsss.gnsstk`` extension (which wrapped
``gnsstk::Position``). Consumers only used ECEF/geodetic/geocentric
conversions and the Madrigal-style :func:`point` geometry, all of which
are closed-form on the WGS84 ellipsoid.
"""
import math

import pyproj


_ECEF = pyproj.CRS(proj='geocent', ellps='WGS84', datum='WGS84')
_LLA = pyproj.CRS(proj='latlong', ellps='WGS84', datum='WGS84')
_TO_LLA = pyproj.Transformer.from_crs(_ECEF, _LLA, always_xy=True)
_TO_ECEF = pyproj.Transformer.from_crs(_LLA, _ECEF, always_xy=True)


class Position(object):
    """
    A point relative to the WGS84 ellipsoid. Constructed in one of four
    coordinate systems (the default is Cartesian):

    - ``'cartesian'``:  Earth-centered, Earth-fixed (x, y, z) [m]
    - ``'geodetic'``:   (geodetic lat [deg], lon [deg], height [m])
    - ``'geocentric'``: (geocentric lat [deg], lon [deg], radius [m])
    - ``'spherical'``:  (zenith theta [deg], azimuth phi [deg], radius [m])
    """

    CoordinateSystem = {'geodetic':   1,
                        'geocentric': 2,
                        'cartesian':  3,
                        'spherical':  4}

    def __init__(self, a, b, c, s=3):
        if s == self.CoordinateSystem['cartesian']:
            self._x, self._y, self._z = float(a), float(b), float(c)
        elif s == self.CoordinateSystem['geodetic']:
            self._x, self._y, self._z = _TO_ECEF.transform(float(b),
                                                           float(a),
                                                           float(c))
        elif s == self.CoordinateSystem['geocentric']:
            lat, lon, r = math.radians(a), math.radians(b), c
            self._x = r * math.cos(lat) * math.cos(lon)
            self._y = r * math.cos(lat) * math.sin(lon)
            self._z = r * math.sin(lat)
        elif s == self.CoordinateSystem['spherical']:
            theta, phi, r = math.radians(a), math.radians(b), c
            self._x = r * math.sin(theta) * math.cos(phi)
            self._y = r * math.sin(theta) * math.sin(phi)
            self._z = r * math.cos(theta)
        else:
            raise ValueError('Unknown coordinate system {}'.format(s))

    def __repr__(self):
        return 'Position(x={:.3f}, y={:.3f}, z={:.3f}) [m ECEF]'.format(self._x,
                                                                       self._y,
                                                                       self._z)

    @property
    def x(self):
        """Return the ECEF X coordinate [m]."""
        return self._x

    @property
    def y(self):
        """Return the ECEF Y coordinate [m]."""
        return self._y

    @property
    def z(self):
        """Return the ECEF Z coordinate [m]."""
        return self._z

    @property
    def xyz(self):
        """Return the tuple of ECEF coordinates (x, y, z) all in [m]."""
        return (self._x, self._y, self._z)

    @property
    def radius(self):
        """Return the distance from the center of the Earth [m]."""
        return math.sqrt(self._x**2 + self._y**2 + self._z**2)

    @property
    def geocentricLatitude(self):
        """Return the geocentric latitude [deg N]."""
        return math.degrees(math.asin(self._z / self.radius))

    @property
    def longitude(self):
        """Return the longitude [deg E] in [0, 360)."""
        return math.degrees(math.atan2(self._y, self._x)) % 360

    @property
    def geodeticLatitude(self):
        """Return the geodetic latitude [deg N]."""
        return self.llh[0]

    @property
    def height(self):
        """Return the height above the ellipsoid [m]."""
        return self.llh[2]

    @property
    def llh(self):
        """
        Return the tuple of geodetic coordinates (lat, lon, height) where
        lat and lon are in [deg] (lon in [0, 360)) and height is in [m].
        """
        lon, lat, alt = _TO_LLA.transform(self._x, self._y, self._z)
        return (lat, lon % 360, alt)


PyPosition = Position
"""Deprecated alias for :class:`Position`."""


def point(stn_point, target_az, target_el, target_range):
    """
    Return the :class:`Position` given by the target azimuth *target_az*
    [deg], elevation *target_el* [deg], and range *target_range* [km]
    relative to the station :class:`Position` *stn_point*.

    The function replicates the Fortran POINT subroutine implemented in
    Madrigal (http://madrigal.haystack.edu/madrigal/madDownload.html).

    Example (Millstone Hill IS Radar at latitude=42.619,
    longitude=288.51, altitude=0.146; azimuth=100, elevation=45, 55,
    and 65, range=1000 [km]):

    41.361338,-64.012651,742.038442
    41.645771,-65.488196,841.794521
    41.932738,-67.098375,918.473216
    """
    sr = stn_point.radius / 1e3
    slat = stn_point.geocentricLatitude
    slon = stn_point.longitude
    theta = math.radians(90 - target_el)
    phi = math.radians(180 - target_az)
    rt = target_range * math.sin(theta) * math.cos(phi)
    rp = target_range * math.sin(theta) * math.sin(phi)
    rr = target_range * math.cos(theta)
    theta = math.radians(slat)
    phi = math.radians(slon)
    ct, st, cp, sp = math.cos(theta), math.sin(theta), math.cos(phi), math.sin(phi)
    fx = ct * cp * rr + st * cp * rt - sp * rp
    fy = ct * sp * rr + st * sp * rt + cp * rp
    fz = st * rr - ct * rt
    stn_theta = math.radians(90 - slat)
    stn_phi = math.radians(slon)
    sx = sr * math.sin(stn_theta) * math.cos(stn_phi)
    sy = sr * math.sin(stn_theta) * math.sin(stn_phi)
    sz = sr * math.cos(stn_theta)
    return Position((sx + fx) * 1e3,
                    (sy + fy) * 1e3,
                    (sz + fz) * 1e3)
