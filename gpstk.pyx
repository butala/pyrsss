from libcpp cimport bool
from libcpp.string cimport string
from libcpp.cast cimport dynamic_cast
from cython.operator cimport dereference as deref

from enum import Enum


cdef extern from 'Triple.hpp' namespace 'gpstk':
    cdef cppclass Triple:
        double operator[](const size_t) const


cdef extern from 'Position.hpp' namespace 'gpstk::Position':
    cdef enum CoordinateSystem:
        Unknown = 0,
        Geodetic,
        Geocentric,
        Cartesian,
        Spherical


cdef extern from 'Position.hpp' namespace 'gpstk':
    cdef cppclass Position(Triple):
        Position() except +

        Position(const double &,
                 const double &,
                 const double &) except +

        Position(const double &,
                 const double &,
                 const double &,
                 CoordinateSystem) except +

        Position(const double[3],
                 CoordinateSystem) except +

        string getSystemName() except +
        string asString() except +

        CoordinateSystem getCoordinateSystem() except +

        Position asGeodetic() except +
        Position asECEF() except +

        double X() except +
        double Y() except +
        double Z() except +

        double geodeticLatitude() except +
        double geocentricLatitude() except +

        double longitude() except +

        double radius() except +
        double height() except +

        double theta() except +
        double phi() except +

        double elevation(const Position &) except +
        double azimuth(const Position &) except +

        double elevationGeodetic(const Position &) except +
        double azimuthGeodetic(const Position &) except +

        Position getIonosphericPiercePoint(const double,
                                           const double,
                                           const double) except +

        bool operator==(const Position &) except +

        bool operator!=(const Position &) except +


    cdef double range(const Position &,
                      const Position &) except +


    cdef Position operator-(const Position&,
                            const Position&) except +

    cdef Position operator+(const Position&,
                            const Position&) except +

    cdef Position operator*(const double &,
                            const Position &) except +

    cdef Position operator*(const Position &,
                            const double &) except +

    cdef Position operator*(const int &,
                            const Position &) except +

    cdef Position operator*(const Position &,
                            const int &) except +


cdef class PyPosition:
    # TODO: Add tolerance get / set

    cdef Position *thisptr

    CoordinateSystem = {
        'geodetic'  : 1,
        'geocentric': 2,
        'cartesian' : 3,
        'spherical' : 4
    }

    def __cinit__(self,
                  double a,
                  double b,
                  double c,
                  int s = 3):
        """
        Note that the reference ellipsoid is WGS84 (overridable, but not
        yet implemented.)
        """
        if s == PyPosition.CoordinateSystem['geodetic']:
            self.thisptr = new Position(a, b, c, Geodetic)
        elif s == PyPosition.CoordinateSystem['geocentric']:
            self.thisptr = new Position(a, b, c, Geocentric)
        elif s == PyPosition.CoordinateSystem['cartesian']:
            self.thisptr = new Position(a, b, c, Cartesian)
        elif s == PyPosition.CoordinateSystem['spherical']:
            self.thisptr = new Position(a, b, c, Spherical)
        else:
            raise ValueError('Unknown coordinate system {}'.format(s))

    def __repr__(self):
        return self.thisptr.asString() + ' (' + self.thisptr.getSystemName() + ')'

    @staticmethod
    cdef PyPosition toPyPosition(Position p):
       return PyPosition(p[0], p[1], p[2], p.getCoordinateSystem())

    def __richcmp__(PyPosition self, PyPosition other, int op):
        """
        See http://docs.cython.org/src/userguide/special_methods.html#rich-comparisons

        Return true if the distance between this position and *other*
        is less than the tolerance.
        """
        cdef Position _this = deref(self.thisptr)
        cdef Position _other = deref(other.thisptr)
        if op == 2:  # ==
            return _this == _other
        elif op == 3:  # !=
            return _this != _other
        else:
            raise NotImplementedError('only == and != comparisons are currently supported')

    def asGeodetic(self):
        """Transform to geodetic coordinate system."""
        cdef Position p = deref(self.thisptr).asGeodetic()
        return PyPosition.toPyPosition(p)

    def asECEF(self):
        """Transform to Cartesian coordinate system."""
        cdef Position p = deref(self.thisptr).asECEF()
        return PyPosition.toPyPosition(p)

    @property
    def x(self):
        """Return the ECEF X coordinate [m]."""
        return deref(self.thisptr).X()

    @property
    def y(self):
        """Return the ECEF Y coordinate [m]."""
        return deref(self.thisptr).Y()

    @property
    def z(self):
        """Return the ECEF Z coordinate [m]."""
        return deref(self.thisptr).Z()

    @property
    def geodeticLatitude(self):
        """Return the geodetic latitude [deg N]."""
        return deref(self.thisptr).geodeticLatitude()

    @property
    def geocentricLatitude(self):
        """Return the geocentric latitude [deg N]."""
        return deref(self.thisptr).geocentricLatitude()

    @property
    def longitude(self):
        """Return the longitude [deg E]."""
        return deref(self.thisptr).longitude()

    @property
    def radius(self):
        """Return the distance from the center of the Earth [m]."""
        return deref(self.thisptr).radius()

    @property
    def height(self):
        """Return the height above the ellipsoid [m]."""
        return deref(self.thisptr).height()

    @property
    def theta(self):
        """Return spherical zenith [deg]."""
        return deref(self.thisptr).theta()

    @property
    def phi(self):
        """Return spherical azimuth [deg]."""
        return deref(self.thisptr).phi()

    def elevation(self, PyPosition target):
        """Return the elevation [deg] from this position to *target*."""
        cdef Position _target = deref(target.thisptr)
        return deref(self.thisptr).elevation(_target)

    def elevationGeodetic(self, PyPosition target):
        """
        Return the elevation [deg] from this position to *target* using a
        geodetic system.
        """
        cdef Position _target = deref(target.thisptr)
        return deref(self.thisptr).elevationGeodetic(_target)

    def azimuth(self, PyPosition target):
        """Return the azimuth [deg] from this position to *target*."""
        cdef Position _target = deref(target.thisptr)
        return deref(self.thisptr).azimuth(_target)

    def azimuthGeodetic(self, PyPosition target):
        """
        Return the azimuth [deg] from this position to *target* using a
        geodetic system.
        """
        cdef Position _target = deref(target.thisptr)
        return deref(self.thisptr).azimuthGeodetic(_target)

    def getIPP(self, double elevation, double azimuth, double shell_height=450 * 1e3):
        """
        Return the position at which a signal received at this location at
        *elevation* [deg] and *azimuth* [deg] cross a thin shell
        ionosphere at height *shell_height* [m].
        """
        cdef Position ipp = deref(self.thisptr).getIonosphericPiercePoint(elevation,
                                                                          azimuth,
                                                                          shell_height)
        return PyPosition.toPyPosition(ipp)

    def distance(self, PyPosition target):
        """
        Return the distance [m] between this position and *target*.
        """
        cdef Position _this = deref(self.thisptr)
        cdef Position _target = deref(target.thisptr)
        return range(_this, _target)
