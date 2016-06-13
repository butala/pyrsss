from libc.math cimport sin, cos, M_PI
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.cast cimport dynamic_cast
from libcpp.map cimport map as cppmap
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from cython.operator cimport dereference as deref, preincrement as inc

from collections import namedtuple

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

"""
enum CoordinateSystem
      {
         Unknown=0,  ///< unknown coordinate system
         Geodetic,   ///< geodetic latitude, longitude, and height above ellipsoid
         Geocentric, ///< geocentric (regular spherical coordinates)
         Cartesian,  ///< cartesian (Earth-centered, Earth-fixed)
         Spherical   ///< spherical coordinates (theta,phi,radius)
      };
"""


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


cdef class PyPosition(object):
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
        cdef Position _self = deref(self.thisptr)
        cdef Position _other = deref(other.thisptr)
        if op == 2:  # ==
            return _self == _other
        elif op == 3:  # !=
            return _self != _other
        else:
            raise NotImplementedError('only == and != comparisons are currently supported')

    def __add__(PyPosition self, PyPosition other):
        """
        Return the sum of this position and *other* (in Cartesian
        coordinates).
        """
        cdef Position _self = deref(self.thisptr)
        cdef Position _other = deref(other.thisptr)
        return PyPosition.toPyPosition(_self + _other)

    def __sub__(PyPosition self, PyPosition other):
        """
        Return the difference of this position and *other* (in Cartesian
        coordinates).
        """
        cdef Position _self = deref(self.thisptr)
        cdef Position _other = deref(other.thisptr)
        return PyPosition.toPyPosition(_self - _other)

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
    def xyz(self):
        """
        Return the tuple of EXEC coordinates (x, y, z) all in [m].
        """
        return (self.x, self.y, self.z)

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
        cdef Position _self = deref(self.thisptr)
        cdef Position _target = deref(target.thisptr)
        return range(_self, _target)


################################################################################


def point(PyPosition stn_point, double target_az, double target_el, double target_range):
   """
   Return the :class:`PyPosition` to the point given by the target
   azimuth *target_az* [deg], elevation *target_el* [deg], and range
   *target_range* [km] relative to the station :class:`PyPosition`
   *stn_point*.

   The function replicates the Fortran POINT subroutine implemented in
   Madrigal (http://madrigal.haystack.edu/madrigal/madDownload.html).
   """
   cdef sr = stn_point.radius / 1e3
   cdef slat = stn_point.geocentricLatitude
   cdef slon = stn_point.longitude
   cdef Position pos = Position(90 - target_el, 180 - target_az, target_range, Spherical)
   cdef double rt = pos.X()
   cdef double rp = pos.Y()
   cdef double rr = pos.Z()
   cdef Position pos_vctcnv = Position(90 - slat, slon, sr, Spherical)
   cdef double theta = slat * M_PI / 180
   cdef double phi = slon * M_PI / 180
   cdef double ct = cos(theta)
   cdef double st = sin(theta)
   cdef double cp = cos(phi)
   cdef double sp = sin(phi)
   cdef double fx = ct*cp*rr + st*cp*rt - sp*rp
   cdef double fy = ct*sp*rr + st*sp*rt + cp*rp
   cdef double fz = st*rr - ct*rt
   return PyPosition((pos_vctcnv.X() + fx) * 1e3,
                     (pos_vctcnv.Y() + fy) * 1e3,
                     (pos_vctcnv.Z() + fz) * 1e3)


################################################################################


# cdef extern from "<iostream>" namespace "std::ios":
#     cdef cppclass openmode:
#         pass


cdef extern from 'RinexObsID.hpp' namespace 'gpstk':
    cdef cppclass RinexObsID:
        string asString() const


ctypedef vector[RinexObsID] RinexObsVec


ctypedef cppmap[string, RinexObsVec] RinexObsMap


cdef extern from '<streambuf>' namespace 'std':
    cdef cppclass streambuf:
        streambuf()


cdef extern from '<iostream>' namespace 'std':
    cdef cppclass istream:
        istream(streambuf *)

    cdef cppclass iostream(istream):
        iostream(streambuf *)

    cdef cppclass fstream(iostream):
        fstream(const char *)


cdef extern from 'FFStream.hpp' namespace 'gpstk':
   cdef cppclass FFStream(fstream):
       pass


cdef extern from 'FFTextStream.hpp' namespace 'gpstk':
    cdef cppclass FFTextStream(FFStream):
        pass


cdef extern from 'Rinex3ObsStream.hpp' namespace 'gpstk':
    cdef cppclass Rinex3ObsStream(FFTextStream):
        Rinex3ObsStream()
        Rinex3ObsStream(const char *)
        Rinex3ObsStream(const string &)


cdef extern from 'FFData.hpp' namespace 'gpstk':
    cdef cppclass FFData:
        pass

    # note that adding 'except +' causes problems
    cdef istream & operator>>(istream &,
                              FFData &)


cdef extern from 'Rinex3ObsBase.hpp' namespace 'gpstk':
    cdef cppclass Rinex3ObsBase(FFData):
        pass


cdef extern from 'Rinex3ObsHeader.hpp' namespace 'gpstk':
    cdef cppclass Rinex3ObsHeader(Rinex3ObsBase):
        Rinex3ObsHeader()

        # TODO: Add other members

        double      version
        string      fileType
        string      fileSys
        string      date
        string      recType
        Triple      antennaPosition
        RinexObsMap mapObsTypes
        double      interval

        bool isValid() const


class RinexInfo(namedtuple('RinexInfo',
                           'version '
                           'file_type '
                           'file_sys '
                           'date '
                           'rec_type '
                           'position '
                           'obs_types '
                           'interval')):
    pass


cdef obsTypes2dict(RinexObsMap &obs_map):
    """
    ???
    """
    cdef cppmap[string, RinexObsVec].iterator it = obs_map.begin()
    obs_dict = {}
    while it != obs_map.end():
        obs_dict[deref(it).first] = [x.asString() for x in deref(it).second]
        inc(it)
    return obs_dict


def rinex_header_info(string fname):
    """
    """
    cdef Rinex3ObsHeader Rhead
    cdef Rinex3ObsStream *istrm_ptr = new Rinex3ObsStream(fname)
    try:
        deref(istrm_ptr) >> Rhead
        position = Rhead.antennaPosition
        return RinexInfo(Rhead.version,
                         Rhead.fileType,
                         Rhead.fileSys,
                         Rhead.date,
                         Rhead.recType,
                         [position[i] for i in [0, 1, 2]],
                         obsTypes2dict(Rhead.mapObsTypes),
                         Rhead.interval)
    finally:
        del istrm_ptr
