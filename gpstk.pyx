from libcpp.string cimport string

from enum import Enum


cdef extern from 'Position.hpp' namespace 'gpstk::Position':
    cdef enum CoordinateSystem:
        Unknown,
        Geodetic,
        Geocentric,
        Cartesian,
        Spherical


cdef extern from 'Position.hpp' namespace 'gpstk':
    cdef cppclass Position:
        Position(const double &,
                 const double &,
                 const double &) except +

        Position(const double &,
                 const double &,
                 const double &,
                 CoordinateSystem) except +

        Position(const double[3],
                 CoordinateSystem) except +

        string getSystemName()
        string asString() const


class PyCoordinateSystem(Enum):
    Geodetic = 1
    Geocentric = 2
    Cartesian = 3
    Spherical = 4


cdef class PyPosition:
    cdef Position *thisptr

    def __cinit__(self,
                  double a,
                  double b,
                  double c,
                  CoordinateSystem s=Cartesian):
        self.thisptr = new Position(a, b, c, s)

    def __repr__(self):
        return self.thisptr.asString() + ' (' + self.thisptr.getSystemName() + ')'
