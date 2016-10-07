from __future__ import division
import logging
import sys
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import datetime, timedelta
from collections import OrderedDict, namedtuple, defaultdict
from abc import ABCMeta, abstractmethod, abstractproperty


logger = logging.getLogger('pyrsss.mag.iaga2002')


def fname2date(iaga2002_fname):
    """
    ???
    """
    return datetime.strptime(os.path.basename(iaga2002_fname)[3:11],
                             '%Y%m%d')


class Header(namedtuple('Header',
                        ['Format',
                         'Source_of_Data',
                         'Station_Name',
                         'IAGA_CODE',
                         'Geodetic_Latitude',
                         'Geodetic_Longitude',
                         'Elevation',
                         'Reported',
                         'Sensor_Orientation',
                         'Digital_Sampling',
                         'Data_Interval_Type',
                         'Data_Type',
                         'Comment'])):
    pass


HEADER_TYPES = defaultdict(lambda: str,
                           [('Geodetic_Latitude', float),
                            ('Geodetic_Longitude', float),
                            ('Elevation', float)])


def convert_float(s):
    """
    Convert the string data field *s* to a float. If the value is
    99999 (missing data) or 88888 (not observed), return not a number.
    """
    f = float(s)
    if int(f) in [99999, 88888]:
        return float('nan')
    return f


class IAGARecord(object):
    __metaclass__ = ABCMeta

    @abstractproperty
    def x(self):
        pass

    @abstractproperty
    def y(self):
        pass

    @abstractproperty
    def z(self):
        pass

    @abstractproperty
    def f(self):
        pass

    @abstractproperty
    def H(self):
        pass

    @abstractproperty
    def D(self):
        pass


class XYZRecord(IAGARecord):
    def __init__(self, x, y, z, f):
        self._x = x
        self._y = y
        self._z = z
        self._f = f

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def z(self):
        return self._z

    @property
    def f(self):
        return self._f

    @property
    def H(self):
        return math.hypot(self._x, self._y)

    @property
    def D(self):
        return math.degrees(math.atan2(self._y, self._x))

    def __repr__(self):
        return 'XYZRecord(x={} y={} z={} f={})'.format(self._x,
                                                       self._y,
                                                       self._z,
                                                       self._f)


def HDZRecord(IAGARecord):
    def __init__(self, H, D_arcmin, z, f):
        self._H = H
        self._D_deg = D_arcmin / 60
        self._D_rad = math.radians(self._D_deg)
        self._z = z
        self._f = f

    @property
    def x(self):
        return self._H * math.cos(self._D_rad)

    @property
    def y(self):
        return self._H * math.sin(self._D_rad)

    @property
    def z(self):
        return self._z

    @property
    def f(self):
        return self._f

    @property
    def H(self):
        return self._H

    @property
    def D(self):
        return self._D_deg

    def __repr__(self):
        return 'HDZRecord(H={} D={} [deg] z={} f={})'.format(self._H,
                                                             self._D_deg,
                                                             self._z,
                                                             self._f)


def record_factory(reported, data):
    """
    ???

    What exists in the INTERMAGNET data record:
    Reported               HDZ                                          |
    Reported               HDZF                                         |
    Reported               XYZ                                          |
    Reported               xyzf                                         |
    Reported               XYZF                                         |
    Reported               XYZG                                         |
    """
    if reported.startswith('HDZ'):
        return HDZRecord(*data)
    elif reported.startswith('XYZ'):
        return XYZRecord(*data)
    else:
        raise NotImplementedError('unknown record type {}'.format(reported))


def parse(fname, strict=True):
    """
    Parser the IAGA2002 format file *fname* and return a tuple with a
    :class:`Header` and mapping of date/times to measured values. If
    *strict*, fail when an nonconforming entry is encountered. If
    *strict* is not set, attempt to carry on when simple parse errors
    are encountered.
    """
    with open(fname) as fid:
        # parse header
        header_map = {}
        comment_lines = []
        for line in fid:
            if line[69] != '|' and line.rstrip()[-1] != '|':
                raise RuntimeError('malformed header line in {} ({}) --- expected | but found "{}"'.format(fname,
                                                                                                           line,
                                                                                                           line[69]))
            elif line[1] == '#':
                comment_lines += [line[3:69].rstrip()]
            elif not strict and line[0] == '#':
                comment_lines += [line[3:69].rstrip()]
            elif line.startswith('DATE'):
                break
            else:
                key = line[:24].strip().replace(' ', '_')
                try:
                    header_map[key] = HEADER_TYPES[key](line[24:69].strip())
                except ValueError, e:
                    logger.warning('could not parse header line {} --- skipping'.format(line))
                    header_map[key] = None
                    continue
        header_map['Comment'] = '\n'.join(comment_lines)
        try:
            header = Header(**header_map)
        except TypeError:
            if strict:
                logger.warning('unknown header record(s) found in {} --- setting header to None'.format(fname))
                header = None
            else:
                header = header_map
        # parse data header record
        fields = line[:69].split()
        if len(fields) != 7:
            raise RuntimeError('malformed data header record in {} ({})'.format(fname,
                                                                                line))
        data_map = OrderedDict()
        # parse data records
        for line in fid:
            try:
                dt = datetime.strptime(line[:23], '%Y-%m-%d %H:%M:%S.%f')
            except ValueError, e:
                if strict:
                    raise e
                else:
                    dt = datetime.strptime(line[:19], '%Y-%m-%d %H:%M:%S')
            data1 = convert_float(line[31:40])
            data2 = convert_float(line[41:50])
            data3 = convert_float(line[51:60])
            data4 = convert_float(line[61:70])
            data_map[dt] = record_factory(header_map['Reported'],
                                          [data1,
                                           data2,
                                           data3,
                                           data4])
    return header, data_map


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Dump a IAGA2002 format file.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('iaga2002_fname',
                        type=str,
                        help='IAGA2002 file name')
    args = parser.parse_args(argv[1:])

    header, data_map = parse(args.iaga2002_fname)

    for key, value in header._asdict().iteritems():
        print('{} = {}'.format(key.replace(' ', '-'), value))
    for dt, values in data_map.iteritems():
        print('{:%Y-%m-%d %H:%M:%S.%f}:  {}  {}  {}  {}'.format(dt,
                                                                values[1],
                                                                values[2],
                                                                values[3],
                                                                values[4]))

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
