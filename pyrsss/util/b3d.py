import struct
from datetime import datetime, timedelta
from collections import namedtuple

import numpy as np



"""
Reference epoch (Adam cites IEEE Std. C37.118.2-2011).
"""
EPOCH = datetime(1970, 1, 1)


"""
Container to store B3D file format header information.
"""
class MetaData(namedtuple('Header',
                          'version '
                          'comments '
                          'channels '
                          'lon_0 '
                          'lon_step '
                          'lon_points '
                          'lat_0 '
                          'lat_step '
                          'lat_points '
                          'time_0 '
                          'time_step '
                          'time_points')):
    def __new__(cls, *args, **kwds):
        if args[-2] == 0:
            try:
                t = kwds.pop('t')
            except KeyError:
                raise ValueError('when time_step=0, the keyword argument t must be provided containing an array of times')
        self = super(MetaData, cls).__new__(cls, *args, **kwds)
        self.lon = self.lon_0 + np.arange(self.lon_points) * self.lon_step
        self.lat = self.lat_0 + np.arange(self.lat_points) * self.lat_step
        if self.time_step == 0:
            # read t as given by the file (providing support for
            # unevenly spaced time epochs)
            self.t = []
            for i in range(self.time_points):
                self.t.append(read_uint(fid))
        else:
            # build t from time_0, time_step, and time_points here
            self.t = self.time_0 + np.arange(self.time_points) * self.time_step * 1e-3
        self.dt = [EPOCH + timedelta(seconds=t_i) for t_i in self.t]
        return self


def read_uint(fid):
    """
    Read one 4-byte unsigned int number from *fid* and return the
    value.
    """
    return struct.unpack('<I', fid.read(4))[0]


def read_float(fid):
    """
    Read one 4-byte floating point number from *fid* and return the
    value.
    """
    return struct.unpack('<f', fid.read(4))[0]


def read_string(fid):
    """
    Read a null terminated string from *fid* and return the value.
    """
    s = ''
    while True:
        c_i = struct.unpack('<c', fid.read(1))[0]
        if c_i == '\0':
            break
        s += c_i
    return s


def write_uint(fid, i):
    """
    Write *i* as a one 4-byte unsigned record to *fid*.
    """
    fid.write(struct.pack('<I', i))
    return fid


def write_float(fid, f):
    """
    Write *f* as a one 4-byte floating point number to *fid*.
    """
    fid.write(struct.pack('<f', f))
    return fid


def write_string(fid, s):
    """
    Write *s* as a null terminated string to *fid*.
    """
    fid.write(s + '\0')
    return fid


def read_b3d(b3d_fname):
    """
    Read the B3D data cube formatted E-field data found in
    *b3d_fname*. Return the tuple Ex, Ey, and metadata. Ex and Ey are
    lists where the ith element is a 2-D array (number of latitude
    rows and number of longitude columns) for the electric field at
    the ith time index. The header information are stored in a
    :class:`MetaData` object.
    """
    with open(b3d_fname) as fid:
        # read header
        key = read_uint(fid)
        assert key == 34280
        version = read_uint(fid)
        assert version == 1
        N_meta_strings = read_uint(fid)
        assert N_meta_strings >= 0
        meta_strings = []
        for i in range(N_meta_strings):
            meta_strings.append(read_string(fid))
        channels = read_uint(fid)
        assert channels == 2
        lon_0 = read_float(fid)
        lon_step = read_float(fid)
        lon_points = read_uint(fid)
        lat_0 = read_float(fid)
        lat_step = read_float(fid)
        lat_points = read_uint(fid)
        time_0 = read_uint(fid)
        time_step = read_uint(fid)
        time_points = read_uint(fid)
        metadata = MetaData(version,
                            meta_strings,
                            channels,
                            lon_0,
                            lon_step,
                            lon_points,
                            lat_0,
                            lat_step,
                            lat_points,
                            time_0,
                            time_step,
                            time_points)
        # read data
        Ex = []
        Ey = []
        for dt_i in metadata.dt:
            Ex_i = []
            Ey_i = []
            for lat_j in metadata.lat:
                for lon_k in metadata.lon:
                    Ex_i.append(read_float(fid))
                    Ey_i.append(read_float(fid))
            Ex_i = np.array(Ex_i)
            Ey_i = np.array(Ey_i)
            Ex_i.shape = lat_points, lon_points
            Ey_i.shape = lat_points, lon_points
            Ex.append(Ex_i)
            Ey.append(Ey_i)
        # make sure we are at the end of file
        assert fid.read() == ''
    return Ex, Ey, metadata


def write_b3d(b3d_fname, Ex, Ey, metadata, key=34280, version=1, channels=2):
    """
    Write the electric field data *Ex* and *Ey* with header
    information *metadata* in B3D format to *b3d_fname*. *Ex* and *Ey*
    are lists, one element per time index, and each element is a 2-D
    array record (number of latitude rows and number of longitude
    columns). The *metadata* record is expected to be a
    :class:`MetaData` object. The keyword arguments are: *key* is the
    file magic number, *version* is the format version (only version 1
    is currently supported), and *channels* denotes the number of
    electric field vector component arrays (only 2 are currently
    supported).
    """
    # sanity checks
    assert metadata.version == 1
    assert len(Ex) == len(Ey)
    assert Ex[0].shape == Ey[0].shape
    rows, cols = Ex[0].shape
    assert len(metadata.lat) == rows
    assert len(metadata.lon) == cols
    assert metadata.channels == 2
    with open(b3d_fname, 'w') as fid:
        # write header
        write_uint(fid, key)
        write_uint(fid, version)
        write_uint(fid, len(metadata.comments))
        for comment in metadata.comments:
            write_string(fid, comment)
        write_uint(fid, channels)
        write_float(fid, metadata.lon_0)
        write_float(fid, metadata.lon_step)
        write_uint(fid, metadata.lon_points)
        write_float(fid, metadata.lat_0)
        write_float(fid, metadata.lat_step)
        write_uint(fid, metadata.lat_points)
        write_uint(fid, metadata.time_0)
        write_uint(fid, metadata.time_step)
        write_uint(fid, metadata.time_points)
        if metadata.time_step == 0:
            for t_i in metadata.t:
                write_uint(fid, int((t_i - metadata.time_0) * 1e3))
        # write data
        for i in range(len(Ex)):
            Ex_i = Ex[i]
            Ey_i = Ey[i]
            for Ex_i_jk, Ey_i_jk in zip(Ex_i.flat, Ey_i.flat):
                write_float(fid, Ex_i_jk)
                write_float(fid, Ey_i_jk)
    return b3d_fname
