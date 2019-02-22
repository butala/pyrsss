from datetime import datetime
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys
import logging
import os
import math

import numpy as NP
from netCDF4 import Dataset, date2num

from repository import get_root, TEMPLATE_MAP, PATH_MAP
from iaga2002 import parse
from util import mag_dec


logger = logging.getLogger('intermagnet_level0')


YEAR_START = 1991
"""
Default year to start processing.
"""

YEAR_END = 2013
"""
Default year to complete processing.
"""


def check_consistency(header1,
                      header2):
    """
    Return true if all critical fields of *header1* equal those of
    *header2*.
    """
    return (header1.Station_Name == header2.Station_Name and
            header1.IAGA_CODE == header2.IAGA_CODE and
            header1.Geodetic_Latitude == header2.Geodetic_Latitude and
            header1.Geodetic_Longitude == header2.Geodetic_Longitude and
            header1.Elevation == header2.Elevation and
            header1.Reported == header2.Reported and
            header1.Sensor_Orientation == header2.Sensor_Orientation and
            header1.Digital_Sampling == header2.Digital_Sampling and
            header1.Data_Interval_Type == header2.Data_Interval_Type and
            header1.Data_Type == header2.Data_Type)


def convert_HDZ_to_XYZ(H, D, Z):
    """
    Convert *H* (magnitude tangential to Earth's surface [nT]), *D*
    (declination [arcmin * 10]), and *Z* (downward magnitude [nT])) to
    *X* (magnitude of the geographic north pole component [nT]), *Y*
    (magnitude of the east component [nT]), and *Z* [nT].
    """
    D_rad = math.radians(arcmin_to_deg(D / 10))
    return (H * math.cos(D_rad),
            H * math.sin(D_rad),
            Z)


def arcmin_to_deg(x):
    """
    Convert *x* from [arcmin] to [deg].
    """
    return x / 60


def deg_to_arcmin(x):
    """
    Convert *x* from [deg] to [arcmin].
    """
    return x * 60


def convert_XYZ_to_HDZ(X, Y, Z):
    """
    Inverse of :func:`convert_heZ_to_XYZ`.
    """
    return (math.hypot(X, Y),
            deg_to_arcmin(math.degrees(math.atan2(Y, X))) * 10,
            Z)


NAME_MAP = {'X': 'magnitude of the geographic north pole component',
            'Y': 'magnitude of the east component',
            'Z': 'magnitude of the downward component',
            'F': 'magnetic field vector magnitude',
            'H': 'magnitude of the field tangential to the Earth''s surface',
            'D': 'declination (clockwise angle from the vector pointing to the geographic north pole to the magnetic field vector)'}
"""
Association between symbols and physical names.
"""


UNITS_MAP = {'X': 'nT',
             'Y': 'nT',
             'Z': 'nT',
             'F': 'nT',
             'H': 'nT',
             'D': 'arcmin'}
"""
Association between symbols and physical units.
"""


def setup_netcdf_root(root,
                      header,
                      N):
    """
    Add variables, dimensions (all of size *N), and attributes from
    *config* to :class:`Dataset` *root* and return *root*.
    """
    # add global attributes
    root.source = header.Source_of_Data
    root.station = header.Station_Name
    root.code = header.IAGA_CODE
    root.latitude = header.Geodetic_Latitude
    root.longitude = header.Geodetic_Longitude
    root.elevation = header.Elevation
    root.reported = header.Reported
    root.sensor_orientation = header.Sensor_Orientation
    root.data_interval_type = header.Data_Interval_Type
    root.data_type = header.Data_Type
    root.creation_time_utc = str(datetime.utcnow())
    # add dimensions
    time_dim = root.createDimension('time', N)
    codes = ['X', 'Y', 'Z', 'F', 'H', 'D']
    for code in codes:
        dim = root.createDimension(code, N)
    # add variables
    time = root.createVariable('time', 'f8', ('time',))
    time.standard_name = 'time'
    time.units = 'seconds since 2000-01-01 12:00:00.0'
    time.calendar = 'gregorian'
    for code in codes:
        var = root.createVariable(code,
                                  'f8',
                                  (code,),
                                  zlib=True)
        var.standard_name = NAME_MAP[code]
        var.units = UNITS_MAP[code]
    return root


def fill_missing(data):
    """
    Replace 99999 and 88888 entries with nan.
    """
    return [x if x not in [99999, 88888] else float('nan') for x in data]


def process_iaga2002(output_nc_fname,
                     input_iaga2002_fnames):
    """
    Convert IAGA2002 records in the list *input_iaga2002_fanmes* and
    store in the netCDF format record *output_nc_fname*. Return
    *output_nv_fname*.
    """
    last_header = None
    dt = []
    X = []
    Y = []
    Z = []
    H = []
    D = []
    F = []
    # gather information
    for iaga2002_fname in input_iaga2002_fnames:
        logger.info('processing {}'.format(iaga2002_fname))
        header, data_map = parse(iaga2002_fname)
        if last_header is not None:
            if not check_consistency(header,
                                     last_header):
                raise RuntimeError('header inconsistency detected')
        for dt_i, data_i in data_map.iteritems():
            # accumulate data in lists
            dt.append(dt_i)
            data_i = fill_missing(data_i)
            if header.Reported == 'HDZF':
                H.append(data_i[0])
                D.append(data_i[1])
                Z.append(data_i[2])
                F.append(data_i[3])
                # convert to XYZ
                X_i, Y_i, _ = convert_HDZ_to_XYZ(H[-1], D[-1], Z[-1])
                X.append(X_i)
                Y.append(Y_i)
            elif header.Reported in ['XYZF', 'XYZG']:
                X.append(data_i[0])
                Y.append(data_i[1])
                Z.append(data_i[2])
                F.append(data_i[3])
                # convert to heZ
                H_i, D_i, _ = convert_XYZ_to_HDZ(X[-1], Y[-1], Z[-1])
                H.append(H_i)
                D.append(D_i)
            else:
                raise RuntimeError('unknown Reported type {}'.format(header.Reported))
        last_header = header
    # write information to netCDF record
    root = Dataset(output_nc_fname, 'w')
    # add header, dimensions, and variables
    setup_netcdf_root(root,
                      header,
                      len(dt))
    # write information to file
    time = [date2num(dt_i, units=root['time'].units, calendar=root['time'].calendar) for dt_i in dt]
    root['time'][:] = time
    root['X'][:] = X
    root['Y'][:] = Y
    root['Z'][:] = Z
    root['F'][:] = F
    root['H'][:] = H
    root['D'][:] = D
    root.close()
    return output_nc_fname


def intermagnet_level0(root, years, stations='all'):
    """
    Process raw INTERMAGNET data stored at path *root* to level 0 for
    each year given in the list *years*. Process those sites listed in
    *stations* (`'all'` will process all stations). Return the list of
    generated level 0 files.
    """
    output_fnames = []
    definitive_path = os.path.join(PATH_MAP['definitive'].format(root=root))
    for year in years:
        logger.info('processing {}'.format(year))
        yyyy_path = os.path.join(definitive_path, str(year))
        for stn in sorted(os.listdir(yyyy_path)):
            if stations == 'all' or stn in stations:
                logger.info('processing {}'.format(stn))
                stn_path = os.path.join(yyyy_path, stn, 'IAGA2002')
                output_fnames.append(TEMPLATE_MAP[0].format(path=PATH_MAP[0].format(root=root),
                                                            year=year,
                                                            stn=stn))
                logger.info('writing to {}'.format(output_fnames[-1]))
                if not os.path.isdir(os.path.dirname(output_fnames[-1])):
                    os.makedirs(os.path.dirname(output_fnames[-1]))
                iaga2002_fnames = [os.path.join(stn_path, x) for x in sorted(os.listdir(stn_path))]
                process_iaga2002(output_fnames[-1], iaga2002_fnames)
    return output_fnames


def main(argv=None):
    if argv is None:
        argv = sys.argv

    root = get_root()

    parser = ArgumentParser('Process raw INTERMAGNET data to level 0.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('stations',
                        nargs='*',
                        type=str,
                        default='all')
    parser.add_argument('--year-range',
                        '-y',
                        nargs=2,
                        default=[YEAR_START,
                                 YEAR_END],
                        help='years to process (inclusive)')
    parser.add_argument('--root',
                        type=str,
                        default=root,
                        help='root path to data archive')
    args = parser.parse_args(argv[1:])

    intermagnet_level0(root,
                       range(args.year_range[0],
                             args.year_range[1] + 1),
                       stations=args.stations)



if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
