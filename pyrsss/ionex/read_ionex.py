"""
IONEX Parser & Interpolator v1.1
written by 

Matthew A. Grawe
University of Illinois at Urbana-Champaign
grawe2@illinois.edu
January 23, 2016

Developed in accordance with the format specifications described in
"IONEX: The IONosphere Map EXchange Format Version 1"
written by

Stefan Schaer, Werner Gurtner
Astronomical Institute, University of Berne, Switzerland

and 

Joachim Feltens
ESA/ESOC, Darmstadt, Germany

Published on February 25, 1998.
"""

import numpy as np
import datetime

def next_line(ionex_file):
    """
    next_line
    Function returns the next line in the file
    that is not a blank line, unless the line is
    '', which is a typical EOF marker.
    """

    done   = False
    while not done:
        line = ionex_file.readline()
        if line == '':
            return line
        elif line.strip():
            return line

def read_misc_line(ionex_file):
    line           = next_line(ionex_file)
    if line == '':
        label   = None
        content = None
    else:
        label   = line[60:80]
        content = line[0:60]

    return label, content

def read_lon_slice(ionex_file, n_lons):
    """
    read_lon_slice
    Reads in a single lon slice of data (for a particular epoch, lat, height)
    as a 1D numpy array filled with the data (TEC/RMS/HEIGHT)
    """
    # 16 ENTRIES MAX PER LINE, I5 FORMAT #
    line_count = int(np.ceil(n_lons/16.))
    data       = []

    for i in range(0, line_count):
        remaining = 16
        if i == line_count - 1:
            # LAST LINE: CHECK IF THE LINE IS FULL OR NOT #
            remaining = n_lons % 16
            if remaining == 0:
                remaining = 16

        line = next_line(ionex_file)
        pos  = 0
        for j in range(0, remaining):
            try:
                val = int(line[pos:pos+5])
                if val == 9999:
                    val = np.nan
                data.append(val)
            except Exception as e:
                print e
                raise
            pos = pos + 5

    return np.array(data)

def read_lonlat_block(ionex_file, n_lons, n_lats):
    """
    read_lonlat_block
    Reads in a set of n_lat lon slices and returns the data
    as a 2D numpy array filled with the data (TEC/RMS/HEIGHT)
    """
    data = np.zeros([n_lons, n_lats])
    for i in range(0, n_lats):
        label, content = read_misc_line(ionex_file)
        if 'LAT/LON1/LON2/DLON/H' in label:
           # CONTENT: 2X, 5F6.1, 28X
            try:
                # Not sure if this data is needed (seems redundant), but it is available here.
                lat  = float(content[2:8])
                lon1 = float(content[8:14])
                lon2 = float(content[14:20])
                dlon = float(content[20:26])
                h    = float(content[26:32])
            except Exception as e:
                print e
                raise

            # READ IN THIS PARTICULAR LON SLICE #
            lon_slice = read_lon_slice(ionex_file, n_lons)
            data[:,i] = lon_slice

        else:
            raise Exception("read_lonlat_block encountered an unexpected line:\n %s" % content)

    return data

def read_block_header(ionex_file, exponent):
    """
    read_block_header
    Reads the lines that precede each data block.
    Sometimes, a new exponent may be specified for
    a data block, so there needs to be a check here
    for that. It is unclear from the specification
    whether a newly defined exponent is intended to
    apply to data from subsequent maps. Here, the
    assumption is made that exponent redefinitions
    are local only to the map they are defined in.
    The function returns the epoch and exponent.
    If no epoch was present, None is returned in its
    place. If no exponent was present, the exponent
    defined in the header is returned or the default
    exponent (-1) is returned if one was not defined
    in the header.
    """
    epoch        = None
    new_exponent = exponent

    done = False
    while not done:
        before = ionex_file.tell()
        label, content = read_misc_line(ionex_file)
        if 'LAT/LON1/LON2/DLON/H' in label:
            # END OF BLOCK HEADER, BACK UP ONE LINE #
            # (this is because read_lonlat_block expects this line) #
            ionex_file.seek(before)
            return epoch, new_exponent
        elif 'EXPONENT' in label:
            try:
                new_exponent = int(content[0:6])
            except Exception as e:
                print e
                raise
        elif 'EPOCH OF CURRENT MAP':
            try:
                year    = int(content[0:6])
                month   = int(content[6:12])
                day     = int(content[12:18])
                hour    = int(content[18:24])
                minute  = int(content[24:30])
                sec     = int(content[30:36])

                epoch = datetime.datetime(year, month, day, \
                                                hour, minute, sec)

            except Exception as e:
                print e
                raise

def read_map(ionex_file, n_lons, n_lats, n_heights, exponent):
    """
    read_map
    Reads in an entire map (TEC/RMS/HEIGHT) over
    longitude, latitude, and height. The epoch and 
    data is returned as a datetime and 3D numpy array
    (respectively) in the form of a 3-tuple. The numpy
    array has shape (n_lats, n_lons, n_heights). If there
    was no epoch listed, None is returned in place of
    the epoch.
    """

    epoch, exponent = read_block_header(ionex_file, exponent)

    # FOR EACH HEIGHT, READ IN A BLOCK OF LAT, LON DATA #
    tec_block = np.zeros([n_lons, n_lats, n_heights])
    for k in range(0, n_heights):
        tec_block[:,:,k] = read_lonlat_block(ionex_file, n_lons, n_lats)*10**(exponent)

    return epoch, tec_block

def read_header(ionex_file):
    """
    Header lines are first read in as 2-tuples.
    The first component is a string of the first 60 columns (CONTENT)
    The second component is a string of the next 20 columns (HEADER LABEL)
    """
    # READ IN HEADER INFORMATION ALL AT ONCE #

    header_lines = []

    done = False
    while not done: 

        # CHECK FOR EPOCH LINE #
        header_label, header_content = read_misc_line(ionex_file)

        header_lines.append((header_label, header_content))

        if 'END OF HEADER' in header_label:
            done = True

    # DISTRIBUTE RELEVANT HEADER INFORMATION INTO A DICTIONARY #
    import datetime
    header_info = {}

    satellite_biases = {}
    station_biases   = {}

    satellite_biases['GPS']     = {}
    satellite_biases['GLONASS'] = {}

    station_biases['GPS']     = {}
    station_biases['GLONASS'] = {}

    # default exponent is -1
    header_info['exponent'] = -1

    for label, content in header_lines:

        if 'EPOCH OF FIRST MAP' in label:
            try:
                year    = int(content[0:6])
                month   = int(content[6:12])
                day     = int(content[12:18])
                hour    = int(content[18:24])
                minute  = int(content[24:30])
                sec     = int(content[30:36])

                header_info['start_epoch'] = datetime.datetime(year, month, day, \
                                                hour, minute, sec)
            except Exception as e:
                print e
                raise
        elif 'EPOCH OF LAST MAP' in label:
            try:
                year    = int(content[0:6])
                month   = int(content[6:12])
                day     = int(content[12:18])
                hour    = int(content[18:24])
                minute  = int(content[24:30])
                sec     = int(content[30:36])

                header_info['stop_epoch'] = datetime.datetime(year, month, day, \
                                                hour, minute, sec)
            except Exception as e:
                print e
                raise
        elif '# OF MAPS IN FILE' in label:
            try:
                header_info['map_count'] = int(content[0:6])
            except Exception as e:
                print e
                raise
        elif 'MAP DIMENSION' in label:
            try:
                map_dimension = int(content[0:6])

                if map_dimension not in [2,3]:
                    raise Exception("Map dimension is not in [2,3]")

                header_info['map_dimension'] = map_dimension

            except Exception as e:
                print e
                raise

        elif 'HGT1 / HGT2 / DHGT' in label:
            try:
                height1 = float(content[2:8])
                height2 = float(content[8:14])
                dh      = float(content[14:20])

                header_info['height1'] = height1
                header_info['height2'] = height2
                header_info['dh']      = dh
            except Exception as e:
                print e
                raise
        elif 'LAT1 / LAT2 / DLAT' in label:
            try:
                lat1 = float(content[2:8])
                lat2 = float(content[8:14])
                dlat = float(content[14:20])

                header_info['lat1'] = lat1
                header_info['lat2'] = lat2
                header_info['dlat'] = dlat
            except Exception as e:
                print e
                raise
        elif 'LON1 / LON2 / DLON' in label:
            try:
                lon1 = float(content[2:8])
                lon2 = float(content[8:14])
                dlon = float(content[14:20])

                header_info['lon1'] = lon1
                header_info['lon2'] = lon2
                header_info['dlon'] = dlon
            except Exception as e:
                print e
                raise
        elif 'EXPONENT' in label:
            try:
                exponent = int(content[0:6])
                header_info['exponent'] = exponent
            except Exception as e:
                print e
                raise
        elif 'PRN / BIAS / RMS' in label:
            # DETERMINE IF THE BIAS IS FOR GPS OR GLONASS
            gnss_character = content[3]
            if gnss_character == ' ' or gnss_character == 'G':
                gnss_system = 'GPS'
            elif gnss_character == 'R':
                gnss_system = 'GLONASS'

            try:
                PRN  = int(content[4:6])
            except Exception as e:
                print e
                raise

            try:
                bias = float(content[6:16])
            except Exception as e:
                print e
                raise

            try:
                rms = float(content[16:26])
            except Exception as e:
                print e
                raise

            satellite_biases[gnss_system][PRN] = (bias, rms)
        elif 'STATION / BIAS / RMS' in label:

            gnss_character = content[3]
            if gnss_character == ' ' or gnss_character == 'G':
                gnss_system = 'GPS'
            elif gnss_character == 'R':
                gnss_system = 'GLONASS'

            station = content[6:10]
            domes   = content[11:20]

            if domes.strip() == '':
                domes = None

            try:
                bias = float(content[26:36])
            except Exception as e:
                print e
                raise

            try:
                rms  = float(content[36:46])
            except Exception as e:
                print e
                raise

            station_biases[gnss_system][station] = (bias, rms, domes)

    return header_info, satellite_biases, station_biases

def parser(path_to_file):
    """
    Returns the raw data from the specified IONEX file path_to_file.
    (lons, lats, heights, tec_maps, rms_maps, height_maps, satellite_biases, station_biases)

        lons: 1D array containing the longitude grid used in the IONEX file
        lats: 1D array containing the latitude grid used in the IONEX file
        heights: 1D array containiner the height grid used in the IONEX file
        tec_maps: dictionary, keyed by the map_ids specified in the IONEX file
                  (typically, this starts at "1" and monotonically increases)

                  - each dictionary entry is a 2-tuple where
                        * the first entry is the epoch of the map as a datetime object
                        * the second entry is a 3D array of shape (n_lons, n_lats, n_heights)
                          containing the TEC data.
                  - If an epoch is "None", no epoch was provided for that map.

        rms_maps: dictionary, with the same specification as tec_maps except with RMS data.
        height_maps: dictionary, with the same specification as tec_maps except with height map data.
        satellite_biases: dictionary with keys "GPS" and "GLONASS"
                  - satellite_biases['GPS'] contains another dictionary (keyed by PRN). Each
                  dictionary value is a 2-tuple, where the first element is the satellite bias,
                  and the second element is the rms.
                  - satellite_biases['GLONASS'] is the same as the GPS case, but for GLONASS.
        station_biases: dictionary with keys "GPS" and "GLONASS"
                  - station_biases['GPS'] contains another dictionary (keyed by site ID).
                  Each dictionary value is a 3-tuple, where the first element is the bias,
                  the second element is the rms, and the third element is the DOMES number.
                  - station_biases['GLONASS'] is the same as the GPS case, but for GLONASS.
                  The description of the DOMES numbering can be found
                  in MERIT/COTES JOINT WORKING GROUPS, MERIT CAMPAIGN: CONNECTION
                  OF REFERENCE FRAMES, IMPLEMENTATION PLAN, 1983
    """

    with open(path_to_file, 'rb') as ionex_file:

        header_info, satellite_biases, station_biases = read_header(ionex_file)

        starting_epoch = header_info['start_epoch']
        stopping_epoch = header_info['stop_epoch']

        # import time
        # starting_sec   = time.mktime(starting_epoch.timetuple())
        # stopping_sec   = time.mktime(stopping_epoch.timetuple())

        starting_ht  = header_info['height1']
        stopping_ht  = header_info['height2']
        dh           = header_info['dh']
        starting_lat = header_info['lat1']
        stopping_lat = header_info['lat2']
        dlat         = header_info['dlat']
        starting_lon = header_info['lon1']
        stopping_lon = header_info['lon2']
        dlon         = header_info['dlon']

        exponent     = header_info['exponent']

        if dh == 0:
            n_heights = 1
        else:
            n_heights = int((stopping_ht - starting_ht)/dh) + 1

        if dlat == 0:
            n_lats = 1
        else:
            n_lats = int((stopping_lat - starting_lat)/dlat) + 1

        if dlon == 0:
            n_lons = 1
        else:
            n_lons = int((stopping_lon - starting_lon)/dlon) + 1

        if n_lons > 1:
            lons = np.arange(starting_lon, stopping_lon + dlon, dlon)
        else:
            lons = np.array([starting_lon])

        if n_lats > 1:
            lats = np.arange(starting_lat, stopping_lat + dlat, dlat)
        else:
            lats = [starting_lat]

        if n_heights > 1:
            heights = np.arange(starting_ht, stopping_ht + dh, dh)
        else:
            heights = [starting_ht]

        # maps are keyed by their map_id #
        tec_maps    = {}
        rms_maps    = {}
        height_maps = {}

        done = False
        while not done:

            # CHECK WHICH TYPE OF MAP THE NEXT DATA BLOCK IS FOR #
            # MAPS ARE STORED IN A DICTIONARY THAT IS KEYED BY THE ID OF THE MAP #
            label, content = read_misc_line(ionex_file)

            if label is None:
                done = True
            elif 'START OF TEC MAP' in label:
                split = content.split()
                try:
                    map_id = int(split[0])
                except Exception as e:
                    print e
                    raise
                epoch, data_block = read_map(ionex_file, n_lons, n_lats, n_heights, exponent)
                tec_maps[map_id] = (epoch, data_block) 
            elif 'START OF RMS MAP' in label:
                split = content.split()
                try:
                    map_id = int(split[0])
                except Exception as e:
                    print e
                    raise
                epoch, data_block = read_map(ionex_file, n_lons, n_lats, n_heights, exponent)
                rms_maps[map_id] = (epoch, data_block)
            elif 'START OF HEIGHT MAP' in label:
                split = content.split()
                try:
                    map_id = int(split[0])
                except Exception as e:
                    print e
                    raise
                epoch, data_block = read_map(ionex_file, n_lons, n_lats, n_heights, exponent)
                height_maps[map_id] = (epoch, data_block)

    return lons, lats, heights, tec_maps, rms_maps, height_maps, satellite_biases, station_biases

def raw(path_to_file):
    """
    raw:
    Reads in the IONEX file path_to_file and returns the data.

    The return value will be a 7-tuple. The first two entires are 1D numpy arrays of 
    the longitude and latitude grid, respectively. The third entry in the tuple is
    a list of datetimes that define the time grid. The fourth entry in the tuple is
    a 3D numpy array with shape (n_lons, n_lats, n_times) containing TEC data. The fifth
    entry in the tuple is a 3D numpy array with shape (n_lons, n_lats, n_times)
    containing RMS data (if it was present in the file, otherwise it will be None).
    n_times is the size of the time grid in the IONEX file, and n_lons, n_lats are the
    size of the longitude and latitude grids in the IONEX file. The sixth and seventh
    entries are the satellite biases and station biases, respectively. 
    satellite_biases: dictionary with keys "GPS" and "GLONASS"
              - satellite_biases['GPS'] contains another dictionary (keyed by PRN). Each
              dictionary value is a 2-tuple, where the first element is the satellite bias,
              and the second element is the rms.
              - satellite_biases['GLONASS'] is the same as the GPS case, but for GLONASS.
    station_biases: dictionary with keys "GPS" and "GLONASS"
              - station_biases['GPS'] contains another dictionary (keyed by site ID).
              Each dictionary value is a 3-tuple, where the first element is the bias,
              the second element is the rms, and the third element is the DOMES number.
              - station_biases['GLONASS'] is the same as the GPS case, but for GLONASS.
              The description of the DOMES numbering can be found
              in MERIT/COTES JOINT WORKING GROUPS, MERIT CAMPAIGN: CONNECTION
              OF REFERENCE FRAMES, IMPLEMENTATION PLAN, 1983
    """

    lons, lats, heights, tec_maps, rms_maps, height_maps, satellite_biases, station_biases = parser(path_to_file)

    time_grid = [tec_maps[key][0] for key in tec_maps]

    tec_grids = np.zeros([len(lons), len(lats), len(time_grid)])

    for i, key in enumerate(tec_maps):
        tec_grids[:,:,i] = tec_maps[key][1][:,:,0]

    # rms grids are not required to be specified in the IONEX standard
    if len(rms_maps) > 0:
        rms_grids = np.zeros([len(lons), len(lats), len(time_grid)])

        for i, key in enumerate(tec_maps):
            rms_grids[:,:,i] = rms_maps[key][1][:,:,0]
    else:
        rms_grids = None

    return (lons, lats, time_grid, tec_grids, rms_grids, satellite_biases, station_biases)

def decreasing(array):
    if len(array) > 1 and array[1] < array[0]:
        return True
    else:
        return False

def interpolate2D_spatial(path_to_file, spatial_grid, method = 'linear'):
    """
    interpolate2D_spatial:
    Reads in the IONEX file path_to_file and performs spatial interpolation on the data in 2D.
    spatial_grid is interpreted as a 2-tuple. The first entry in the tuple should be a 1D list
    of longitudes. The second entry in the tuple should be a 1D list of latitudes.

    The return value will be a 5-tuple. The first entry in the tuple is a 1D list of datetimes that
    define the time grid. The last two entries are 3D numpy arrays with shape
    (n_lons, n_lats, n_times) where n_times is the length of the time grid in the IONEX file and
    n_lons, n_lats are the length of the longitude and latitude arrays specified in spatial_grid.
    If the RMS data is not available in the IONEX file, None is returned in place of the 3D numpy
    array. The fourth and fifth entries are the satellite biases and station biases, respectively.

    satellite_biases: dictionary with keys "GPS" and "GLONASS"
              - satellite_biases['GPS'] contains another dictionary (keyed by PRN). Each
              dictionary value is a 2-tuple, where the first element is the satellite bias,
              and the second element is the rms.
              - satellite_biases['GLONASS'] is the same as the GPS case, but for GLONASS.
    station_biases: dictionary with keys "GPS" and "GLONASS"
              - station_biases['GPS'] contains another dictionary (keyed by site ID).
              Each dictionary value is a 3-tuple, where the first element is the bias,
              the second element is the rms, and the third element is the DOMES number.
              - station_biases['GLONASS'] is the same as the GPS case, but for GLONASS.
              The description of the DOMES numbering can be found
              in MERIT/COTES JOINT WORKING GROUPS, MERIT CAMPAIGN: CONNECTION
              OF REFERENCE FRAMES, IMPLEMENTATION PLAN, 1983 

    spatial_method specifies the spatial interpolation type:
    'nearest' or 'linear'
    default: 'linear'
    """

    # TODO: unify outputs, write documentation, generate movie
    lons, lats, heights, tec_maps, rms_maps, height_maps, satellite_biases, station_biases = parser(path_to_file)

    time_grid = [tec_maps[key][0] for key in tec_maps]

    si_grid_tec   = np.zeros([len(spatial_grid[0]), len(spatial_grid[1]), len(time_grid)])

    if len(rms_maps) > 0:
        si_grid_rms   = np.zeros([len(spatial_grid[0]), len(spatial_grid[1]), len(time_grid)])
    else:
        si_grid_rms = None

    if spatial_grid is not None:

        # SPATIAL INTERPOLATION #
        from scipy import interpolate

        # FORCE THE DATA TO BE STRICTLY INCREASING #
        if decreasing(lons):
            lons = lons[::-1]
            for key in tec_maps:
                tec_maps[key][1][:,:,0] = tec_maps[key][1][::-1,:,0]
                if len(rms_maps) > 0:
                    rms_maps[key][1][:,:,0] = rms_maps[key][1][::-1,:,0]

        if decreasing(lats):
            lats = lats[::-1]
            for key in tec_maps:
                tec_maps[key][1][:,:,0] = tec_maps[key][1][:,::-1,0]
                if len(rms_maps) > 0:
                    rms_maps[key][1][:,:,0] = rms_maps[key][1][:,::-1,0]

        for i, key in enumerate(tec_maps):
            tec_map            = tec_maps[key][1][:,:,0]
            s_interpolator_tec = interpolate.RegularGridInterpolator((lons, lats), tec_map, method = method, bounds_error = False)

            if len(rms_maps) > 0:
                rms_map            = rms_maps[key][1][:,:,0]
                s_interpolator_rms = interpolate.RegularGridInterpolator((lons, lats), rms_map, method = method, bounds_error = False)

            DENSE_LATS, DENSE_LONS = np.meshgrid(spatial_grid[1], spatial_grid[0])
            points                 = zip(DENSE_LONS.flatten(), DENSE_LATS.flatten())

            s_interpolated_tec = s_interpolator_tec(points).reshape(DENSE_LATS.shape)
            si_grid_tec[:,:,i] = s_interpolated_tec

            if len(rms_maps) > 0:
                s_interpolated_rms = s_interpolator_rms(points).reshape(DENSE_LATS.shape)
                si_grid_rms[:,:,i] = s_interpolated_rms


    return time_grid, si_grid_tec, si_grid_rms, satellite_biases, station_biases


def interpolate2D_temporal(path_to_file, temporal_grid, method = "linear", data = 'tec'):
    """
    interpolate2D_temporal:
    Reads in the IONEX file path_to_file and performs temporal interpolation on the data.
    temporal_grid is interpreted as a list of datetime objects and defines the new grid on which
    to interpolate the data.

    The return value will be a 6-tuple, where the first two entries are 1D numpy arrays containing
    the longitude and latitude grid points (respectively). The next two entries are 3D numpy arrays
    with shape (n_lons, n_lats, n_times) of the TEC and RMS data, respectively. n_times is the length
    of temporal_grid. n_lons and n_lats are the length of the longitude and latitude grids in the 
    IONEX file. If the RMS data is not available in the IONEX file, None is returned in place of the
    3D numpy array. The fifth and sixth entries are the satellite biases and station biases, respectively.

    satellite_biases: dictionary with keys "GPS" and "GLONASS"
              - satellite_biases['GPS'] contains another dictionary (keyed by PRN). Each
              dictionary value is a 2-tuple, where the first element is the satellite bias,
              and the second element is the rms.
              - satellite_biases['GLONASS'] is the same as the GPS case, but for GLONASS.
    station_biases: dictionary with keys "GPS" and "GLONASS"
              - station_biases['GPS'] contains another dictionary (keyed by site ID).
              Each dictionary value is a 3-tuple, where the first element is the bias,
              the second element is the rms, and the third element is the DOMES number.
              - station_biases['GLONASS'] is the same as the GPS case, but for GLONASS.
              The description of the DOMES numbering can be found
              in MERIT/COTES JOINT WORKING GROUPS, MERIT CAMPAIGN: CONNECTION
              OF REFERENCE FRAMES, IMPLEMENTATION PLAN, 1983

    temporal_method specifies the temporal interpolation type:
    'linear', 'nearest', 'zero', 'slinear', 'quadratic', or 'cubic'
    default: 'linear'
    """

    lons, lats, heights, tec_maps, rms_maps, height_maps, satellite_biases, station_biases = parser(path_to_file)

    import calendar
    from scipy import interpolate

    k      = tec_maps.keys()[0]
    n_lons = tec_maps[k][1].shape[0]
    n_lats = tec_maps[k][1].shape[1]

    timestamps = [calendar.timegm(dt.timetuple()) for dt in temporal_grid]

    ti_grid_tec    = np.zeros([n_lons, n_lats, len(timestamps)])

    if len(rms_maps) > 0:
        ti_grid_rms    = np.zeros([n_lons, n_lats, len(timestamps)])

    for i in range(n_lons):
        for j in range(n_lats):
            pixel_times    = [calendar.timegm(tec_maps[key][0].timetuple()) for key in tec_maps]

            pixel_list_tec     = [tec_maps[key][1][i,j,0] for key in tec_maps]
            t_interpolator_tec = interpolate.interp1d(pixel_times, pixel_list_tec, kind = method, bounds_error = False)
            ti_grid_tec[i,j,:] = t_interpolator_tec(timestamps)

            if len(rms_maps) > 0:
                pixel_list_rms     = [rms_maps[key][1][i,j,0] for key in rms_maps]
                t_interpolator_rms = interpolate.interp1d(pixel_times, pixel_list_rms, kind = method, bounds_error = False)
                ti_grid_rms[i,j,:] = t_interpolator_rms(timestamps)      

    return lons, lats, ti_grid_tec, ti_grid_rms, satellite_biases, station_biases

def interpolate2D_spatiotemporal(path_to_file, temporal_grid, spatial_grid, temporal_method = "linear", spatial_method = "linear", data = 'tec'):
    """
    interpolate2D_spatiotemporal:
    Reads in the IONEX file path_to_file and performs spatiotemporal interpolation on the data in 2D.
    temporal_grid is interpreted as a list of datetime objects and defines the new grid on which
    to interpolate the data. spatial_grid is interpreted as a 2-tuple. The first entry in the tuple
    should be a 1D list of longitudes. The second entry in the tuple should be a 1D list of latitudes.

    The return value will be a 4-tuple. The first entry is a 3D numpy array with shape
    (n_lons, n_lats, n_times) containing the TEC data. The second entry is a 3D numpy array with shape
    (n_lons, n_lats, n_times) containing the RMS data. n_times is length of the temporal grid, and 
    n_lons, n_lats are the size of the longitude and latitude arrays specified in spatial_grid. If the
    RMS data is not available in the IONEX file, None is returned in place of the 3D numpy array.
    The third and fourth entries are the satellite biases and station biases, respectively.

    satellite_biases: dictionary with keys "GPS" and "GLONASS"
              - satellite_biases['GPS'] contains another dictionary (keyed by PRN). Each
              dictionary value is a 2-tuple, where the first element is the satellite bias,
              and the second element is the rms.
              - satellite_biases['GLONASS'] is the same as the GPS case, but for GLONASS.
    station_biases: dictionary with keys "GPS" and "GLONASS"
              - station_biases['GPS'] contains another dictionary (keyed by site ID).
              Each dictionary value is a 3-tuple, where the first element is the bias,
              the second element is the rms, and the third element is the DOMES number.
              - station_biases['GLONASS'] is the same as the GPS case, but for GLONASS.
              The description of the DOMES numbering can be found
              in MERIT/COTES JOINT WORKING GROUPS, MERIT CAMPAIGN: CONNECTION
              OF REFERENCE FRAMES, IMPLEMENTATION PLAN, 1983

    temporal_method specifies the temporal interpolation type:
    'linear', 'nearest', 'zero', 'slinear', 'quadratic', or 'cubic'
    default: 'linear'

    spatial_method specifies the spatial interpolation type:
    'nearest' or 'linear'
    default: 'linear'
    """

    lons, lats, heights, tec_maps, rms_maps, height_maps, satellite_biases, station_biases = parser(path_to_file)

    si_grid_tec  = np.zeros([len(spatial_grid[0]), len(spatial_grid[1]), len(tec_maps.keys())])

    if len(rms_maps) > 0:
        si_grid_rms   = np.zeros([len(spatial_grid[0]), len(spatial_grid[1]), len(tec_maps.keys())])
    else:
        si_grid_rms = None

   # SPATIAL INTERPOLATION #
    from scipy import interpolate

    # FORCE THE DATA TO BE STRICTLY INCREASING #
    if decreasing(lons):
        lons = lons[::-1]
        for key in tec_maps:
            tec_maps[key][1][:,:,0] = tec_maps[key][1][::-1,:,0]
            if len(rms_maps) > 0:
                rms_maps[key][1][:,:,0] = rms_maps[key][1][::-1,:,0]

    if decreasing(lats):
        lats = lats[::-1]
        for key in tec_maps:
            tec_maps[key][1][:,:,0] = tec_maps[key][1][:,::-1,0]
            if len(rms_maps) > 0:
                rms_maps[key][1][:,:,0] = rms_maps[key][1][:,::-1,0]

    for i, key in enumerate(tec_maps):
        tec_map            = tec_maps[key][1][:,:,0]
        s_interpolator_tec = interpolate.RegularGridInterpolator((lons, lats), tec_map, method = spatial_method, bounds_error = False)

        if len(rms_maps) > 0:
            rms_map            = rms_maps[key][1][:,:,0]
            s_interpolator_rms = interpolate.RegularGridInterpolator((lons, lats), rms_map, method = spatial_method, bounds_error = False)

        DENSE_LATS, DENSE_LONS = np.meshgrid(spatial_grid[1], spatial_grid[0])
        points                 = zip(DENSE_LONS.flatten(), DENSE_LATS.flatten())

        s_interpolated_tec = s_interpolator_tec(points).reshape(DENSE_LATS.shape)
        si_grid_tec[:,:,i] = s_interpolated_tec

        if len(rms_maps) > 0:
            s_interpolated_rms = s_interpolator_rms(points).reshape(DENSE_LATS.shape)
            si_grid_rms[:,:,i] = s_interpolated_rms

    # TEMPORAL INTERPOLATION #
    import calendar

    sti_grid_tec   = np.zeros([len(spatial_grid[0]), len(spatial_grid[1]), len(temporal_grid)])
    timestamps = [calendar.timegm(dt.timetuple()) for dt in temporal_grid]

    if len(rms_maps) > 0:
        sti_grid_rms   = np.zeros([len(spatial_grid[0]), len(spatial_grid[1]), len(temporal_grid)])
    else:
        sti_grid_rms   = None

    for i in range(0, len(spatial_grid[0])):
        for j in range(0, len(spatial_grid[1])):
            pixel_times         = [calendar.timegm(tec_maps[key][0].timetuple()) for key in tec_maps]

            pixel_list_tec      = si_grid_tec[i,j,:]
            t_interpolator_tec  = interpolate.interp1d(pixel_times, pixel_list_tec, kind = temporal_method, bounds_error = False)
            sti_grid_tec[i,j,:] = t_interpolator_tec(timestamps)

            if len(rms_maps) > 0:
                pixel_list_rms      = si_grid_rms[i,j,:]
                t_interpolator_rms  = interpolate.interp1d(pixel_times, pixel_list_rms, kind = temporal_method, bounds_error = False)
                sti_grid_rms[i,j,:] = t_interpolator_rms(timestamps)               

    return sti_grid_tec, sti_grid_rms, satellite_biases, station_biases