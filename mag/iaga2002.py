import logging
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import datetime, timedelta
from collections import OrderedDict, namedtuple, defaultdict

logger = logging.getLogger('pyrsss.mag.iaga2002')


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


def parse(fname):
    """
    Parser the IAGA2002 format file *fname* and return a tuple with a
    :class:`Header` and mapping of date/times to measured values.
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
            # print(header_map.keys())
            header = Header(**header_map)
        except TypeError:
            logger.warning('unknown header record found in {} --- setting header to None'.format(fname))
            header = None
        # parse data header record
        fields = line[:69].split()
        if len(fields) != 7:
            raise RuntimeError('malformed data header record in {} ({})'.format(fname,
                                                                                line))
        DataRecord = namedtuple('DataRecord', ' '.join(fields[3:]))
        data_map = OrderedDict()
        # parse data records
        for line in fid:
            dt = datetime.strptime(line[:23], '%Y-%m-%d %H:%M:%S.%f')
            data1 = convert_float(line[31:40])
            data2 = convert_float(line[41:50])
            data3 = convert_float(line[51:60])
            data4 = convert_float(line[61:70])
            data_map[dt] = DataRecord(data1,
                                      data2,
                                      data3,
                                      data4)
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
