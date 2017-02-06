import sys
import logging

import pandas as PD

from pyhdf.HDF import HDF, HDF4Error
from pyhdf import VS


def parse_ace_data(hdf5_fname, hdf4_fname, N=1000):
    """
    """
    key = key_from_fname(hdf4_fname)
    with closing(HDF(hdf4_fname)) as hdf:
        vs = hdf.vstart()
        vdata = vs.attach(key)
        fieldinfo = vdata.fieldinfo()
        loop_divmod = divmod(vdata.inquire()[0], N)
        fields = [x[0] for x in fieldinfo]
        data_map = defaultdict(list)
        for i in range(loop_divmod[0] + 1):
            try:
                data = vdata.read(N if i < loop_divmod[0] else loop_divmod[1])
            except HDF4Error:
                break
            for data_i in data:
                for data_ii, field in zip(data_i, fields):
                    data_map[field].append(data_ii)
    # convert to DataFrame
    remove_set = set(['year',
                      'fp_year',
                      'day',
                      'fp_doy',
                      'hr',
                      'min',
                      'sec',
                      'ACEepoch'])
    dt = []
    for year, day, hr, minute, sec in zip(*[data_map[x] for x in ['year',
                                                                  'day',
                                                                  'hr',
                                                                  'min',
                                                                  'sec']]):
        dt.append(datetime(year, 1, 1) + timedelta(days=day - 1,
                                                   hours=hr,
                                                   minutes=minute,
                                                   seconds=sec))
    data = {k: v for k, v in data_map.iteritems() if k not in remove_set}
    df = PD.DataFrame(index=dt,
                      data=data)
    return df


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Convert ACE HDF4 file to pandas HDF5 record.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('hdf5_fname',
                        type=str,
                        help='output HDF5 file')
    parser.add_argument('hdf4_fname',
                        type=str,
                        help='input ACE HDF4 data record')
    args = parser.parse_args(argv[1:])




if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
