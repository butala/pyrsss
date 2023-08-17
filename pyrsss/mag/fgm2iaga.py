import sys
import logging
import os
from datetime import datetime
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import numpy as np
import pandas as pd

logger = logging.getLogger('pyrsss.mag.fm2iaga')


HEADER_TEMPLATE = """\
 Format                 IAGA-2002                                    |
 Source of Data         CARISMA                                      |
 IAGA CODE              {stn}                                          |
 Geodetic Latitude      {lat:<8.3f}                                     |
 Geodetic Longitude     {lon:<8.3f}                                     |
 Elevation              {el:<8.3f}                                     |
 Reported               XYZF                                         |
DATE       TIME         DOY     {stn}X       {stn}Y     {stn}Z      {stn}F   |
"""



def parse(fname):
    """
    Parse FGM format data *fname*. Return :class:`DataFrame`
    containing all information found in the file.

    The FGM file format is used by CARISMA to store data and is
    described here:
    http://www.carisma.ca/carisma-data/fgm-data-format.
    """
    with open(fname) as fid:
        siteid, lat, lon, date, pos_format, units, sample_rate = fid.next().split()
        dt = []
        x = []
        y = []
        z = []
        flag = []
        for line in fid:
            cols = line.split()
            dt.append(datetime.strptime(cols[0], '%Y%m%d%H%M%S'))
            x.append(float(cols[1]))
            y.append(float(cols[2]))
            z.append(float(cols[3]))
            if cols[4] == '.':
                flag.append(False)
            elif cols[4] == 'x':
                flag.append(True)
            else:
                raise ValueError('unknown flag value {} encountered in {}'.format(cols[4], fname))
        f = np.hypot(x, np.hypot(y, z))
        df = pd.DataFrame(data={'x': x, 'y': y, 'z': z, 'f': f, 'flag': flag},
                          index=dt)
        df.siteid = siteid
        df.lat = float(lat)
        df.lon = float(lon)
        df.date = datetime.strptime(date, '%Y%m%d')
        df.pos_format = pos_format
        df.units = units
        df.sample_rate = sample_rate
    return df




def fgm2iaga(path,
             fgm_fname,
             ftype='v',
             output_template='{stn}{date:%Y%m%d}{ftype}{interval}.{interval}'):
    """
    Parse FGM format file *fgm_fname* and reformat it to IAGA2002 and
    save at *path* (using *output_tempalte* to form the file
    name). Return the file name. The *ftype* denotes the file type: p
    - provisional, d - definitive, q - quasi-definitive, or v -
    variation.
    """
    df = parse(fgm_fname)
    delta = (df.index[1] - df.index[0]).total_seconds()
    if delta == 1.0:
        interval = 'sec'
    elif delta == 60.0:
        interval = 'min'
    else:
        raise ValueError('unknown data interval found in {}'.format(fgm_fname))
    stn = df.siteid[:3].upper()
    out_fname = os.path.join(path,
                             output_template.format(stn=stn.lower(),
                                                    date=df.date,
                                                    ftype=ftype,
                                                    interval=interval))
    with open(out_fname, 'w') as fid:
        fid.write(HEADER_TEMPLATE.format(stn=stn.upper(),
                                         lat=df.lat,
                                         lon=df.lon,
                                         el=0))
        for row in df.itertuples():
            dt = row.Index
            if row.flag:
                X = Y = Z = F = 99999
            else:
                X = row.x
                Y = row.y
                Z = row.z
                F = np.linalg.norm([X, Y, Z])
            fid.write('{date:%Y-%m-%d %H:%M:%S.000} {date:%j}'
                      '    {X:>9.2f} {Y:>9.2f} {Z:>9.2f} {F:>9.2f}\n'.format(date=dt,
                                                                             X=X,
                                                                             Y=Y,
                                                                             Z=Z,
                                                                             F=F))
    return out_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Convert FGM format data (CARISMA) to IAGA2002 format.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_path',
                        type=str,
                        help='path to store daily IAGA2002 format files')
    parser.add_argument('fgm_fnames',
                        type=str,
                        nargs='+',
                        metavar='fgm_fname',
                        help='FGM format file')
    args = parser.parse_args(argv[1:])


    for fgm_fname in args.fgm_fnames:
        iaga_fname = fgm2iaga(args.output_path, fgm_fname)
        logger.info('{} -> {}'.format(fgm_fname, iaga_fname))


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
