from datetime import datetime

import pandas as PD


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
        df = PD.DataFrame(data={'x': x, 'y': y, 'z': z, 'flag': flag},
                          index=dt)
        df.siteid = siteid
        df.lat = float(lat)
        df.lon = float(lon)
        df.date = datetime.strptime(date, '%Y%m%d')
        df.pos_format = pos_format
        df.units = units
        df.sample_rate = sample_rate
    return df
