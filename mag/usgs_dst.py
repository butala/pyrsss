import os
from itertools import izip
from datetime import datetime, timedelta

import pandas as PD

from repository import get_root


def date_parser(*cols):
    def get_int(x):
        return int(x[:-10])
    return [datetime(*map(get_int, [year, month, day])) + \
            timedelta(minutes=get_int(minute) - 1) for (month,
                                                        day,
                                                        year,
                                                        minute) in izip(*cols)]


def read_dst(fname=os.path.join(get_root(),
                                'USGS_dst',
                                'Dst_definitive_fix_minute.out')):
    df = PD.read_csv(fname,
                     header=None,
                     names=['month', 'day', 'year', 'minute', 'doy', 'dst'],
                     delim_whitespace=True,
                     engine='c',
                     parse_dates=[[0, 1, 2, 3]],
                     usecols=['month', 'day', 'year', 'minute', 'dst'],
                     date_parser=date_parser)
    df.columns = ['dt', 'dst']
    return PD.Series(data=usgs_dst['dst'].values,
                     index=usgs_dst['dt'].values)
