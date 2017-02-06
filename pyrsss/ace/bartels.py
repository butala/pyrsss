import os
import logging
from urllib2 import urlopen
from contextlib import closing
from collections import OrderedDict, namedtuple
from datetime import datetime, timedelta

from intervals import DateTimeInterval

logger = logging.getLogger('pyrsss.ace.bartels')


BARTELS_URL = 'http://www.srl.caltech.edu/ACE/ASC/DATA/bartels/bartels.txt'
""" ??? """


BARTELS_FNAME = os.path.join(os.path.dirname(__file__),
                           'bartels.txt')
""" ??? """


def update(local_fname=BARTELS_FNAME,
           bartels_url=BARTELS_URL):
    """
    """
    with open(local_fname, 'w') as fid, closing(urlopen(bartels_url)) as fid_in:
        logger.info('fetching {}'.format(bartels_url))
        fid.write(fid_in.read())
    return local_fname


def initialize(local_fname=BARTELS_FNAME):
    """
    ???
    """
    if os.path.isfile(local_fname):
       logger.warning('{} exists --- skipping initialization'.format(local_fname))
       return local_fname
    return update(local_fname=local_fname)


def efloat(x):
    """
    """
    if x.startswith('e'):
        return float(x[1:]), True
    elif x.strip() == '-':
        return None, False
    else:
        return float(x), False


class BartelsField(namedtuple('BartelsField',
                              ' '.join(['dt',
                                        'start',
                                        'doy',
                                        'ut_sec',
                                        'ace_epoch',
                                        'spacecraft_clock',
                                        'spacecraft_clock_estimated']))):
                    pass


def parse(local_fname=BARTELS_FNAME,
          types=[int,
                 lambda x: datetime.strptime(x, '%m/%d/%Y'),
                 int,
                 float,
                 float,
                 efloat],
          data_map=None):
    """
    """
    if data_map is None:
        data_map = OrderedDict()
    with open(local_fname) as fid:
        for line in fid:
            if line.startswith('----'):
                break
        for line in fid:
            if line.startswith('Pre L1 Insertion:'):
                continue
            if line.startswith('Post L1 Insertion:'):
                continue
            if len(line.strip()) == 0:
                continue
            if line.startswith(' *-Seconds since 00:00:00 UT on 1/1/96'):
                continue
            if line.startswith(' e-Estimated Times'):
                continue
            data = [convert(x) for x, convert in zip(line.split(), types)]
            rot_n = data[0]
            data_out = [data[1] + timedelta(seconds=data[3])]
            data_out.extend(data[1:-1])
            data_out.extend(data[-1])
            data_map[rot_n] = BartelsField(*data_out)
    return data_map


class Bartels(dict):
    def __new__(cls, local_fname=BARTELS_FNAME):
        """
        """
        obj = super(Bartels, cls).__new__(cls)
        parse(local_fname=local_fname, data_map=obj)
        obj._interval_map = OrderedDict()
        times_1 = [x.dt for x in obj.values()[:-1]]
        times_2 = [x.dt for x in obj.values()[1:]]
        for b_id, (t1, t2) in zip(obj,
                                  zip(times_1,
                                      times_2)):
            interval = DateTimeInterval.closed_open(t1, t2)
            obj._interval_map[interval] = b_id
        return obj


    def __call__(self, dt):
        """
        """
        for interval, b_id in self._interval_map.iteritems():
            if dt in interval:
                return b_id
        return None
