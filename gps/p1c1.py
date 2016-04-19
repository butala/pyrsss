import logging
import os
from datetime import datetime, timedelta
from ftplib import FTP
from contextlib import closing
from tempfile import gettempdir
from gzip import GzipFile
from collections import OrderedDict

import numpy as NP

from constants import M_TO_TECU
from ..util.path import replace_path

logger = logging.getLogger('pyrsss.gps.p1c1')


P1C1_FNAME = os.path.join(os.path.dirname(__file__),
                          'CA-P')
"""
Full path to the receiver type file.
"""


P1C1_SERVER = 'sideshow.jpl.nasa.gov'
"""
FTP server for the P1-C1 DCB file.
"""


P1C1_SERVER_FNAME = '/pub/gipsy_files/gipsy_params/CA-P.gz'
"""
Path to the remote P1-C1 DCB file.
"""


def update_p1c1(p1c1_fname=P1C1_FNAME,
                p1c1_server=P1C1_SERVER,
                p1c1_server_fname=P1C1_SERVER_FNAME,
                temp_path=gettempdir()):
    """
    Update the receiver type table stored at *p1c1_fname*. The remote
    file is accessed via FTP on *p1c1_server* at
    *p1c1_server_fname*. The path *temp_path* is used to store
    intermediate files.
    """
    dest_fname = replace_path(temp_path, p1c1_server_fname)
    logger.info('opening connection to {}'.format(p1c1_server))
    with closing(FTP(p1c1_server)) as ftp, open(dest_fname, 'w') as fid:
        logger.info('logging in')
        ftp.login()
        logger.info('writing to {}'.format(dest_fname))
        ftp.retrbinary('RETR ' + p1c1_server_fname, fid.write)
    logger.info('uncompressing file to {}'.format(p1c1_fname))
    with GzipFile(dest_fname) as gzip_fid, open(p1c1_fname, 'w') as fid:
        fid.write(gzip_fid.read())
    return p1c1_fname


class P1C1Table(OrderedDict):
    def __init__(self, p1c1_fname=P1C1_FNAME):
        """
        Parse *p1c1_fname* and store DCB (in [TECU]) in the mapping
        :class:`datetime` -> ['svn', 'prn'] -> integer ID.
        """
        super(P1C1Table, self).__init__()
        with open(p1c1_fname) as fid:
            for line in fid:
                if line.startswith('#'):
                    continue
                cols = line.split()
                date = datetime.strptime(cols[0], '%Y-%m-%d')
                prn = int(cols[1])
                svn = int(cols[2])
                CA_P_m = float(cols[3])
                CA_P_tecu = CA_P_m * M_TO_TECU
                self.setdefault(date, {}).setdefault('prn', {})[prn] = CA_P_tecu
                self[date].setdefault('svn', {})[svn] = CA_P_tecu

    def __call__(self, date, delta=timedelta(days=32)):
        """
        Return the table entry closest to *date*. Check that the closest
        table entry is no greater than *delta* away.
        """
        diff = NP.array([abs((x - date).total_seconds()) for x in self.iterkeys()])
        I = NP.argmin(diff)
        closest_date = self.keys()[I]
        # make sure date argument is no further than 1 month away from
        # a table entry
        assert abs(date - closest_date) < delta
        return self[closest_date]


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    update_p1c1()

    p1c1_table = P1C1Table()

    date = datetime(2016, 1, 25)
    print(p1c1_table(date))
