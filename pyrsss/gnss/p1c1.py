import logging
import os
from datetime import datetime, timedelta
from collections import OrderedDict

import numpy as np

from ..util.path import replace_path
from .sideshow import update_sideshow_file


P1C1_FNAME = os.path.join(os.path.dirname(__file__),
                          'CA-P')
"""
Full path to the receiver type file.
"""


P1C1_SERVER_FNAME = '/pub/gipsy_files/gipsy_params/CA-P.gz'
"""
Path to the remote P1-C1 DCB file.
"""


def update_p1c1(p1c1_fname=P1C1_FNAME,
                p1c1_server_fname=P1C1_SERVER_FNAME):
    """
    Update the receiver type table stored at *p1c1_fname*. The remote
    file is accessed via FTP at *p1c1_server_fname*.
    """
    return update_sideshow_file(p1c1_fname,
                                p1c1_server_fname)


class P1C1Table(OrderedDict):
    def __init__(self, p1c1_fname=P1C1_FNAME):
        """
        Parse *p1c1_fname* and store DCB (in [TECU]) in the mapping
        :class:`datetime` -> ['svn', 'prn'] -> integer ID.
        """
        super(P1C1Table, self).__init__()
        if p1c1_fname == P1C1_FNAME and not os.path.isfile(p1c1_fname):
            update_p1c1()
        with open(p1c1_fname) as fid:
            for line in fid:
                if line.startswith('#'):
                    continue
                cols = line.split()
                date = datetime.strptime(cols[0], '%Y-%m-%d')
                prn = int(cols[1])
                svn = int(cols[2])
                CA_P_m = float(cols[3])
                self.setdefault(date, {}).setdefault('prn', {})[prn] = CA_P_m
                self[date].setdefault('svn', {})[svn] = CA_P_m

    def __call__(self, date, delta=timedelta(days=32)):
        """
        Return the table entry closest to *date*. Check that the closest
        table entry is no greater than *delta* away.
        """
        diff = np.array([abs((x - date).total_seconds()) for x in self.iterkeys()])
        I = np.argmin(diff)
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
