import logging
import os
from datetime import datetime, date
from collections import namedtuple

import portion as P

from .sideshow import update_sideshow_file


PRN_GPS_FNAME = os.path.join(os.path.dirname(__file__),
                               'PRN_GPS')
"""
Full path to the SVN/PRN mapping file.
"""


PRN_GPS_SERVER_FNAME = '/pub/gipsy_files/gipsy_params/PRN_GPS.gz'
"""
Path to the remote SVN/PRN mapping file.
"""


def update_prn_gps(prn_gps_fname=PRN_GPS_FNAME,
                   prn_gps_server_fname=PRN_GPS_SERVER_FNAME):
    """
    Update the SVN/PRN mapping table stored at *prn_gps_fname*. The
    remote file is accessed via FTP at *prn_gps_server_fname*.
    """
    return update_sideshow_file(prn_gps_fname,
                                prn_gps_server_fname)


class Info(namedtuple('Info',
                      'launch '
                      'deactivation '
                      'svn '
                      'prn '
                      'block '
                      'orbit '
                      'clock')):
    pass


class Table(list):
    def __init__(self, prn_gps_fname=PRN_GPS_FNAME):
        """
        ???
        """
        super(Table, self).__init__()
        if prn_gps_fname == PRN_GPS_FNAME and not os.path.isfile(prn_gps_fname):
            update_prn_gps()
        with open(prn_gps_fname) as fid:
            for i, line in enumerate(fid):
                if i == 0:
                    # skip header line
                    continue
                cols = line.split()
                launch = datetime.strptime(cols[0], '%Y-%m-%d').date()
                if cols[1] == '0000':
                    deactivate = datetime.utcnow().date()
                else:
                    deactivate = datetime.strptime(cols[1], '%Y-%m-%d').date()
                svn = int(cols[2])
                prn = int(cols[3])
                block = cols[4]
                if len(cols) >= 6:
                    orbit = cols[5]
                else:
                    orbit = None
                if len(cols) >= 7:
                    clock = cols[6:]
                else:
                    clock = None
                self.append(Info(launch,
                                 deactivate,
                                 svn,
                                 prn,
                                 block,
                                 orbit,
                                 clock))

    def prn(self, svn, date):
        """
        ???
        """
        for candidate in filter(lambda x: x.svn == svn, self):
            if date in P.closed(candidate.launch, candidate.deactivation):
                return candidate.prn
        raise RuntimeError('could not find PRN associated with SVN={} on '
                           '{:%Y-%m-%d}'.format(svn,
                                                date))


    def svn(self, prn, date):
        """
        ???
        """
        for candidate in filter(lambda x: x.prn == prn, self):
            if date in P.closed(candidate.launch, candidate.deactivation]):
                return candidate.svn
        raise RuntimeError('could not find SVN associated with PRN={} on '
                           '{:%Y-%m-%d}'.format(prn,
                                                date))


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    update_prn_gps()
    prn_gps = Table()

    print(prn_gps.svn(7, date(2007, 1, 1)))
