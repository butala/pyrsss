from datetime import datetime, date
from collections import namedtuple

from intervals import DateInterval


URL = 'ftp://sideshow.jpl.nasa.gov/pub/gipsy_files/gipsy_params/CA-P.gz'


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
    def __init__(self, prn_gps_fname):
        """
        ???
        """
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
            if date in DateInterval([candidate.launch,
                                     candidate.deactivation]):
                return candidate.prn
        raise RuntimeError('could not find PRN associated with SVN={} on '
                           '{:%Y-%m-%d}'.format(svn,
                                                date))


    def svn(self, prn, date):
        """
        ???
        """
        for candidate in filter(lambda x: x.prn == prn, self):
            if date in DateInterval([candidate.launch,
                                     candidate.deactivation]):
                return candidate.svn
        raise RuntimeError('could not find SVN associated with PRN={} on '
                           '{:%Y-%m-%d}'.format(prn,
                                                date))


if __name__ == '__main__':
    prn_gps = Table('/Users/butala/Desktop/PRN_GPS')

    print(prn_gps.svn(7, date(2007, 1, 1)))
