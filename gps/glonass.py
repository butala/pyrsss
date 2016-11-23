import os
from collections import namedtuple, OrderedDict
from datetime import datetime

from intervals import DateTimeInterval

from sideshow import update_sideshow_file


GLO_STATUS_FNAME = os.path.join(os.path.dirname(__file__),
                                'GLO_STATUS')
"""
Full path to the GLONASS status file.
"""


GLO_STATUS_SERVER_FNAME = '/pub/gipsy_files/gipsy_params/GLO_STATUS.gz'
"""
Path to the remote GLONASS status file.
"""


def update_glo_status(glo_status_fname=GLO_STATUS_FNAME,
                      glo_status_server_fname=GLO_STATUS_SERVER_FNAME):
    """
    Update the GLONASS status table stored at *glo_status_fname*. The
    remote file is accessed via FTP at *glo_status_server_fname*.
    """
    return update_sideshow_file(glo_status_fname,
                                glo_status_server_fname)


class StatusInfo(namedtuple('StatusInfo', 'launch slot freq plane GLONASS cosmos')):
    pass


class GLONASS_Status(dict):
    def __init__(self, glo_status_fname=GLO_STATUS_FNAME):
        """
        Parse *glo_status_fname* and store GLONASS status information.
        """
        super(GLONASS_Status, self).__init__()
        if glo_status_fname == GLO_STATUS_FNAME and not os.path.isfile(glo_status_fname):
            update_glo_status()
        def parse_dt(date, time):
            if date == '0000-00-00' and time == '00:00':
                return None
            else:
                return datetime.strptime(date + ' ' + time,
                                         '%Y-%m-%d %H:%M')
        with open(glo_status_fname) as fid:
            for line in fid:
                if line.startswith('#'):
                    continue
                toks = line.split()
                launch_dt = parse_dt(toks[0], toks[1])
                start_dt = parse_dt(toks[2], toks[3])
                end_dt = parse_dt(toks[4], toks[5])
                slot, freq, plane, GLONASS, cosmos = map(int, toks[6:])
                interval = DateTimeInterval.closed_open(start_dt, end_dt)
                info = StatusInfo(launch_dt, slot, freq, plane, GLONASS, cosmos)
                self.setdefault(slot, OrderedDict())[interval] = info

    def __call__(self, slot, dt):
        """
        Return the :class:`StatusInfo` associated with GLONASS satellite
        with ID *slot* at :class:`datetime` *dt*.
        """
        for interval, info in self[slot].iteritems():
            if dt in interval:
                return info
        raise KeyError('no record for {} at {:%Y-%m-%d %H:%M} found'.format(slot, dt))


if __name__ == '__main__':
    glonass_status = GLONASS_Status()

    print(glonass_status(10, datetime(2014, 1, 1)))
