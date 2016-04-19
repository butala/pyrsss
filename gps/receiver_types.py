import logging
import os
from collections import namedtuple
from ftplib import FTP
from contextlib import closing
from tempfile import gettempdir
from gzip import GzipFile

from ..util.path import replace_path

logger = logging.getLogger('pyrsss.gps.receiver_types')


RECEIVER_TYPES_FNAME = os.path.join(os.path.dirname(__file__),
                                   'GPS_Receiver_Types')
"""
Full path to the receiver types file.
"""


RECEIVER_TYPES_SERVER = 'sideshow.jpl.nasa.gov'
"""
FTP server for the receiver types file.
"""


RECEIVER_TYPES_SERVER_FNAME = '/pub/gipsy_files/gipsy_params/GPS_Receiver_Types.gz'
"""
Path to the remote receiver types file.
"""


def update_receiver_types(receiver_types_fname=RECEIVER_TYPES_FNAME,
                          receiver_types_server=RECEIVER_TYPES_SERVER,
                          receiver_types_server_fname=RECEIVER_TYPES_SERVER_FNAME,
                          temp_path=gettempdir()):
    """
    Update the receiver types table stored at
    *receiver_types_fname*. The remote file is accessed via FTP on
    *receiver_types_server* at *receiver_types_server_fname*. The path
    *temp_path* is used to store intermediate files.
    """
    dest_fname = replace_path(temp_path, receiver_types_server_fname)
    logger.info('opening connection to {}'.format(receiver_types_server))
    with closing(FTP(receiver_types_server)) as ftp, open(dest_fname, 'w') as fid:
        logger.info('logging in')
        ftp.login()
        logger.info('writing to {}'.format(dest_fname))
        ftp.retrbinary('RETR ' + receiver_types_server_fname, fid.write)
    logger.info('uncompressing file to {}'.format(receiver_types_fname))
    with GzipFile(dest_fname) as gzip_fid, open(receiver_types_fname, 'w') as fid:
        fid.write(gzip_fid.read())
    return receiver_types_fname


class ReceiverTypeInfo(namedtuple('ReceiverTypesInfo', 'c1p1 fixtags igs')):
    pass


class ReceiverTypes(dict):
    def __init__(self, fname=RECEIVER_TYPES_FNAME):
        """
        Parse the receiver types file *fname and store the mapping between
        receiver type and :class:`ReceiverTypeInfo` classifying the
        receiver.
        """
        with open(fname) as fid:
            for line in fid:
                if line.startswith('#'):
                    continue
                self[line[:20].rstrip()] = ReceiverTypeInfo(map(int, line[20:].split()[:3]))


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    update_receiver_types()

    receiver_types = ReceiverTypes()
    print(len(receiver_types))
