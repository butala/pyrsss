import logging
import os
from collections import namedtuple
from ftplib import FTP
from contextlib import closing
from tempfile import gettempdir
from gzip import GzipFile

from ..util.path import replace_path

logger = logging.getLogger('pyrsss.gps.receiver_type')


RECEIVER_TYPE_FNAME = os.path.join(os.path.dirname(__file__),
                                   'GPS_Receiver_Types')
"""
Full path to the receiver type file.
"""


RECEIVER_TYPE_SERVER = 'sideshow.jpl.nasa.gov'
"""
FTP server for the receiver type file.
"""


RECEIVER_TYPE_SERVER_FNAME = '/pub/gipsy_files/gipsy_params/GPS_Receiver_Types.gz'
"""
Path to the remote receiver type file.
"""


def update_receiver_type(receiver_type_fname=RECEIVER_TYPE_FNAME,
                         receiver_type_server=RECEIVER_TYPE_SERVER,
                         receiver_type_server_fname=RECEIVER_TYPE_SERVER_FNAME,
                         temp_path=gettempdir()):
    """
    Update the receiver type table stored at
    *receiver_type_fname*. The remote file is accessed via FTP on
    *receiver_type_server* at *receiver_type_server_fname*. The path
    *temp_path* is used to store intermediate files.
    """
    dest_fname = replace_path(temp_path, receiver_type_server_fname)
    logger.info('opening connection to {}'.format(receiver_type_server))
    with closing(FTP(receiver_type_server)) as ftp, open(dest_fname, 'w') as fid:
        logger.info('logging in')
        ftp.login()
        logger.info('writing to {}'.format(dest_fname))
        ftp.retrbinary('RETR ' + receiver_type_server_fname, fid.write)
    logger.info('uncompressing file to {}'.format(receiver_type_fname))
    with GzipFile(dest_fname) as gzip_fid, open(receiver_type_fname, 'w') as fid:
        fid.write(gzip_fid.read())
    return receiver_type_fname


class ReceiverTypeInfo(namedtuple('ReceiverTypeInfo', 'c1p1 fixtags igs')):
    pass


def parse_receiver_types(fname=RECEIVER_TYPE_FNAME):
    """
    Parse the receiver type file *fname*. Return a mapping between
    receiver type and :class:`ReceiverTypeInfo` classifying the
    receiver.
    """
    receiver_map = {}
    with open(fname) as fid:
        for line in fid:
            if line.startswith('#'):
                continue
            receiver_map[line[:20].rstrip()] = map(int, line[20:].split()[:3])
    return receiver_map


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    update_receiver_type()

    receiver_types = parse_receiver_types()
    print(len(receiver_types))
