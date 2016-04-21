import logging
from ftplib import FTP
from contextlib import closing
from tempfile import gettempdir
from gzip import GzipFile

from ..util.path import replace_path

logger = logging.getLogger('pyrsss.gps.sideshow')


SIDESHOW_SERVER = 'sideshow.jpl.nasa.gov'
"""
Name of JPL sideshow FTP server.
"""


def update_sideshow_file(fname,
                         server_fname,
                         server=SIDESHOW_SERVER,
                         temp_path=gettempdir()):
    """
    Update the JPL side show file stored locally at *fname*. The
    remote file is accessed via FTP on *server* at *server_fname*. The
    path *temp_path* is used to store intermediate files. Return
    *fname*.
    """
    dest_fname = replace_path(temp_path, server_fname)
    logger.info('opening connection to {}'.format(server))
    with closing(FTP(server)) as ftp, open(dest_fname, 'w') as fid:
        logger.info('logging in')
        ftp.login()
        logger.info('writing to {}'.format(dest_fname))
        ftp.retrbinary('RETR ' + server_fname, fid.write)
    logger.info('uncompressing file to {}'.format(fname))
    with GzipFile(dest_fname) as gzip_fid, open(fname, 'w') as fid:
        fid.write(gzip_fid.read())
    return fname
