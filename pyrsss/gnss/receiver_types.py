import logging
import os
from collections import namedtuple

from ..util.path import replace_path
from sideshow import update_sideshow_file


RECEIVER_TYPES_FNAME = os.path.join(os.path.dirname(__file__),
                                   'GPS_Receiver_Types')
"""
Full path to the receiver types file.
"""


RECEIVER_TYPES_SERVER_FNAME = '/pub/gipsy_files/gipsy_params/GPS_Receiver_Types.gz'
"""
Path to the remote receiver types file.
"""


def update_receiver_types(receiver_types_fname=RECEIVER_TYPES_FNAME,
                          receiver_types_server_fname=RECEIVER_TYPES_SERVER_FNAME):
    """
    Update the receiver types table stored at
    *receiver_types_fname*. The remote file is accessed via FTP at
    *receiver_types_server_fname*.
    """
    return update_sideshow_file(receiver_types_fname,
                                receiver_types_server_fname)


class ReceiverTypeInfo(namedtuple('ReceiverTypeInfo', 'c1p1 fixtags igs')):
    pass


class ReceiverTypes(dict):
    def __init__(self, fname=RECEIVER_TYPES_FNAME):
        """
        Parse the receiver types file *fname and store the mapping between
        receiver type and :class:`ReceiverTypeInfo` classifying the
        receiver.
        """
        if fname == RECEIVER_TYPES_FNAME and not os.path.isfile(fname):
            update_receiver_types()
        with open(fname) as fid:
            for line in fid:
                if line.startswith('#') or len(line.strip()) == 0:
                    continue
                self[line[:20].rstrip()] = ReceiverTypeInfo(*map(int, line[20:].split()[:3]))


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    update_receiver_types()
    receiver_types = ReceiverTypes()

    print(len(receiver_types))
