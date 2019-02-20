import sys
import os
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from zipfile import ZipFile

from .index import initialize
from ..util.path import touch_path


EPILOG = """Use the following procedure to mirror the 3-D EMTFs stored at the
SPUDS repository (http://ds.iris.edu/spud/emtf):

1. Download sets (via Download SPUDS bundle) of EMTF records, 200 each
maximum and store at the same location.\n

2. Run this script to extract the pertinent information from the .zip
records and validate the expected number of EMTF records were found.
"""

logger = logging.getLogger('pyrsss.usarray_emtf.mirror')


def mirror(destination_path,
           zip_path,
           expected):
    """
    Process SPUDS bundles found at *zip_path* and store the XML and
    PNG contents to *destination_path*. Verify that the number of EMTF
    records processed is equal to *expected*. Return the full path to
    the pickle file containing the index of EMTF records.
    """
    zip_fnames = [x for x in os.listdir(zip_path) if x.startswith('SPUD') and x.endswith('.zip')]
    logger.info('found {} SPUDS .zip files at {}:'.format(len(zip_fnames),
                                                          zip_path))
    for zip_fname in zip_fnames:
        logger.info(zip_fname)
    xml_target_path = os.path.join(destination_path, 'xml')
    png_target_path = os.path.join(destination_path, 'png')
    touch_path(xml_target_path)
    touch_path(png_target_path)
    xml_set = set()
    png_set = set()
    # check that number found is equal to the number expected
    for zip_fname in zip_fnames:
        with ZipFile(os.path.join(zip_path, zip_fname)) as zipfile:
            for name in zipfile.namelist():
                if name.endswith('xml'):
                    xml_set.add(name.split('/')[-1])
                elif name.endswith('png'):
                    png_set.add(name.split('/')[-1])
    assert len(xml_set) == expected
    # unpack XML and PNG files
    xml_set_copy = set(xml_set)
    png_set_copy = set(png_set)
    for zip_fname in zip_fnames:
        logger.info('processing {}'.format(zip_fname))
        with ZipFile(os.path.join(zip_path, zip_fname)) as zipfile:
            for name in zipfile.namelist():
                if name.endswith('xml'):
                    xml_fname = name.split('/')[-1]
                    if xml_fname in xml_set_copy:
                        xml_set_copy.remove(xml_fname)
                        target_fname = os.path.join(xml_target_path, xml_fname)
                        with zipfile.open(name) as fid, open(target_fname, 'w') as out_fid:
                            out_fid.buffer.write(fid.read())
                elif name.endswith('png'):
                    png_fname = name.split('/')[-1]
                    if png_fname in png_set_copy:
                        png_set_copy.remove(png_fname)
                        target_fname = os.path.join(png_target_path, png_fname)
                        with zipfile.open(name) as fid, open(target_fname, 'w') as out_fid:
                            out_fid.buffer.write(fid.read())
    assert len(xml_set_copy) == len(png_set_copy) == 0
    # index the repository
    return initialize(xml_target_path)


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Prepare EMTF material from SPUDS.',
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            epilog=EPILOG)
    parser.add_argument('destination_path',
                        type=str,
                        help='root path to store XML and PNG files')
    parser.add_argument('zip_path',
                        type=str,
                        help='root path where SPUDS EMTF .zip files are to be found')
    parser.add_argument('expected',
                        type=int,
                        help='number of EMTF records expected across the .zip files')
    args = parser.parse_args(argv[1:])

    mirror(args.destination_path,
           args.zip_path,
           args.expected)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
