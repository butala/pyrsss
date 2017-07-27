import os
import logging
import cPickle
from glob import glob

from ..emtf.calc_e_3d import parse_xml_header
from ..util.distance import distance

logger = logging.getLogger('pyrsss.usarray_emtf.index')


PKL_FNAME = 'index.pkl'
"""
File name to associate the EMTF repository index.
"""


class EMTFIndex(dict):
    def __init__(self, repository_path):
        """
        Initialize the EMTF repository index based on the .xml files found
        at *repository_path*.
        """
        key = os.path.join(repository_path, '*.xml')
        logger.info('search for USArray 3-D EMTF files {}'.format(key))
        for fname in sorted(glob(key)):
            fname_key = '_'.join(os.path.basename(fname).split('.')[:-1])
            info = parse_xml_header(fname)
            self[fname_key] = (info, fname)
        logger.info('index a total of {} .xml file'.format(len(self)))

    def quality_subset(self, min_quality=5):
        """
        Return the :class:`EMTFIndex` corresponding to those EMTFs with a
        quality index of at least *min_quality*.
        """
        index = super(EMTFIndex, self).__new__(EMTFIndex)
        for key, (info, fname) in self.iteritems():
            if info['rating'] >= min_quality:
                index[key] = (info, fname)
        return index

    def by_distance(self, lat, lon, d_km):
        """
        Return the :class:`EMTFIndex` corresponding to those EMTFs within
        *d_km* (in km) of geodetic *lat* and *lon*.
        """
        usarray_lat, usarray_lon = [], []
        for key, (info, fname) in self.iteritems():
            usarray_lat.append(info['lat'])
            usarray_lon.append(info['lon'])
        lat_list = [lat] * len(usarray_lat)
        lon_list = [lon] * len(usarray_lon)
        d = distance(lat_list, lon_list,
                     usarray_lat, usarray_lon)
        index = super(EMTFIndex, self).__new__(EMTFIndex)
        for d_i, (k, v) in zip(d, self.iteritems()):
            if d_i <= d_km:
                index[k] = v
        return index


def update(repository_path, pkl_fname=PKL_FNAME):
    """
    Update the EMTF index file at *repository_path* named *pkl_fname*
    based on the .xml files found at *repository_path*. Return name of
    the pickle file storing the index.
    """
    index = EMTFIndex(repository_path)
    pkl_fullpath = os.path.join(repository_path,
                                pkl_fname)
    logger.info('updating index file {}'.format(pkl_fullpath))
    with open(pkl_fullpath, 'w') as fid:
        cPickle.dump(index, fid, -1)
    return pkl_fullpath


def initialize(repository_path, pkl_fname=PKL_FNAME):
    """
    Initialize the EMTF index file at *repository_path* named
    *pkl_fname* if it does not exist.
    """
    pkl_fullpath = os.path.join(repository_path,
                                pkl_fname)
    if os.path.isfile(pkl_fullpath):
        logger.info('{} exists --- not updating'.format(pkl_fullpath))
        did_update = False
    else:
        logger.info('{} does not exist --- initialize'.format(pkl_fullpath))
        update(repository_path, pkl_fname=pkl_fname)
        did_update = True
    return pkl_fullpath, did_update


def get_index(repository_path, pkl_fname=PKL_FNAME):
    """
    Return the index information for the EMTF repository located at
    *repository_path*.
    """
    pkl_fname, _ = initialize(repository_path, pkl_fname=pkl_fname)
    with open(pkl_fname) as fid:
        return cPickle.load(fid)
