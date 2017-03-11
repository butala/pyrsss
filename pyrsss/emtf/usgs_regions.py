import logging
import os
import shutil
from urllib2 import urlopen
from contextlib import closing

import sh
import cartopy.io.shapereader as shpreader
from shapely.geometry import Point

from ..util.path import SmartTempDir
from ..util.angle import convert_lon

logger = logging.getLogger('pyrsss.emtf.usgs_regions')


KMZ_URL = 'https://geomag.usgs.gov/conductivity/ConductivityRegions.kmz'
"""
URL where USGS publishes the conductivity region KMZ file.
"""


REGION_PATH = os.path.join(os.path.dirname(__file__),
                           'ConductivityRegions')
"""
Local path where processed conductivity region information is
cached.
"""


def update(path, kmz_url=KMZ_URL):
    """
    Update the local conductivity region cache at *path* using the KMZ
    file found at *kmz_url*.
    """
    with SmartTempDir() as tmp_path:
        # fetch KMZ file
        local_fname = os.path.join(tmp_path,
                                   kmz_url.split('/')[-1])
        logger.info('fetching {}'.format(kmz_url))
        with open(local_fname, 'w') as out_fid, closing(urlopen(kmz_url)) as in_fid:
            out_fid.write(in_fid.read())
        # convert KMZ file
        logger.info('converting KMZ file')
        sh.ogr2ogr('ConductivityRegions',
                   'ConductivityRegions.kmz',
                   '-nlt',
                   'POLYGON',
                   f='ESRI Shapefile',
                   _cwd=tmp_path)
        # save result
        dst = os.path.join(path, 'ConductivityRegions')
        shutil.copytree(os.path.join(tmp_path, 'ConductivityRegions'),
                        dst)
    assert os.path.isdir(dst)
    return dst


def initialize(region_path=REGION_PATH):
    """
    Initialize the local conductivity region cache if it does not
    already exist at *region_path*.
    """
    if os.path.isdir(region_path):
        logger.warning('{} exists --- skipping initialization'.format(region_path))
        return region_path
    return update(os.path.split(region_path)[0])


def get_region(lat, lon, region_path=REGION_PATH):
    """
    Return the 2 character, 1 number code (e.g., IP-3) associated with
    the physiographic region that encloses the point specified by
    *lat* and *lon*. Return `None` if the point is not interior to any
    region. Look for the required files at *region_path*.
    """
    try:
        regions = shpreader.Reader(os.path.join(region_path, 'ConductivityRegions'))
    except:
        initialize()
        regions = shpreader.Reader(os.path.join(region_path, 'ConductivityRegions'))
    geos = regions.geometries()
    names = [x.attributes['Name'] for x in regions.records()]
    point = Point(convert_lon(lon),
                  lat)
    encloses = [name for name, geo in zip(names, geos) if geo.contains(point)]
    if encloses:
        assert len(encloses) == 1
        return encloses[0]
    else:
        return None
