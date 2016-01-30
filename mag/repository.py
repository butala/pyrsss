import os
import logging

import pandas as PD
from netCDF4 import Dataset

from ..util.date import fromJ2000

logger = logging.getLogger('pyrsss.mag.repository')


TEMPLATE_MAP = {'definitive': '/rdata/airglow/butala/data/intermagnet/{date:%Y}/ams/IAGA2002/{stn}{date:%Y%m%d}d.min',
                0:            '/rdata/airglow/butala/data/intermagnet_level0/{year}/{stn}.{year}.nc',
                1:            '/rdata/airglow/butala/data/intermagnet_level1/{year}/{stn}.{year}.nc'}
"""
Mapping between INTERMAGNET repository data types and file name
template.
"""


def nc_to_dataframe(nc_fname):
    """
    Return a pandas data frame containing the information in the
    netCDF file *nc_fname*. Return it and a mapping to the header
    metadata.
    """
    root = Dataset(nc_fname)
    data = {}
    data.update({dim: root[dim][:] for dim in root.dimensions if dim != 'time'})
    index = data['time'] = map(fromJ2000, root['time'][:])
    return (PD.DataFrame(data=data, index=index),
            {x: getattr(root, x) for x in root.ncattrs()})


def repository_years(level):
    """
    Given a repository data product *level*, return a sorted list of
    integer years found within the repository. Assumes that the
    directory structure is /some/root/{year}/{files for this year}.nc.
    """
    base_path = os.path.split(os.path.split(TEMPLATE_MAP[level])[0])[0]
    year_list = []
    for x in os.listdir(base_path):
        if os.path.isdir(os.path.join(base_path, x)) and len(x) == 4:
            year_list.append(x)
    return sorted(year_list)


def complete_time_series(stn, level):
    """
    Return a pandas data frame and mapping header information (one
    entry per year) for all available data for a given *stn* and data
    product *level*.
    """
    frames = []
    info_map = {}
    for year in repository_years(level):
        nc_fname = TEMPLATE_MAP[level].format(year=year, stn=stn)
        logger.info('loading {}'.format(nc_fname))
        dataframe_i, info_i = nc_to_dataframe(nc_fname)
        frames.append(dataframe_i)
        info_map[year] = info_i
    return (PD.concat(frames),
            info_map)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    # combine all netCDF records for a given station
    dataframe, info_map = complete_time_series('frd', 0)
    dataframe.info()
