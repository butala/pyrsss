import os
import json
import logging

import pandas as PD
from netCDF4 import Dataset

from ..util.date import fromJ2000

logger = logging.getLogger('pyrsss.mag.repository')


CONFIG_JSON = os.path.expanduser('~/etc/intermagnet.json')
"""
Path to INTERMAGNET repository JSON configuration file.
"""


def get_config(fname=CONFIG_JSON):
    """
    Return the configuration mapping specified in the JSON file
    *fname*.
    """
    if not os.path.isfile(fname):
        raise RuntimeError('could not fined repository configuration file {}'.format(fname))
    logger.info('loading repository root from {}'.format(fname))
    with open(fname) as fid:
        try:
            config = json.load(fid)
        except:
            logger.error('could not parse JSON configuration file {}'.format(fname))
            raise
    return config


def get_root(fname=CONFIG_JSON):
    """
    Return the repository root path found in the JSON configuration
    file *fname*.
    """
    return get_config(fname=fname)['root']


def get_login(fname=CONFIG_JSON):
    """
    Return the a tuple containing the INTERMAGNET ftp login user name
    and password stored in the JSON configuration file *fname*.
    """
    config = get_config(fname=fname)
    return config['user'], config['password']


PATH_MAP = {'definitive': '{root}/intermagnet',
            0:            '{root}/intermagnet_level0',
            1:            '{root}/intermagnet_level1'}
"""
Mapping between processing levels and repository sub-directories.
"""


TEMPLATE_MAP = {'definitive': '{path}/{date:%Y}/{stn}/IAGA2002/{stn}{date:%Y%m%d}d.min',
                0:            '{path}/{year}/{stn}.{year}.nc',
                1:            '{path}/{year}/{stn}.{year}.nc'}
"""
Mapping between INTERMAGNET repository data types and file name
template.
"""


def nc_to_dataframe(nc_fname,
                    columns=slice(None)):
    """
    Return a pandas data frame containing the information in the
    netCDF file *nc_fname*. Return it and a mapping to the header
    metadata. Use *columns* to select columns (via a list of column
    names).
    """
    root = Dataset(nc_fname)
    data = {}
    data.update({dim: root[dim][:] for dim in root.dimensions if dim != 'time'})
    index = data['time'] = map(fromJ2000, root['time'][:])
    return (PD.DataFrame(data=data, index=index)[columns],
            {x: getattr(root, x) for x in root.ncattrs()})


def repository_years(root, level):
    """
    Given a repository at path *root* and data product *level*, return
    a sorted list of integer years found within the repository.
    """
    base_path = os.path.join(root, 'intermagnet_level{}'.format(level))
    year_list = []
    for x in os.listdir(base_path):
        if os.path.isdir(os.path.join(base_path, x)) and len(x) == 4:
            year_list.append(x)
    return sorted(year_list)


def get_data_frame(root,
                   stn,
                   level,
                   columns=slice(None),
                   years='all'):
    """
    Return a pandas data frame and mapping header information (one
    entry per year) for all available data for a given *stn* and data
    product *level* at repository at path *root*. Use *columns* to
    select columns (via a list of column names). Return data only for
    those years specified in the list *years* (the default is to
    return all available data).
    """
    frames = []
    info_map = {}
    if years == 'all':
        years = repository_years(root, level)
    for year in years:
        nc_fname = TEMPLATE_MAP[level].format(path=PATH_MAP[level].format(root=root),
                                              year=year,
                                              stn=stn)
        logger.info('loading {}'.format(nc_fname))
        dataframe_i, info_i = nc_to_dataframe(nc_fname,
                                              columns=columns)
        frames.append(dataframe_i)
        info_map[year] = info_i
    try:
        return (PD.concat(frames),
                info_map)
    except ValueError:
        raise RuntimeError('could find no level={} data for {}'.format(level,
                                                                       stn,
                                                                       root))


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    # combine all netCDF records for a given station
    root = get_root()
    dataframe, info_map = get_data_frame(root,
                                         'frd',
                                         0,
                                         columns=['time',
                                                  'H'],
                                         years=[1991, 1992])
    dataframe.info()
