import os
import json
import logging

import pandas as PD
from netCDF4 import Dataset

from ..util.date import fromJ2000

logger = logging.getLogger('pyrsss.mag.repository')


ROOT = None
"""
Data repository root path.
"""


def load_config(fname=os.path.expanduser('~/etc/intermagnet.json')):
    """
    Load repository information from *fname* and return the
    configuration mapping.
    """
    global ROOT
    if not os.path.isfile(fname):
        raise RuntimeError('could not fined repository configuration file {}'.format(fname))
    logger.info('loading repository information from {}'.format(fname))
    with open(fname) as fid:
        config = json.load(fid)
    ROOT = config['root']
    return config


def check_root():
    """
    Check that the global variable *ROOT* has been specified.
    """
    global ROOT
    if ROOT is None:
        raise RuntimeError('repository root path not specified (remember to call load_config)')
    return ROOT


TEMPLATE_MAP = {'definitive': '{root}/intermagnet/{date:%Y}/{stn}/IAGA2002/{stn}{date:%Y%m%d}d.min',
                0:            '{root}/intermagnet_level0/{year}/{stn}.{year}.nc',
                1:            '{root}/intermagnet_level1/{year}/{stn}.{year}.nc'}
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


def repository_years(level):
    """
    Given a repository data product *level*, return a sorted list of
    integer years found within the repository.
    """
    check_root()
    base_path = os.path.join(ROOT, 'intermagnet_level{}'.format(level))
    year_list = []
    for x in os.listdir(base_path):
        if os.path.isdir(os.path.join(base_path, x)) and len(x) == 4:
            year_list.append(x)
    return sorted(year_list)


def get_data_frame(stn,
                   level,
                   columns=slice(None),
                   years='all'):
    """
    Return a pandas data frame and mapping header information (one
    entry per year) for all available data for a given *stn* and data
    product *level*. Use *columns* to select columns (via a list of
    column names). Return data only for those years specified in the
    list *years* (the default is to return all available data).
    """
    check_root()
    frames = []
    info_map = {}
    if years == 'all':
        years = repository_years(level)
    for year in years:
        nc_fname = TEMPLATE_MAP[level].format(root=ROOT,
                                              year=year,
                                              stn=stn)
        logger.info('loading {}'.format(nc_fname))
        dataframe_i, info_i = nc_to_dataframe(nc_fname,
                                              columns=columns)
        frames.append(dataframe_i)
        info_map[year] = info_i
    return (PD.concat(frames),
            info_map)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    load_config()
    # combine all netCDF records for a given station
    dataframe, info_map = get_data_frame('frd', 0, columns=['time',
                                                            'H'])
    dataframe.info()
