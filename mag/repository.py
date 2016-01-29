import pandas as PD
from netCDF4 import Dataset

from ..util.date import fromJ2000


TEMPLATE_MAP = {'definitive': '/rdata/airglow/butala/data/intermagnet/{date:%Y}/ams/IAGA2002/{stn}{date:%Y%m%d}d.min',
                0:            '/rdata/airglow/butala/data/intermagnet_level0/{year}/{stn}.{year}.nc',
                1:            '/rdata/airglow/butala/data/intermagnet_level1/{year}/{stn}.{year}.nc'}


def nc_to_dataframe(nc_fname):
    root = Dataset(nc_fname)
    data = {}
    data.update({dim: root[dim][:] for dim in root.dimensions if dim != 'time'})
    index = data['time'] = map(fromJ2000, root['time'][:])
    return (PD.DataFrame(data=data, index=index),
            {x: getattr(root, x) for x in root.ncattrs()})


if __name__ == '__main__':
    dataframe, info = nc_to_dataframe(TEMPLATE_MAP[0].format(year=1991, stn='frd'))
