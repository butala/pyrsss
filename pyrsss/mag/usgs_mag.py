import os

import geomagio
from obspy.core import UTCDateTime
from geomagio.algorithm import XYZAlgorithm
import numpy as np


def fetch_usgs(out_path,
               date,
               stn,
               type='variation',
               interval='second',
               he=False,
               out_template='{stn}{date:%Y%m%d}v{suffix}.{suffix}'):
    """
    Fetch USGS magnetometer data for *date* at *stn*. Store the data
    in IAGA2002 format and return the file name (*out_template* serves
    as a template). Limit to data *type* and *interval*. If *he*, then
    include the H and E channels in the output (that is, local
    magnetic north and east components).
    """
    out_fname = os.path.join(out_path,
                             out_template.format(stn=stn.lower(),
                                                 date=date,
                                                 suffix=interval[:3]))

    input_factory = geomagio.edge.EdgeFactory()
    timeseries = input_factory.get_timeseries(
        observatory = stn,
        channels = ('H', 'E', 'Z', 'F'),
        type = type,
        interval = interval,
        starttime = UTCDateTime('{date:%Y-%m-%d}T00:00:00Z'.format(date=date)),
        endtime = UTCDateTime('{date:%Y-%m-%d}T23:59:59Z'.format(date=date)))

    if all([np.isnan(trace).all() for trace in timeseries.traces]):
        raise ValueError('no data for {} on {:%Y-%m-%d} found'.format(stn, date))

    # convert from HEZF channels to XYZF channels
    algorithm = XYZAlgorithm(informat='obs', outformat='geo')
    xyzf = algorithm.process(timeseries)

    with open(out_fname, 'w') as fid:
        output_factory = geomagio.iaga2002.IAGA2002Factory()
        output_factory.write_file(
            channels = ('H', 'E', 'X', 'Y', 'Z', 'F') if he else ('X', 'Y', 'Z', 'F'),
            fh = fid,
            timeseries = xyzf)
    return out_fname
