import os

import geomagio
from obspy.core import UTCDateTime
from geomagio.algorithm import XYZAlgorithm



def fetch_usgs(out_path,
               date,
               stn,
               type='variation',
               interval='second',
               out_template='{stn}{date:%Y%m%d}vsec.sec'):
    """
    Fetch USGS magnetometer data for *date* at *stn*. Store the data
    in IAGA2002 format and return the file name (*out_template* serves
    as a template). Limit to data *type* and *interval*.
    """
    out_fname = os.path.join(out_path,
                             out_template.format(stn=stn.lower(),
                                                 date=date))

    input_factory = geomagio.edge.EdgeFactory()
    timeseries = input_factory.get_timeseries(
        observatory = stn,
        channels = ('H', 'E', 'Z', 'F'),
        type = type,
        interval = interval,
        starttime = UTCDateTime('{date:%Y-%m-%d}T00:00:00Z'.format(date=date)),
        endtime = UTCDateTime('{date:%Y-%m-%d}T23:59:59Z'.format(date=date)))

    # convert from HEZF channels to XYZF channels
    algorithm = XYZAlgorithm(informat='obs', outformat='geo')
    xyzf = algorithm.process(timeseries)

    with open(out_fname, 'w') as fid:
        output_factory = geomagio.iaga2002.IAGA2002Factory()
        output_factory.write_file(
            channels = ('X', 'Y', 'Z', 'F'),
            fh = fid,
            timeseries = xyzf)
    return out_fname
