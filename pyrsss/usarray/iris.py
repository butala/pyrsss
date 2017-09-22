import logging
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import pandas as PD

from resp import get_station_resp
from ..util.date import dt_parser, UNIX_EPOCH

logger = logging.getLogger('pyrsss.usarray.iris')


def fetch(stn, dt1, dt2, location=0, resp=None):
    """
    Request USArray MT data from IRIS (http://ds.iris.edu/ds) at
    station *stn* starting at UTC time *dt1* and ending at
    *dt2*. Return as a :class:`DataFrame`.
    """
    d1 = UTCDateTime((dt1 - UNIX_EPOCH).total_seconds())
    d2 = UTCDateTime((dt2 - UNIX_EPOCH).total_seconds())
    client = Client('IRIS')
    lfe = client.get_waveforms('EM', stn, location, 'LFE', d1, d2, longestonly=True)
    lfn = client.get_waveforms('EM', stn, location, 'LFN', d1, d2, longestonly=True)
    lfz = client.get_waveforms('EM', stn, location, 'LFZ', d1, d2, longestonly=True)
    lqe = client.get_waveforms('EM', stn, location, 'LQE', d1, d2, longestonly=True)
    lqn = client.get_waveforms('EM', stn, location, 'LQN', d1, d2, longestonly=True)
    # time sanity checks
    assert lfe.traces[0].meta.starttime == lfn.traces[0].meta.starttime \
        == lfe.traces[0].meta.starttime == lqe.traces[0].meta.starttime \
        == lqn.traces[0].meta.starttime
    assert (lfe.traces[0].times() == lfn.traces[0].times()).all()
    assert (lfe.traces[0].times() == lfz.traces[0].times()).all()
    assert (lfe.traces[0].times() == lqe.traces[0].times()).all()
    assert (lfe.traces[0].times() == lqn.traces[0].times()).all()
    dt = [(lfe.traces[0].meta.starttime + x).datetime for x in lfe.traces[0].times()]
    # log information about data
    logger.info('time of first record = {}'.format(dt[0]))
    logger.info('time of last record = {}'.format(dt[-1]))
    logger.info('total data points = {}'.format(len(dt)))
    # get instrument response if needed
    if resp is None:
        try:
            resp = get_station_resp(stn, dt1)
        except:
            raise ValueError('could not find instrument response for {} on {:%Y-%m-%d}'.format(stn, dt1))
    if dt2 not in resp['interval']:
        raise NotImplementedError('date range {} -- {} spans multiple station response records'.format(dt1, dt2))
    assert resp['station'] == stn
    assert resp['network'] == 'EM'
    # build DataFrame and apply calibration values
    Bx = lfn.traces[0].data / resp['LFN']
    By = lfe.traces[0].data / resp['LFE']
    Bz = lfz.traces[0].data / resp['LFZ']
    Ex = lqn.traces[0].data / resp['LQN']
    Ey = lqe.traces[0].data / resp['LQE']
    return PD.DataFrame(index=dt, data={'B_X': Bx,
                                        'B_Y': By,
                                        'B_Z': Bz,
                                        'E_X': Ex,
                                        'E_Y': Ey})


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Fetch USArray MT data from IRIS.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_fname',
                        type=str,
                        help='name of output HDF5 file')
    parser.add_argument('station_id',
                        type=str,
                        help='USArray station identifier')
    parser.add_argument('dt1',
                        type=dt_parser,
                        help='start date and time (UTC)')
    parser.add_argument('dt2',
                        type=dt_parser,
                        help='end date and time (UTC)')
    args = parser.parse_args(argv[1:])

    df = fetch(args.station_id,
               args.dt1,
               args.dt2)

    df.to_hdf(args.output_fname, 'iris')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
