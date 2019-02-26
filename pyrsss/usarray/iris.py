import logging
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import numpy as NP
import pandas as PD

from .resp import get_station_resp
from ..util.date import dt_parser, UNIX_EPOCH

logger = logging.getLogger('pyrsss.usarray.iris')


def fetch(stn, dt1, dt2, location=0, resp=None):
    """
    Request USArray MT data from IRIS (http://ds.iris.edu/ds) at
    station *stn* starting at UTC time *dt1* and ending at
    *dt2*. Return as a :class:`DataFrame`. Magnetic field components
    are stored with units nT. Electric field components are stored in
    mV / km.

    Note that USArray sensors are aligned to geomagnetic coordinates.
    """
    # get instrument response if needed
    if resp is None:
        try:
            resp = get_station_resp(stn, dt1, dt2)
        except:
            raise RuntimeError('could not find instrument response for {} in time range {:%Y-%m-%d %H:%M:%S} -- {:%Y-%m-%d %H:%M:%S}'.format(stn, dt1, dt2))
    d1 = UTCDateTime((dt1 - UNIX_EPOCH).total_seconds())
    d2 = UTCDateTime((dt2 - UNIX_EPOCH).total_seconds())
    client = Client('IRIS')
    lfe = client.get_waveforms('EM', stn, location, 'LFE', d1, d2)
    lfn = client.get_waveforms('EM', stn, location, 'LFN', d1, d2)
    lfz = client.get_waveforms('EM', stn, location, 'LFZ', d1, d2)
    lqe = client.get_waveforms('EM', stn, location, 'LQE', d1, d2)
    lqn = client.get_waveforms('EM', stn, location, 'LQN', d1, d2)
    assert len(lfe) == len(lfn) == len(lfz) == len(lqe) == len(lqn)
    dt_list = []
    Bn_list = []
    Be_list = []
    Bz_list = []
    En_list = []
    Ee_list = []
    for i in range(len(lfe)):
        # time sanity checks
        assert lfe.traces[i].meta.starttime == lfn.traces[i].meta.starttime \
            == lfe.traces[i].meta.starttime == lqe.traces[i].meta.starttime \
            == lqn.traces[i].meta.starttime
        assert (lfe.traces[i].times() == lfn.traces[i].times()).all()
        assert (lfe.traces[i].times() == lfz.traces[i].times()).all()
        assert (lfe.traces[i].times() == lqe.traces[i].times()).all()
        assert (lfe.traces[i].times() == lqn.traces[i].times()).all()
        dt = [(lfe.traces[i].meta.starttime + x).datetime for x in lfe.traces[i].times()]
        # log information about data
        logger.info('time of first record = {}'.format(dt[0]))
        logger.info('time of last record = {}'.format(dt[-1]))
        logger.info('total data points = {}'.format(len(dt)))
        # if dt[-1] not in resp['interval']:
            #raise NotImplementedError('date range {} -- {} spans multiple station response records'.format(dt[0], dt[-1]))
            # logger.warning('date range {} -- {} spans multiple station response records --- skipping'.format(dt[0], dt[-1]))
            # continue
        assert resp.stn == stn
        assert resp.network == 'EM'
        # apply calibration values and store data from trace
        dt_list.extend(dt)
        Bn_list.append(lfn.traces[i].data / resp.sensitivity(dt[0], 'LFN') * 1e9)
        Be_list.append(lfe.traces[i].data / resp.sensitivity(dt[0], 'LFE') * 1e9)
        Bz_list.append(lfz.traces[i].data / resp.sensitivity(dt[0], 'LFZ') * 1e9)
        En_list.append(lqn.traces[i].data / resp.sensitivity(dt[0], 'LQN') * 1e6)
        Ee_list.append(lqe.traces[i].data / resp.sensitivity(dt[0], 'LQE') * 1e6)
    return PD.DataFrame(index=dt_list,
                        data={'B_N': NP.hstack(Bn_list),
                              'B_E': NP.hstack(Be_list),
                              'B_Z': NP.hstack(Bz_list),
                              'E_N': NP.hstack(En_list),
                              'E_E': NP.hstack(Ee_list)})


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
