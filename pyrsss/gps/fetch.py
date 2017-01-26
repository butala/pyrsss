import logging
import os
from ftplib import FTP
from contextlib import closing
from collections import defaultdict

from ..util.path import touch_path, decompress


logger = logging.getLogger('pyrsss.gps.fetch')


TEMPLATE_MAP = {'CORS': {'obs': ('geodesy.noaa.gov',
                                 '/cors/rinex/{date:%Y}/{date:%j}/{stn}/{stn}{date:%j}0.{date:%y}',
                                 'o.gz')}}


def fetch(source,
          dates,
          stns,
          rinex_type='obs',
          template_map=TEMPLATE_MAP,
          local_path='./',
          local_template='{stn}{date:%j}0.{date:%y}{suffix}'):
    """
    ???
    """
    server, template, suffix = TEMPLATE_MAP[source][rinex_type]
    fname_map = defaultdict(dict)
    logger.info('opening connection to {}'.format(server))
    with closing(FTP(server)) as ftp:
        ftp.login()
        for date in dates:
            for stn in stns:
                remote_fname = template.format(date=date, stn=stn) + suffix
                local_fname = os.path.join(local_path.format(date=date, stn=stn, suffix=suffix),
                                           local_template.format(date=date, stn=stn, suffix=suffix))
                logger.info('fetching {} and storing to {}'.format(remote_fname,
                                                                   local_fname))
                touch_path(os.path.dirname(local_fname))
                with open(local_fname, 'w') as fid:
                    try:
                        ftp.retrbinary('RETR {}'.format(remote_fname),
                                       fid.write)
                        fname_map[date][stn] = local_fname
                    except Exception as e:
                        logger.warning('could not fetch {} ({}) --- skipping'.format(remote_fname,
                                                                                     e))
                        os.remove(local_fname)
                        continue
    for date in sorted(fname_map):
        for stn in sorted(fname_map[date]):
            fname_map[date][stn] = decompress(fname_map[date][stn])
    return fname_map
