import logging
from collections import namedtuple

from level import LeveledArc
from util import shell_mapping
from ipp import cnv_azel2latlon
from teqc import rinex_info

"""
Instead, call this unbias or debias?

Remove transmitter and receiver biases from phase-connected arcs. Use
Attila's method to estimate the receiver bias (using IGS IONEX records
containing satellite biases and modeled VTEC maps).
"""


class AugLeveledArc(namedtuple('AugLeveledArc',
                               ' '.join(LeveledArc._fields + ('el_map',
                                                              'ipp_lat',
                                                              'ipp_lon')))):
    def __new__(cls,
                leveled_arc,
                site_lat,
                site_lon,
                height=450):
        fields = leveled_arc._asdict()
        fields['el_map'] = [shell_mapping(el_i, h=height) for el_i in fields['el']]
        ipp_lat, ipp_lon = zip(*[cnv_azel2latlon(az_i,
                                                 el_i,
                                                 (site_lat, site_lon),
                                                 ht=height) for (az_i, el_i) in zip(fields['az'],
                                                                                    fields['el'])])
        fields['ipp_lat'] = ipp_lat
        fields['ipp_lon'] = ipp_lon
        return cls._make(fields.values())


if __name__ == '__main__':
    from phase_edit import phase_edit, filter_obs_map
    from rinex import read_rindump, dump_rinex
    from level import level_phase_to_code

    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)

    rinex_fname = '/Users/butala/src/absolute_tec/jplm0010.14o'
    nav_fname = '/Users/butala/src/absolute_tec/jplm0010.14n'
    interval = 30

    (time_reject_map,
     phase_adjust_map) = phase_edit(rinex_fname, interval)

    rinex_dump_fname = '/tmp/jplm0010.14o.dump'
    dump_rinex(rinex_dump_fname,
               rinex_fname,
               nav_fname)

    obs_map = read_rindump(rinex_dump_fname)

    obs_map = filter_obs_map(obs_map,
                             time_reject_map,
                             phase_adjust_map)

    arc_map = level_phase_to_code(obs_map)

    info = rinex_info(rinex_fname,
                      nav_fname)

    leveled_arc = arc_map['G26'][0]
    aug_leveled_arc = AugLeveledArc(leveled_arc,
                                    info['lat'],
                                    info['lon'])

    import ipdb; ipdb.set_trace()
