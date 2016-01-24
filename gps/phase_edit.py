import math
from collections import namedtuple

from ..util.stats import Stats

"""
Called/used by process.py.

Implement cycle-slip detection and repair: Liu (and perhaps also
Blewitt).
"""

class WidelaneInfo(namedtuple('WidelaneInfo', 'arc_id N mean sigma')):
    pass


# GATHER CONFIG INTO CLASS

MIN_EL = 10

DELTA = 4

OUTLIER_DELTA = 10

N_OUTLIERS = 2

GAP_LENGTH = 5  # in [minutes]


def widelane_edit(obs_time_series,
                  min_el=MIN_EL,
                  delta=DELTA,
                  n_outliers=N_OUTLIERS):
    # NEED TIME --- REVISE!
    obs_list = list(obs_time_series.values())
    stats = Stats()
    index = 0
    phase_arc_id = 0
    outlier_count = 0
    info_list = []
    while index < len(obs_list):
        # SANITY CHECK --- REMOVE WHEN DEBUGGED
        assert len(info_list) == index

        obs_i = obs_list[index]
        if obs_i.el < min_el:
            info_list.append(WidelaneInfo(-1, None, None, None))
            index += 1
        elif obs_i.L1 == 0 or obs_i.L2 == 0:
            info_list.append(WidelaneInfo(-2, None, None, None))
            index += 1
        else:
            if stats.N > 2 and math.fabs(obs_i.N_WL - stats.mean) > delta * stats.sigma:
                outlier_count += 1
            else:
                outlier_count = 0
            if outlier_count == n_outliers:
                stats.reset()
                phase_arc_id += 1
                # BACKTRACK
                index -= outlier_count - 1
                if index < 0:
                    # PREVENT INFINITE LOOP
                    raise RuntimeError('CORNER CASE')
                if outlier_count > 1:
                    info_list = info_list[:-(outlier_count - 1)]
            else:
                stats(obs_i.N_WL)
                info_list.append(WidelaneInfo(phase_arc_id,
                                              stats.N,
                                              stats.mean,
                                              stats.sigma))
                index += 1
    return info_list


if __name__ == '__main__':
    from rinex import read_rindump
    # obs_map = read_rindump('/Users/butala/src/absolute_tec/jplm0010.14o.dump')
    obs_map = read_rindump('/Users/butala/src/absolute_tec/pie12600.07o.dump')

    # JPL, G18: delete short arcs?
    # JPL, G19, G20, G22: new phase arc after time gap?
    # JPL, G20: missed cycle slip in WL detection?

    # PIE1, G02: sanity check?
    # PIE1, G03: verify cycle slip is repaired
    # PIE1, G05, G09, G12, G13(?), G14(?), G23(caught in rate of change?), G29(see G23) : clear outlier fouls up the statistics
    # PIE1, G16, G20: what do here (rapid sign change in WL)?
    # PIE1, G17, G26: new phase arc after time gap?
    # PIE1, G30: delete short isolated arcs?

    sat = 'G32'

    MIN_EL = 10

    PLOT_MAP = {0:  {'c':      'b',
                     'marker': 'o'},
                1:  {'c':      'r',
                     'marker': '^'},
                2:  {'c':      'g',
                     'marker': 's'},
                3:  {'c':      'darkorange',
                     'marker': '<'},
                4:  {'c':      'darkmagenta',
                     'marker': 'H'},
                5:  {'c':      'saddlebrown',
                     'marker': '>'}}

    from itertools import groupby

    def plot_phase_arcs(dts,
                        values,
                        arc_ids,
                        arc_marker=True,
                        **kwds):
        for arc_id, z in groupby(zip(dts,
                                     values,
                                     arc_ids),
                                 lambda x: x[-1]):
            if arc_id < 0:
                # arcs missing L1 or L2
                continue
            dts_i, values_i, _ = zip(*z)
            PL.plot_date(dts_i,
                         values_i,
                         mec='none',
                         c=PLOT_MAP[arc_id]['c'],
                         marker=PLOT_MAP[arc_id]['marker'] if arc_marker else 'None',
                         **kwds)

    widelane_edit_info = widelane_edit(obs_map[sat])

    # print([x.arc_id for x in widelane_edit_info[-30:]])

    arc_ids = [x.arc_id for x in widelane_edit_info]

    import pylab as PL


    # fig = PL.figure()
    # PL.plot_date(obs_map[sat].keys(),
    #              [x.el for x in obs_map[sat].itervalues()],
    #              c='k')
    # PL.axhline(MIN_EL)
    # PL.xlabel('GPS time')
    # PL.ylabel('Elevation [deg]')

    LI = [obs.L1m - obs.L2m for obs in obs_map[sat].itervalues()]

    BETA = 1.5
    fig = PL.figure(figsize=(BETA * 11, BETA * 8.5))

    ax = PL.subplot(311)
    plot_phase_arcs(obs_map[sat].keys(),
                    LI,
                    arc_ids)
    PL.xlabel('GPS time')
    PL.ylabel('LI = L1 - L2 [m]')

    # ax = PL.subplot(311)
    # plot_phase_arcs(obs_map[sat].keys(),
    #                 [obs.L1m for obs in obs_map[sat].itervalues()],
    #                 arc_ids)
    # PL.xlabel('GPS time')
    # PL.ylabel('L1 [m]')

    # ax = PL.subplot(311)
    # plot_phase_arcs(obs_map[sat].keys(),
    #                 [obs.L1m for obs in obs_map[sat].itervalues()],
    #                 arc_ids)
    # PL.xlabel('GPS time')
    # PL.ylabel('L1 [m]')


    PI = [obs.P2 - obs.P1 for obs in obs_map[sat].itervalues()]

    ax = PL.subplot(312)
    plot_phase_arcs(obs_map[sat].keys(),
                    PI,
                    arc_ids)
    PL.xlabel('GPS time')
    PL.ylabel('PI = P2 - P1 [m]')

    N_WL = [obs.N_WL for obs in obs_map[sat].itervalues()]

    N_WL_mean = [info.mean for info in widelane_edit_info]
    N_WL_sigma = [info.sigma for info in widelane_edit_info]

    sigma_plus = [mean_i + DELTA * sigma_i if sigma_i is not None else None for mean_i, sigma_i in zip(N_WL_mean, N_WL_sigma)]
    sigma_minus = [mean_i - DELTA * sigma_i if sigma_i is not None else None for mean_i, sigma_i in zip(N_WL_mean, N_WL_sigma)]

    ax = PL.subplot(313)
    plot_phase_arcs(obs_map[sat].keys(),
                    N_WL,
                    arc_ids)
    plot_phase_arcs(obs_map[sat].keys(),
                    N_WL_mean,
                    arc_ids,
                    arc_marker=False,
                    ls='-')
    plot_phase_arcs(obs_map[sat].keys(),
                    sigma_plus,
                    arc_ids,
                    arc_marker=False,
                    ls='-')
    plot_phase_arcs(obs_map[sat].keys(),
                    sigma_minus,
                    arc_ids,
                    arc_marker=False,
                    ls='-')
    PL.xlabel('GPS time')
    PL.ylabel('N_WL [m]')
    # PL.grid(axis='y')

    PL.suptitle('sat = {}'.format(sat))

    PL.show()
