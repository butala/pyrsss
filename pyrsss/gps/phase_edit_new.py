import logging
from collections import namedtuple

from phase_edit import phase_edit
from rinex_new import week_sec2dt, RinexDump


""" ??? """
class ArcInfo(namedtuple('ArcInfo', 'gap tot sat ok s start end dt obs_types')):
    pass


def parse_discfix_log(log_fname):
    """
    """
    phase_breaks = []
    with open(log_fname) as fid:
        for line in fid:
            if line.startswith('Fine'):
                toks = line.split()[2:]
                info = ArcInfo(int(toks[0]),
                               int(toks[1]),
                               toks[2],
                               int(toks[3]),
                               int(toks[4]),
                               week_sec2dt(int(toks[5]),
                                           float(toks[6])),
                               week_sec2dt(int(toks[7]),
                                           float(toks[8])),
                               float(toks[9]),
                               ' '.join(toks[10:]))
                phase_breaks.append(info)
    return phase_breaks


def label_phase_arcs(rinex_dump, phase_breaks):
    """
    """
    rinex_dump.loc[:, 'arc'] = -1
    for i, phase_break in enumerate(phase_breaks):
        I = (rinex_dump.sat == phase_break.sat) & \
            (rinex_dump.gps_time >= phase_break.start) & \
            (rinex_dump.gps_time <= phase_break.end)
        rinex_dump.loc[I, 'arc'] = i
    # why are there data not assigned to an arc?
    I = rinex_dump[rinex_dump.arc == -1].index
    # rinex_dump.drop(I, inplace=True)
    return rinex_dump


def apply_rejections(rinex_dump, time_reject_map):
    """
    """
    # total = 0
    for sat, rejections in time_reject_map.iteritems():
        # count = 0
        for rejection in rejections:
            I = (rinex_dump.sat == sat) & \
                (rinex_dump.gps_time >= rejection.lower) & \
                (rinex_dump.gps_time <= rejection.upper)
            rinex_dump.drop(rinex_dump[I].index, inplace=True)
            # count += sum(I)
            # total += sum(I)
        # print(sat, count)
    # print(total)
    return rinex_dump


def apply_phase_adjustments(rinex_dump, phase_adjust_map):
    """
    """
    for sat, adjustments in phase_adjust_map.iteritems():
        for dt, col, offset in adjustments:
            I = (rinex_dump.sat == sat) & \
                (rinex_dump.gps_time >= dt)
            # print(rinex_dump.loc[I, col].iloc[0])
            rinex_dump.loc[I, col] -= offset
            # print(rinex_dump.loc[I, col].iloc[0])
    return rinex_dump


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    normalized_rinex_fname = '/tmp/jplm0010.14o'

    work_path = '/tmp/work'
    discfix_args = {}

    (time_reject_map,
     phase_adjust_map) = phase_edit(normalized_rinex_fname,
                                    work_path=work_path,
                                    discfix_args=discfix_args)

    # print(time_reject_map.items()[0])
    # print(phase_adjust_map.items()[0])

    # print(time_reject_map)
    # print(phase_adjust_map)

    log_fname = '/tmp/work/jplm0010.14o.df.log'

    phase_breaks = parse_discfix_log(log_fname)


    dump_fname = '/tmp/jplm0010.14o.dump'

    rinex_dump = RinexDump.load(dump_fname)

    count = 0
    for p in phase_breaks:
        count += p.tot

    ok_count = 0
    for p in phase_breaks:
        ok_count += p.ok

    print(rinex_dump.shape)
    print(ok_count)
    print(count)
    print(count + 64)

    print(phase_breaks[0])

    from collections import Counter

    counts = Counter()
    total = 0
    for p in phase_breaks:
        counts[p.sat] += p.tot - p.ok
        total += p.tot - p.ok

    print(counts)
    print(total)

    label_phase_arcs(rinex_dump, phase_breaks)

    df = rinex_dump[rinex_dump.arc == -1]

    # I  = rinex_dump.P1.isnull()
    # I |= rinex_dump.P2.isnull()
    I  = rinex_dump.P2.isnull()
    # I |= rinex_dump.L1.isnull()
    I |= rinex_dump.L2.isnull()
    # I |= rinex_dump.C1.isnull()

    print(df.shape)

    # print(df)

    df2 = rinex_dump.loc[I]
    print(df2.shape)

    print(sum(rinex_dump.P2.isnull()))
    print(sum(rinex_dump.L2.isnull()))

    df3 = rinex_dump[rinex_dump.arc != -1]

    print(sum(df3.C1.isnull()))
    print(sum(df3.P1.isnull()))
    print(sum(df3.L1.isnull()))

    print(sum(df3.P2.isnull()))
    print(sum(df3.L2.isnull()))


    # assert False

    apply_phase_adjustments(rinex_dump, phase_adjust_map)

    # assert False

    print(rinex_dump.shape)
    apply_rejections(rinex_dump, time_reject_map)
    print(rinex_dump.shape)

    # phase_adjust(rinex_dump, phase_adjust_map)

    # print(rinex_dump.iloc[:10])
    print(rinex_dump[rinex_dump.arc == -1].shape)


    print(rinex_dump.shape[0] - rinex_dump[rinex_dump.arc == -1].shape[0])

    pkl_fname = '/tmp/jplm0010.14o.pkl'
    rinex_dump.to_pickle(pkl_fname)
