import logging
import sys
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import defaultdict, namedtuple, OrderedDict
from datetime import datetime

import sh

import pandas as pd
from intervals import DateTimeInterval

from ..util.path import SmartTempDir, replace_path
from .path import get_gnsstk_build_path
from .constants import week_sec2dt
from .rinex import RinexDump, dump_rinex
from .preprocess import normalize_rinex

logger = logging.getLogger('pyrsss.gps.phase_edit')


"""
Called/used by process.py.

Implement cycle-slip detection and repair.
"""

DISC_FIX_RELPATH = os.path.join('ext',
                                'apps',
                                'geomatics',
                                'cycleslips',
                                'DiscFix')
"""Path to the DiscFix tool relative to the GNSSTk build directory."""


def get_disc_fix():
    """Return the full path to the DiscFix tool."""
    return os.path.join(get_gnsstk_build_path(), DISC_FIX_RELPATH)

"""
ADD CONFIG CLASS
"""

"""
GPSTk Discontinuity Corrector (GDC) v.6.3 12/15/2015 configuration:
 DT=-1              : nominal timestep of data (seconds) [required - no default!]
 Debug=0            : level of diagnostic output to log, from 0(none) to 7(extreme)
 GFVariation=16     : expected maximum variation in GF phase in time DT (meters)
 MaxGap=180         : maximum allowed time gap within a segment (seconds)
 MinPts=13          : minimum number of good points in phase segment ()
 OutputDeletes=1    : if non-zero, include delete commands in the output cmd list
 OutputGPSTime=0    : if 0, output Y,M,D,H,M,S else: W,SoW in edit cmds (log uses SatPass fmt)
 ResetUnique=0      : if non-zero, reset the unique number to zero
 WLSigma=1.5        : expected WL sigma (WL cycle) [NB = ~0.83*p-range noise(m)]
 useCA1=0           : use L1 C/A code pseudorange (C1) ()
 useCA2=0           : use L2 C/A code pseudorange (C2) ()

For DiscFix, GDC commands are of the form --DC<GDCcmd>, e.g. --DCWLSigma=1.5
"""

"""
# Data config:
 --decimate <dt>    Decimate data to time interval (sec) dt (0.00)
 --gap <t>          Minimum gap (sec) between passes [same as --DCMaxGap] (600) (600.00)
 --noCA1            Fail if L1 P-code is missing, even if L1 CA-code is present (don't)
 --noCA2            Fail if L2 P-code is missing, even if L2 CA-code is present (don't)
 --forceCA1         Use C/A L1 range, even if L1 P-code is present (don't)
 --forceCA2         Use C/A L2 range, even if L2 P-code is present (don't)
 --onlySat <sat>    Process only satellite <sat> (a SatID, e.g. G21 or R17) ()
 --exSat <sat>      Exclude satellite(s) [e.g. --exSat G22,R] [repeat] ()
 --doGLO            Process GLONASS satellites as well as GPS (don't)
 --GLOfreq <sat:n>  GLO channel #s for each sat [e.g. R17:-4] [repeat] ()
"""


def phase_edit(rinex_fname,
               work_path=None,
               disc_fix=None,
               discfix_args=[],
               glonass=False):
    """
    Apply GNSSTk DiscFix cycle-slip detection and repair to *rinex_fname*.
    Return the (time_reject_map, phase_adjust_map) parsed from the
    generated edit commands.
    """
    if disc_fix is None:
        disc_fix = get_disc_fix()
    logger.info('applying GNSSTk DiscFix to {}'.format(rinex_fname))
    command = sh.Command(disc_fix)
    with SmartTempDir(work_path) as work_path:
        basename = os.path.basename(rinex_fname)
        log_fname = os.path.join(work_path, basename + '.df.log')
        stdout_fname = os.path.join(work_path, basename + '.df.stdout')
        stderr_fname = os.path.join(work_path, basename + '.df.stderr')
        cmd_fname = os.path.join(work_path, basename + '.df.out')
        args = ['--obs', rinex_fname,
                '--log', log_fname,
                '--cmd', cmd_fname]
        if glonass:
            args += ['--doGLO']
        if discfix_args:
            logger.info('passing options to DiscFix: {}'.format(' '.join(discfix_args)))
            args += discfix_args
        command(*args,
                _out=stdout_fname,
                _err=stderr_fname)
        return parse_edit_commands(cmd_fname)


"""
# ------ Editing commands ------
# RINEX header modifications (arguments with whitespace must be quoted)
 --HDp <p>         Set header 'PROGRAM' field to <p> ()
 --HDr <rb>        Set header 'RUN BY' field to <rb> ()
 --HDo <obs>       Set header 'OBSERVER' field to <obs> ()
 --HDa <a>         Set header 'AGENCY' field to <a> ()
 --HDx <x,y,z>     Set header 'POSITION' field to <x,y,z> (ECEF, m) ()
 --HDm <m>         Set header 'MARKER NAME' field to <m> ()
 --HDn <n>         Set header 'MARKER NUMBER' field to <n> ()
 --HDj <n>         Set header 'REC #' field to <n> ()
 --HDk <t>         Set header 'REC TYPE' field to <t> ()
 --HDl <v>         Set header 'REC VERS' field to <v> ()
 --HDs <n>         Set header 'ANT #' field to <n> ()
 --HDt <t>         Set header 'ANT TYPE' field to <t> ()
 --HDh <h,e,n>     Set header 'ANTENNA OFFSET' field to <h,e,n> (Ht,East,North) ()
 --HDc <c>         Add 'COMMENT' <c> to the output header [repeat] ()
 --HDdc            Delete all comments [not --HDc] from input header (don't)
 --HDda            Delete all auxiliary header data (don't)
# Time related [t,f are strings, time t conforms to format f; cf. gpstk::Epoch.]
# Default t(f) is 'week,sec-of-week'(%F,%g) OR 'y,m,d,h,m,s'(%Y,%m,%d,%H,%M,%S)
 --OF <f,t>        At RINEX time <t>, close output file and open another named <f> ()
 --TB <t[:f]>      Start time: Reject data before this time ([Beginning of dataset])
 --TE <t[:f]>      Stop  time: Reject data after this time ([End of dataset])
 --TT <dt>         Tolerance in comparing times, in seconds (0.00)
 --TN <dt>         If dt>0, decimate data to times = TB + N*dt [sec, w/in tol] (0.00)
# In the following <SV> is a RINEX satellite identifier, e.g. G17 R7 E22 R etc.
#              and <OT> is a 3- or 4-char RINEX observation code e.g. C1C GL2X S2N
# Delete cmds; for start(stop) cmds. stop(start) time defaults to end(begin) of data
#     and 'deleting' data for a single OT means it is set to zero - as RINEX requires.
 --DA <t>          Delete all data at a single time <t> [repeat] ()
 --DA+ <t>         Delete all data beginning at time <t> [repeat] ()
 --DA- <t>         Stop deleting at time <t> [repeat] ()
 --DO <OT>         Delete RINEX obs type <OT> entirely (incl. header) [repeat] ()
 --DS <SV>         Delete all data for satellite <SV> [SV may be char]
 --DS <SV,t>       Delete all data for satellite <SV> at single time <t> [repeat] ()
 --DS+ <SV,t>      Delete data for satellite <SV> beginning at time <t> [repeat] ()
 --DS- <SV,t>      Stop deleting data for sat <SV> beginning at time <t> [repeat] ()
 --DD <SV,OT,t>    Delete a single RINEX datum(SV,OT) at time <t> [repeat] ()
 --DD+ <SV,OT,t>   Delete all RINEX data(SV,OT) starting at time <t> [repeat] ()
 --DD- <SV,OT,t>   Stop deleting RINEX data(SV,OT) at time <t> [repeat] ()
 --SD <SV,OT,t,d>  Set data(SV,OT) to value <d> at single time <t> [repeat] ()
 --SS <SV,OT,t,s>  Set SSI(SV,OT) to value <s> at single time <t> [repeat] ()
 --SL <SV,OT,t,l>  Set LLI(SV,OT) to value <l> at single time <t> [repeat] ()
 --SL+ <SV,OT,t,l> Set all LLI(SV,OT) to value <l> starting at time <t> [repeat] ()
 --SL- <SV,OT,t,l> Stop setting LLI(SV,OT) to value <l> at time <t> [repeat] ()
# Bias cmds: (BD cmds apply only when data is non-zero, unless --BZ)
 --BZ              Apply BD command even when data is zero (i.e. 'missing') (don't)
 --BS <SV,OT,t,s>  Add the value <s> to SSI(SV,OT) at single time <t> [repeat] ()
 --BL <SV,OT,t,l>  Add the value <l> to LLI(SV,OT) at single time <t> [repeat] ()
 --BD <SV,OT,t,d>  Add the value <d> to data(SV,OT) at single time <t> [repeat] ()
 --BD+ <SV,OT,t,d> Add the value <d> to data(SV,OT) beginning at time <t> [repeat] ()
 --BD- <SV,OT,t,d> Stop adding the value <d> to data(SV,OT) at time <t> [repeat] ()
"""


"""
-DSG14,2014,1,1,0,26,30.000000
-DS+G22,2014,1,1,0,11,30.000000 # begin delete of 11 points
-DS-G22,2014,1,1,0,15,30.000000 # end delete of 11 points
-BD+G22,L1,2014,1,1,0,16,0.000000,0 # WL
-BD+G22,L2,2014,1,1,0,16,0.000000,23 # WL
-DSG31,2014,1,1,0,4,0.000000
-DSG31,2014,1,1,0,24,30.000000
"""


def remove_comment(command):
    """
    Return the contents of *command* appearing before #.
    """
    return command.split('#')[0].strip()


def parse_date_fields(fields):
    """
    ???
    """
    assert len(fields) == 6
    date_str = '{}-{}-{} {}:{}:{}'.format(fields[0],
                                          fields[1].zfill(2),
                                          fields[2].zfill(2),
                                          fields[3].zfill(2),
                                          fields[4].zfill(2),
                                          fields[5])
    return datetime.strptime(date_str,
                             '%Y-%m-%d %H:%M:%S.%f')


def parse_delete_command(command):
    """
    ???
    """
    original_command = command
    command = remove_comment(command)
    if command.startswith('-DS+') or command.startswith('-DS-'):
        prefix = command[:4]
        sv_time = command[4:]
    elif command.startswith('-DS'):
        prefix = command[:3]
        sv_time = command[3:]
    else:
        raise RuntimeError('unrecognized delete command '
                           '{}'.format(original_command))
    fields = sv_time.split(',')
    if len(fields) != 7:
        raise RuntimeError('unrecognized fields found in command '
                           '{}'.format(original_command))
    sv = fields[0]
    dt = parse_date_fields(fields[1:])
    return prefix, sv, dt


def parse_bias_command(command):
    """
    ???
    """
    original_command = command
    command = remove_comment(command)
    if command.startswith('-BD+'):
        prefix = command[:4]
        fields = command[4:].split(',')
        if len(fields) != 9:
            raise RuntimeError('unrecognized fields found in command '
                               '{}'.format(original_command))
        sv = fields[0]
        obs_type = fields[1]
        dt = parse_date_fields(fields[2:8])
        offset = float(fields[8])
    else:
        raise RuntimeError('unrecognized bias command '
                           '{}'.format(original_command))
    return prefix, sv, obs_type, dt, offset


def parse_edit_commands(df_fname):
    """
    ???
    """
    time_reject_map = defaultdict(list)
    phase_adjust_map = defaultdict(list)
    start_command = None
    with open(df_fname) as fid:
        for line in fid:
            if line.startswith('-DS+'):
                if start_command is not None:
                    raise RuntimeError('adjacent start time ranges detected in '
                                       '{} ({})'.format(df_fname,
                                                        line))
                start_command = line
            elif line.startswith('-DS-'):
                if start_command is None:
                    raise RuntimeError('found time range end prior to start in '
                                       '{} ({})'.format(df_fname,
                                                        line))
                _, sv_start, dt_start = parse_delete_command(start_command)
                _, sv_end, dt_end = parse_delete_command(line)
                if sv_start != sv_end:
                    raise RuntimeError('time range start is for {} but end is '
                                       'for {}'.format(sv_start, sv_end))
                time_reject_map[sv_start].append(DateTimeInterval([dt_start,
                                                                   dt_end]))
                start_command = None
            elif line.startswith('-DS'):
                _, sv, dt = parse_delete_command(line)
                time_reject_map[sv].append(DateTimeInterval([dt, dt]))
            elif line.startswith('-BD+'):
                _, sv, obs_type, dt, offset = parse_bias_command(line)
                phase_adjust_map[sv].append((dt, obs_type, offset))
            else:
                raise NotImplementedError('unhandled edit command {}'.format(line))
    return time_reject_map, phase_adjust_map


# MAJOR REWRITE --- I OUTPUT AN ArcMap!!!

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)
    sys.exit(main())


class ArcInfo(namedtuple('ArcInfo', 'gap tot sat ok s start end dt obs_types')):
    """
    Record from the DiscFix "Fine" arc summary table: gap and tot point
    counts, satellite, ok count, solution flag, start/stop times, length
    [s], and observation types.
    """
    pass


def parse_discfix_log(log_fname):
    """
    Parse the DiscFix arc summary ("Fine" lines) from *log_fname*. Return
    the list of :class:`ArcInfo`.
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
    """
    Label each observation of *rinex_dump* with the DiscFix arc index
    given *phase_breaks* (from :func:`parse_discfix_log`). Observations
    outside all arcs are dropped (they cannot be leveled).
    """
    rinex_dump.loc[:, 'arc'] = -1
    for i, phase_break in enumerate(phase_breaks):
        I = (rinex_dump.sat == phase_break.sat) & \
            (rinex_dump.gps_time >= phase_break.start) & \
            (rinex_dump.gps_time <= phase_break.end)
        rinex_dump.loc[I, 'arc'] = i
    unassigned = rinex_dump.arc == -1
    if unassigned.any():
        logger.warning('dropping {} observations outside DiscFix '
                       'arcs'.format(int(unassigned.sum())))
        rinex_dump.drop(rinex_dump.index[unassigned], inplace=True)
    return rinex_dump


def apply_rejections(rinex_dump, time_reject_map):
    """
    Drop the rejected time intervals (sat -> list of intervals) from
    *rinex_dump*.
    """
    # total = 0
    for sat, rejections in time_reject_map.items():
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
    Apply phase clock offset adjustments (sat -> list of (dt, column,
    offset)) to *rinex_dump*.
    """
    for sat, adjustments in phase_adjust_map.items():
        for dt, col, offset in adjustments:
            I = (rinex_dump.sat == sat) & \
                (rinex_dump.gps_time >= dt)
            # print(rinex_dump.loc[I, col].iloc[0])
            rinex_dump.loc[I, col] -= offset
            # print(rinex_dump.loc[I, col].iloc[0])
    return rinex_dump


def phase_edit_rinex(rinex_fname,
                     nav_fname,
                     work_path=None,
                     discfix_args=[],
                     glonass=False,
                     preprocess=True,
                     p1c1=True):
    """
    Run the phase edit front end on *rinex_fname*: normalize the RINEX
    file, apply DiscFix, dump the observables, and return the edited
    :class:`RinexDump` with labeled arcs (see :func:`label_phase_arcs`).
    """
    with SmartTempDir(work_path) as work_path:
        if preprocess:
            normalized_rinex_fname = replace_path(work_path, rinex_fname)
            normalize_rinex(normalized_rinex_fname,
                            rinex_fname)
        else:
            normalized_rinex_fname = rinex_fname
        (time_reject_map,
         phase_adjust_map) = phase_edit(normalized_rinex_fname,
                                        work_path=work_path,
                                        discfix_args=discfix_args,
                                        glonass=glonass)
        log_fname = os.path.join(work_path,
                                 os.path.basename(normalized_rinex_fname) + '.df.log')
        dump_fname = replace_path(work_path, normalized_rinex_fname + '.dump')
        dump_rinex(dump_fname,
                   normalized_rinex_fname,
                   nav_fname)
        rinex_dump = RinexDump.load(dump_fname, p1c1=p1c1)
        phase_breaks = parse_discfix_log(log_fname)
        label_phase_arcs(rinex_dump, phase_breaks)
        apply_phase_adjustments(rinex_dump, phase_adjust_map)
        apply_rejections(rinex_dump, time_reject_map)
        return rinex_dump


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Preprocess and apply the GNSSTk phase editor '
                            '(DiscFix) to an input RINEX record and produce '
                            'an edited, arc-labeled dump record suitable for '
                            'subsequent processing (phase leveling).',
                            formatter_class=ArgumentDefaultsHelpFormatter,
                            epilog='Unrecognized arguments are passed on to DiscFix. See the DiscFix usage message for accepted options.')
    parser.add_argument('output_fname',
                        type=str,
                        help='output pickle file containing the edited RinexDump')
    parser.add_argument('rinex_fname',
                        type=str,
                        help='input RINEX observation file')
    parser.add_argument('nav_fname',
                        type=str,
                        help='input RINEX navigation file')
    parser.add_argument('--work-path',
                        '-w',
                        type=str,
                        default=None,
                        help='path to store intermediate files (use an '
                             'automatically cleaned up area if not specified)')
    parser.add_argument('--no-preprocess',
                        action='store_true',
                        help='disable RINEX preprocess step (i.e., '
                             'normalization)')
    parser.add_argument('--glonass',
                        action='store_true',
                        help='enable GLONASS processing')
    args, discfix_args = parser.parse_known_args(argv[1:])

    rinex_dump = phase_edit_rinex(args.rinex_fname,
                                  args.nav_fname,
                                  work_path=args.work_path,
                                  discfix_args=discfix_args,
                                  glonass=args.glonass,
                                  preprocess=not args.no_preprocess)
    pd.to_pickle(rinex_dump, args.output_fname)
    return args.output_fname


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)
    sys.exit(main())
