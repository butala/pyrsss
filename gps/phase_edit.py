import logging
import sys
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import defaultdict, OrderedDict
from datetime import datetime

import sh
from intervals import DateTimeInterval

from ..util.path import SmartTempDir, replace_path
from path import GPSTK_BUILD_PATH
from rinex import read_rindump, Observation, dump_rinex
from observation import ObsMap
from preprocess import normalize_rinex
from teqc import rinex_info

logger = logging.getLogger('pyrsss.gps.phase_edit')


"""
Called/used by process.py.

Implement cycle-slip detection and repair.
"""

DISC_FIX = os.path.join(GPSTK_BUILD_PATH,
                        'ext',
                        'apps',
                        'geomatics',
                        'cycleslips',
                        'DiscFix')
"""
???
"""

"""
ADD CONFIG CLASS
"""

"""
GPSTk Discontinuity Corrector (GDC) v.5.3 7/14/2008 configuration:
 DT=-1              : nominal timestep of data (seconds) [required - no default!]
 Debug=0            : level of diagnostic output to log, from none(0) to extreme(7)
 GFVariation=16     : expected maximum variation in GF phase in time DT (meters)
 MaxGap=180         : maximum allowed time gap within a segment (seconds)
 MinPts=13          : minimum number of good points in phase segment ()
 OutputDeletes=1    : if non-zero, include delete commands in the output cmd list
 OutputGPSTime=0    : if 0: Y,M,D,H,M,S  else: W,SoW (GPS) in editing commands
 WLSigma=1.5        : expected WL sigma (WL cycle) [NB = ~0.83*p-range noise(m)]
 useCA=0            : use C/A code pseudorange (C1) ()
For DiscFix, GDC commands are of the form --DC<GDCcmd>, e.g. --DCWLSigma=1.5
"""

"""
# Data configuration
 --decimate <dt>     Decimate data to time interval (sec) dt
 --forceCA           Use C/A code range, NOT P code (default: only if P absent)
 --gap <t>           Minimum data gap (sec) separating satellite passes (600)
 --onlySat <sat>     Process only satellite <sat> (a GPS SatID, e.g. G21)
 --exSat <sat>       Exclude satellite(s) [e.g. --exSat G22]
"""


def phase_edit(rinex_fname,
               interval,
               work_path=None,
               disc_fix=DISC_FIX):
    """
    ???
    """
    logger.info('applying GPSTk DiscFix to {} (interval={})'.format(rinex_fname,
                                                                    interval))
    command = sh.Command(disc_fix)
    with SmartTempDir(work_path) as work_path:
        log_fname = os.path.join(work_path, 'df.log')
        stdout_fname = os.path.join(work_path, 'df.stdout')
        stderr_fname = os.path.join(work_path, 'df.stderr')
        cmd_fname = os.path.join(work_path, 'df.out')
        command('--inputfile', rinex_fname,
                '--dt', str(interval),
                '--logOut', log_fname,
                '--cmdOut', cmd_fname,
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


def filter_obs_map(obs_map,
                   time_reject_map,
                   phase_adjust_map):
    """
    ???
    """
    edited_obs_map = ObsMap()
    for sat in sorted(obs_map):
        # add C1_delta, P1_delta, P2_delta: cc2noncc happens here
        L1_delta = 0
        L2_delta = 0
        reject_list = list(time_reject_map[sat])
        offset_list = list(phase_adjust_map[sat])
        for dt, obs in obs_map[sat].iteritems():
            # time rejection
            while reject_list and dt > reject_list[0].upper:
                reject_list.pop(0)
            if reject_list and dt in reject_list[0]:
                # delete observation from stream
                continue
            # phase adjustment
            while offset_list and dt >= offset_list[0][0]:
                _, obs_type, offset = offset_list.pop(0)
                if obs_type == 'L1':
                    L1_delta += offset
                elif obs_type == 'L2':
                    L2_delta += offset
                else:
                    # impossible
                    assert False
            edited_obs_map[sat][dt] = [obs.C1,
                                       obs.P1,
                                       obs.P2,
                                       None if obs.L1 is None else obs.L1 + L1_delta,
                                       None if obs.L2 is None else obs.L2 + L2_delta,
                                       obs.az,
                                       obs.el,
                                       obs.satx,
                                       obs.saty,
                                       obs.satz]
    return edited_obs_map


def phase_edit_process(h5_fname,
                       rinex_fname,
                       nav_fname,
                       work_path=None,
                       preprocess=True):
    """ ??? """
    info = rinex_info(rinex_fname,
                      nav_fname)
    interval = info['interval']
    with SmartTempDir(work_path) as work_path:
        # preprocess
        if preprocess:
            logger.info('preprocessing {}'.format(rinex_fname))
            unprocessed_rinex = rinex_fname
            rinex_fname = replace_path(work_path, rinex_fname)
            normalize_rinex(rinex_fname,
                            unprocessed_rinex)
        # phase edit
        logger.info('phase edit {}'.format(rinex_fname))
        (time_reject_map,
         phase_adjust_map) = phase_edit(rinex_fname,
                                        interval,
                                        work_path=work_path)
        # dump RINEX and read in ObsMap
        logger.info('dumping {}'.format(rinex_fname))
        rinex_dump_fname = replace_path(work_path, rinex_fname + '.dump')
        dump_rinex(rinex_dump_fname,
                   rinex_fname,
                   nav_fname)
        obs_map = read_rindump(rinex_dump_fname)
        # apply phase edit adjustments to ObsMap
        logger.info('applying phase edit adjustments')
        edited_obs_map = filter_obs_map(obs_map,
                                        time_reject_map,
                                        phase_adjust_map)
        # store ObsMap to file
        logger.info('storing output to {}'.format(h5_fname))
        edited_obs_map.dump(h5_fname, title='pyrsss.gps.phase_edit output')
    return h5_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Preprocess and apply the GPSTk phase editor '
                            '(DiscFix) to an input RINEX record and produce '
                            'a dump record suitable for subsequent processing '
                            '(phase leveling).',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('dump_fname',
                        type=str,
                        help='output dump file')
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
    # DECIMATE OPTION?
    args = parser.parse_args(argv[1:])

    phase_edit_process(args.dump_fname,
                       args.rinex_fname,
                       args.nav_fname,
                       work_path=args.work_path,
                       preprocess=not args.no_preprocess)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)
    sys.exit(main())
