import sys
import os
import logging

import sh

from ..util.path import SmartTempDir, replace_path

logger = logging.getLogger('pyrsss.gps.teqc')


def rinex_info(rinex_fname,
               nav_fname,
               work_path=None):
    """
    Query RINEX file *rinex_fname* and RINEX nav file *nav_fname* for
    useful information and return in a key/value mapping. Store
    intermediate files in *work_path* (a temporary, automatically
    cleaned up area if not specified).
    """
    # information mapping
    info = {}
    def process_output(line):
        if line.startswith('Receiver type'):
            info['receiver'] = line.split(':')[1].split('(')[0].strip()
        elif line.lstrip().startswith('antenna WGS 84 (xyz)'):
            # make sure units are [m]
            assert line.rstrip().endswith('(m)')
            info['xyz'] = map(float, line.split(':')[1].split('(')[0].split())
        elif line.lstrip().startswith('antenna WGS 84 (geo)'):
            if line.split(':')[1].lstrip()[0] in ['N', 'S']:
                # skip arcmin, arcsec line
                pass
            else:
                lat, _, lon, _ = line.split(':')[1].split(None, 3)
                info['lat'] = float(lat)
                lon = float(lon)
                while lon > 180:
                    lon -= 360
                info['lon'] = lon
        elif line.lstrip().startswith('WGS 84 height'):
            assert line.rstrip().endswith('m')
            info['height'] = float(line.split(':')[1].rstrip()[:-1])
        elif line.startswith('|qc - header| position'):
            # make sure units are [m]
            assert line.rstrip()[-1] == 'm'
            info['xyz error'] = float(line.split(':')[1].rstrip()[:-1])
        elif line.startswith('Observation interval'):
            info['interval'] = float(line.split(':')[1].split()[0])
        elif line.startswith('Moving average MP12'):
            info['MP12'] = float(line.split(':')[1].rstrip()[:-1])
        elif line.startswith('Moving average MP21'):
            info['MP21'] = float(line.split(':')[1].rstrip()[:-1])
    # query the RINEX file via teqc quality check --- process in given
    # work area to avoid intermediate file pollution
    with SmartTempDir(work_path) as work_path:
        intermediate_rinex_fname = replace_path(work_path, rinex_fname)
        os.symlink(os.path.abspath(rinex_fname),
                   intermediate_rinex_fname)
        intermediate_nav_fname = replace_path(work_path, nav_fname)
        os.symlink(os.path.abspath(nav_fname),
                   intermediate_nav_fname)
        sh.teqc('+qc',
                '+quiet',
                '-R',
                '-S',
                '-E',
                '-C',
                '-J',
                '-nav', intermediate_nav_fname,
                intermediate_rinex_fname,
                _cwd=work_path,
                _out=process_output,
                _err=sys.stderr)
        os.remove(intermediate_rinex_fname)
        os.remove(intermediate_nav_fname)
    return info


def rinex_merge(output_fname, rinex_fnames, _err=sys.stderr):
    """
    Using teqc, merge *rinex_fnames* and store to the file
    *output_fname*. Returns *output_fname*. Redirect error output to
    *_err*.
    """
    args = ['-pch'] + rinex_fnames
    sh.teqc(*args,
            _out=output_fname,
            _err=_err)
    return output_fname


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.getLogger('sh').setLevel(logging.WARNING)

    rinex_fname = '/Users/butala/src/absolute_tec/jplm0010.14o'
    nav_fname = '/Users/butala/src/absolute_tec/jplm0010.14n'

    info = rinex_info(rinex_fname,
                      nav_fname)

    for key in sorted(info):
        print('{:10s}: {}'.format(key, info[key]))
