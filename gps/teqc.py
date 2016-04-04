import sys
import logging

import sh

logger = logging.getLogger('pyrsss.gps.teqc')


def rinex_info(rinex_fname,
               nav_fname):
    """
    Query RINEX file *rinex_fname* and RINEX nav file *nav_fname* for
    useful information and return in a key/value mapping.
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
    # error handling
    error = ''
    def process_error(line):
        error += line
    # query the RINEX file via teqc quality check
    sh.teqc('+qc',
            '+quiet',
            '-R',
            '-S',
            '-E',
            '-C',
            '-J',
            '-nav', nav_fname,
            rinex_fname,
            _out=process_output,
            _err=process_error)
    if error:
        raise RuntimeError('failure running teqc on '
                           '{} nav={} ({})'.format(rinex_fname,
                                                   nav_fname))
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
