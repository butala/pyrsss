from datetime import datetime

import sh


def output_clk_file(clk_fname,
                    Ts,
                    first_sample = datetime(2000, 1, 1),
                    epoch = datetime(2000, 1, 1)):
    """
    Ts: sample period (in seconds)
    """
    with open(clk_fname, 'w') as fid:
        fid.write(f'{Ts:d}.\n')
        fid.write(f'{first_sample:%y %m %d %H %M %S}\n')
        fid.write(f'{epoch:%y %m %d %H %M %S}\n')
    return clk_fname


def output_asc_file(asc_fname, Bx, By, Bz, Ex, Ey):
    """Units?

    Important: Bx, By, Bz, Ex, and Ey must be integer valued. The
    internal data structure to store this information in emtf is an
    integer array.
    """
    with open(asc_fname, 'w') as fid:
        for vals in zip(Bx, By, Bz, Ex, Ey):
            line = ('{:7d}' * len(vals)).format(*vals)
            fid.write(line + '\n')
    return asc_fname


def rfasc(bin_fname,
          asc_fname,
          clk_fname,
          station_id='NONE',
          header='',
          bin='~/src/emtf/EMTF/RF/RFASC/rfasc'):
    """
    """
    cmd = sh.Command(bin)
    input_str = f'{station_id}\n{bin_fname}\n{header}\n{asc_fname}\n{clk_fname}\nn\nn\n'
    output = cmd(_in=input_str)
    return bin_fname
