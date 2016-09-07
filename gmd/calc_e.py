from __future__ import division

import logging
import sys
import os
import math
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import numpy as NP
from scipy.constants import mu_0
from scipy.io import savemat

from conductivity import surface_impedance_1D
from usgs_conductivity import USGS_MODEL_MAP
from ..util.signal import nextpow2
from ..mag.iaga2002 import parse
from ..util.nan import nan_interp
from ..util.date import toJ2000


def calc_e(Bx, By, Zw_function, interval):
    """
    Calculate the electric field induced by the magnetic field (*Bx*,
    eastward oriented magnetic field intensity, and *By*, northward
    oriented magnetic field intensity) and surface impedance function
    *Zw_function* (input is angular frequency [rad] and output is
    impedance [Ohms]). Return the tuple *Ex*, *Ey* (same orientation
    as the magnetic field).

    Units:
       (input)  Bx, By: [T]
       (input)  Zw_function: [rad] -> [Ohm]
       (output) Ex, Ey: [V/m]
    """
    # This function is the reimplementation of Prof. Zhu's matlab code
    # which is based on the algorithm detailed in the NERC Application
    # Guide, Dec. 2013. Note that the NERC Guide contains a sign error
    # that has been corrected in the code below.

    # function global variables
    lp = 60  # length at beginning and end of time series to fit
             # detrend parameters
    p  = 60  # length at beginning and end of time series for
             # frequency domain window
    # get the FFT size
    assert len(Bx) == len(By)
    N = len(Bx)
    Nfft = nextpow2(N)
    # zero-mean and remove linear trend from the magnetic time series
    # data
    ax = NP.mean(Bx[:lp])
    bx = NP.mean(Bx[-lp-1:])
    N_range = NP.arange(1, N+1)
    cx = ax * (N - N_range) / N + bx * N_range / N
    Bx0 = Bx - cx

    ay = NP.mean(By[:lp])
    by = NP.mean(By[-lp-1:])
    cy = ay * (N - N_range) / N + by * N_range / N
    By0 = By - cy
    # window the magnetic data time series
    wp = NP.hstack((0.5 * (1 - NP.cos(2 * math.pi * NP.arange(0, p/2) / p)),
                    NP.ones(N - p),
                    0.5 * (1 - NP.cos(2 * math.pi * NP.arange(p/2 - 1, -1, -1) / p))))
    Bx1 = Bx0 * wp
    By1 = By0 * wp
    # compute FFT
    Sbx = NP.fft.fft(Bx1, n=Nfft) / N
    Sby = NP.fft.fft(By1, n=Nfft) / N
    freq = NP.arange(0, Nfft) / Nfft / interval
    omega = 2 * math.pi * freq

    Zw_positive = Zw_function(omega[1:])
    Zw = NP.hstack(([0], Zw_positive))

    Zw2 = NP.hstack((Zw[:(int(Nfft/2)+1)],
                     NP.conj(Zw[1:int(Nfft/2)])[::-1]))

    Se_x =  Zw2 * Sby / mu_0
    Se_y = -Zw2 * Sbx / mu_0

    Ex = NP.real(NP.fft.ifft(Se_x, Nfft) * N)
    Ey = NP.real(NP.fft.ifft(Se_y, Nfft) * N)

    return (Ex[:N],
            Ey[:N])


def apply_transfer_function(Bx,
                            By,
                            interval,
                            model,
                            model_map=USGS_MODEL_MAP):
    """
    Filter *Bx* and *By* (in [T]) with 3-D transfer function where
    *interval* is the sample period (in [s]). Uses the USGS 1-D
    conductivity model with key *model*. The model information are
    stored in the *model* key map *model_map*. Ex and Ey (in [V/m]).
    """
    # setup surface impedance function
    usgs_model = model_map[model]
    Zw_function = lambda omega: surface_impedance_1D(usgs_model, omega)
    # calculate E field
    Ex, Ey = calc_e(nan_interp(Bx),
                    nan_interp(By),
                    Zw_function,
                    interval)
    return Ex, Ey


def process(output_mat_fname,
            input_iaga2002_fname,
            model,
            model_map=USGS_MODEL_MAP,
            save_B=False):
    """
    End-to-end processing of an IAGA2002 magnetometer data record
    *input_iaga2002_fname* to the output file *output_mat_fname*
    containing the calculated E-field (units are [V/m]). Use the 1-D
    USGS conductivity model with ID *model* identifying the
    conductivity model to use in *conductivity_map* (keys are the IDs
    and values are 1-D conductivity models in the USGS format).
    """
    # gather Bx and By magnetometer measurements
    _, data_map = parse(iaga2002_fname)
    interval = int((data_map.keys()[1] - data_map.keys()[0]).total_seconds())
    Bx = nan_interp([record.x * 1e-9 for record in data_map.itervalues()])
    By = nan_interp([record.y * 1e-9 for record in data_map.itervalues()])
    # filter with transfer function
    Ex, Ey = apply_transfer_function(Bx,
                                     By,
                                     interval,
                                     model,
                                     model_map=model_map)
    # save E field
    stn_name = os.path.basename(input_iaga2002_fname)[:3]
    j2000 = map(toJ2000, data_map.iterkeys())
    mdict = {'Ex': Ex,
             'Ey': Ey,
             'j2000': j2000,
             'input_fname': os.path.abspath(input_iaga2002_fname),
             'model': model}
    if save_B:
        mdict['Bx'] = Bx
        mdict['By'] = By
        mdict['Bx_raw'] = [getattr(x, stn_name.upper() + 'X') * 1e-9 for x in data_map.itervalues()]
        mdict['By_raw'] = [getattr(x, stn_name.upper() + 'Y') * 1e-9 for x in data_map.itervalues()]
    savemat(output_mat_fname, mdict)
    return output_mat_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Compute E-field from B-field using USGS 1-D model.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_mat_fname',
                        type=str,
                        help='output, containing modeled E-field, in .mat format')
    parser.add_argument('input_iaga2002_fname',
                        type=str,
                        help='input IAGA2002 magnetometer data file')
    parser.add_argument('model',
                        type=str,
                        choices=sorted(USGS_MODEL_MAP),
                        help='process use the given 1-D conductivity model')
    args = parser.parse_args(argv[1:])

    process(args.output_mat_fname,
            args.input_iaga2002_fname,
            args.model)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
