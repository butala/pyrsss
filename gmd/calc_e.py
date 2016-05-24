from __future__ import division

import logging
import sys
import os
import math
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from cStringIO import StringIO

import numpy as NP
import scipy.constants
from scipy.io import savemat

from pyrsss.util.signal import nextpow2
from pyrsss.gmd.conductivity import (NAME_MAP,
                                     parse_conductivity,
                                     surface_impedance_1D)
from pyrsss.mag.iaga2002 import parse
from pyrsss.util.nan import nan_interp
from pyrsss.util.date import toJ2000


def calc_e(Bx, By, Zw_function, interval):
    """
    Calculate the electric field induced by the magnetic field (*Bx*,
    eastward oriented magnetic field intensity, and *By*, northward
    oriented magnetic field intensity) and surface impedance function
    *Zw_function* (input is angular frequency [rad] and output is
    impedance [Ohms]). Return the tuple *Ex*, *Ey* (same orientation
    as the magnetic field).
    """
    # This function is the reimplementation of Prof. Zhu's matlab code
    # which is based on the algorithm detailed in the NERC Application
    # Guide, Dec. 2013. Note that the NERC Guide contains a sign error
    # that has been corrected in the code below.

    # function global variables
    lp = 60 * int(60 / interval) # length at beginning and end of time
                                 # series to fit detrend parameters
    p  = 60 * int(60 / interval) # length at beginning and end of time
                                 # series for frequency domain window
    mu0 = scipy.constants.mu_0
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

    Se_x =  Zw2 * Sby / mu0
    Se_y = -Zw2 * Sbx / mu0

    Ex = NP.real(NP.fft.ifft(Se_x, Nfft) * N)
    Ey = NP.real(NP.fft.ifft(Se_y, Nfft) * N)

    return (Ex[:N],
            Ey[:N])


def process(output_mat_fname,
            input_iaga2002_fname,
            model,
            conductivity_map=NAME_MAP):
    """
    End-to-end processing of an IAGA2002 magnetometer data record
    *input_iaga2002_fname* to the output file *output_mat_fname*
    containing the calculated E-field. Use the 1-D USGS conductivity
    model with ID *model* identifying the conductivity model to use in
    *conductivity_map* (keys are the IDs and values are 1-D
    conductivity models in the USGS format).
    """
    # gather Bx and By magnetometer measurements
    _, data_map = parse(input_iaga2002_fname)
    stn_name = os.path.basename(input_iaga2002_fname)[:3]
    interval = int((data_map.keys()[1] - data_map.keys()[0]).total_seconds())
    Bx = nan_interp([getattr(x, stn_name.upper() + 'X') * 1e-9 for x in data_map.itervalues()])
    By = nan_interp([getattr(x, stn_name.upper() + 'Y') * 1e-9 for x in data_map.itervalues()])
    # setup surface impedance function
    fid = StringIO(conductivity_map[model])
    usgs_model = parse_conductivity(fid)
    Zw_function = lambda omega: surface_impedance_1D(usgs_model, omega)
    # calculate E field
    Ex, Ey = calc_e(nan_interp(Bx),
                    nan_interp(By),
                    Zw_function,
                    interval)
    # save E field
    j2000 = map(toJ2000, data_map.iterkeys())
    savemat(output_mat_fname,
            {'Ex': Ex,
             'Ey': Ey,
             'j2000': j2000,
             'input_fname': os.path.abspath(input_iaga2002_fname),
             'model': model})
    return output_mat_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_mat_fname',
                        type=str,
                        help='')
    parser.add_argument('input_iaga2002_fname',
                        type=str,
                        help='input IAGA2002 magnetometer data file')
    parser.add_argument('model',
                        type=str,
                        choices=sorted(NAME_MAP),
                        help='process use the given 1-D conductivity model')
    args = parser.parse_args(argv[1:])

    process(args.output_mat_fname,
            args.input_iaga2002_fname,
            args.model)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
