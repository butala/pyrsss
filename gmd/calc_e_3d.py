from __future__ import division

import logging
import sys
import os
import math
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import xml.etree.ElementTree as ET
from collections import OrderedDict

import numpy as NP
from scipy.constants import mu_0
from scipy.io import savemat
from scipy.interpolate import CubicSpline

from ..util.signal import nextpow2
from ..mag.iaga2002 import parse
from ..util.nan import nan_interp
from ..util.date import toJ2000


def calc_e_3d(Bx,
              By,
              Zxx_function,
              Zxy_function,
              Zyx_function,
              Zyy_function,
              interval):
    """
    Calculate the electric field induced by the magnetic field (*Bx*,
    eastward oriented magnetic field intensity, and *By*, northward
    oriented magnetic field intensity). *Zxx_function*,
    *Z_xy_function*, *Z_yx_function*, and *Z_yy_function* are the
    surface impedance functions (input is angular frequency [rad] and
    output is impedance [Ohms]). Return the tuple *Ex*, *Ey* (same
    orientation as the magnetic field).

    Units:
       (input)  Bx, By: [T]
       (input)  Zxx_function: [rad] -> [Ohm]
       (input)  Zxy_function: [rad] -> [Ohm]
       (input)  Zyx_function: [rad] -> [Ohm]
       (input)  Zyy_function: [rad] -> [Ohm]
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

    Zxx_positive = Zxx_function(omega[1:])
    Zxx = NP.hstack(([0], Zxx_positive))

    Zxx2 = NP.hstack((Zxx[:(int(Nfft/2)+1)],
                      NP.conj(Zxx[1:int(Nfft/2)])[::-1]))

    Zxy_positive = Zxy_function(omega[1:])
    Zxy = NP.hstack(([0], Zxy_positive))

    Zxy2 = NP.hstack((Zxy[:(int(Nfft/2)+1)],
                      NP.conj(Zxy[1:int(Nfft/2)])[::-1]))

    Zyx_positive = Zyx_function(omega[1:])
    Zyx = NP.hstack(([0], Zyx_positive))

    Zyx2 = NP.hstack((Zyx[:(int(Nfft/2)+1)],
                      NP.conj(Zyx[1:int(Nfft/2)])[::-1]))

    Zyy_positive = Zyy_function(omega[1:])
    Zyy = NP.hstack(([0], Zyy_positive))

    Zyy2 = NP.hstack((Zyy[:(int(Nfft/2)+1)],
                      NP.conj(Zyy[1:int(Nfft/2)])[::-1]))

    Se_x = []
    Se_y = []
    for Zxx_i, Zxy_i, Zyx_i, Zyy_i, Sbx_i, Sby_i in zip(Zxx2, Zxy2, Zyx2, Zyy2, Sbx, Sby):
        Z = NP.array([[Zxx_i, Zxy_i], [Zyx_i, Zyy_i]])
        Sb_i = NP.array([[Sbx_i / mu_0], [Sby_i / mu_0]])
        Se_i = NP.dot(Z, Sb_i)
        Se_x.append(Se_i[0, 0])
        Se_y.append(Se_i[1, 0])

    Ex = NP.real(NP.fft.ifft(Se_x, Nfft) * N)
    Ey = NP.real(NP.fft.ifft(Se_y, Nfft) * N)

    return (Ex[:N],
            Ey[:N])


def parse_xml(xml_fname):
    """
    Parse the E-M transfer function file *xml_fname* returning a
    mapping between period keys (in [s]) and the 2x2 matrices
    containing the associated Zxx, Zxy, Zyx, and Zyy parameters (in
    [mV / km] / [nT]).

    These files are available at http://ds.iris.edu/spud/emtf.
    """
    tree = ET.parse(xml_fname)
    root = tree.getroot()

    data_list = root.findall('Data')
    assert len(data_list) == 1
    data = data_list[0]

    Z_map = OrderedDict()
    for period in data.findall('Period'):
        Z_list = period.findall('Z')
        assert len(Z_list) == 1
        Z = Z_list[0]
        values = []
        for value, name in zip(Z, ['Zxx', 'Zxy', 'Zyx', 'Zyy']):
            assert value.attrib['name'] == name
            values.append(complex(*map(float, value.text.split())))
        Z_array = NP.array(values)
        Z_array.shape = 2, 2
        Z_map[float(period.attrib['value'])] = Z_array
    return Z_map


class Zw_interpolator(object):
    def __init__(self, Z_map):
        """
        Construct a cubic-spline 3-D E-M transfer function interpolater
        using the information in *Z_map* returned from
        :func:`parse_xml` as the function samples.
        """
        self.Z_map = Z_map
        periods = Z_map.keys()
        self.f = NP.array([1/x for x in periods[::-1]])
        self.omega = 2 * math.pi * self.f
        self.Zxx_interp = CubicSpline(self.omega, [x[0, 0] for x in Z_map.values()[::-1]])
        self.Zxy_interp = CubicSpline(self.omega, [x[0, 1] for x in Z_map.values()[::-1]])
        self.Zyx_interp = CubicSpline(self.omega, [x[1, 0] for x in Z_map.values()[::-1]])
        self.Zyy_interp = CubicSpline(self.omega, [x[1, 1] for x in Z_map.values()[::-1]])
        self.key_map = {'xx': self.Zxx_interp,
                        'xy': self.Zxy_interp,
                        'yx': self.Zyx_interp,
                        'yy': self.Zyy_interp}

    def __call__(self, omega, key):
        """
        Return the interpolated value of the transfer function (*key* is
        xx, xy, yx, or yy) at angular frequency *omega*. The units of
        *omega* and the output are the same as *Z_map* ([rad] and
        [mV/km]/[nT] by default).
        """
        return self.key_map[key](omega)


def apply_transfer_function(Bx, By, interval, xml_fname):
    """
    Filter *Bx* and *By* (in [T]) with 3-D transfer function where
    *interval* is the sample period (in [s]). Uses the 3-D transfer
    function model given in *xml_fname*. Return Ex and Ey (in [V/m]).
    """
    # setup surface impedance function
    Z_map = parse_xml(xml_fname)
    interp = Zw_interpolator(Z_map)
    # mu_0 * 1e3 converts from [mv / km] / [nT] to [Ohm]
    Zxx_function = lambda omega: interp(omega, 'xx') * mu_0 * 1e3
    Zxy_function = lambda omega: interp(omega, 'xy') * mu_0 * 1e3
    Zyx_function = lambda omega: interp(omega, 'yx') * mu_0 * 1e3
    Zyy_function = lambda omega: interp(omega, 'yy') * mu_0 * 1e3
    # calculate E field
    Ex, Ey = calc_e_3d(nan_interp(Bx),
                       nan_interp(By),
                       Zxx_function,
                       Zxy_function,
                       Zyx_function,
                       Zyy_function,
                       interval)
    return Ex, Ey


def process(output_mat_fname,
            input_iaga2002_fname,
            xml_fname,
            save_B=False):
    """
    End-to-end processing of an IAGA2002 magnetometer data record
    *input_iaga2002_fname* to the output file *output_mat_fname*
    containing the calculated E-field (units are [V/m]). Use the 3-D
    transfer function model given in *xml_fname.* Also save the
    B-field information if *save_B*.
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
                                     xml_fname)
    # save E field
    stn_name = os.path.basename(iaga2002_fname)[:3]
    j2000 = map(toJ2000, data_map.iterkeys())
    mdict = {'Ex': Ex,
             'Ey': Ey,
             'j2000': j2000,
             'input_fname': os.path.abspath(input_iaga2002_fname),
             'xml_fname': os.path.abspath(xml_fname)}
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

    parser = ArgumentParser('Compute E-field from B-field using 3-D transfer function model.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_mat_fname',
                        type=str,
                        help='')
    parser.add_argument('input_iaga2002_fname',
                        type=str,
                        help='input IAGA2002 magnetometer data file')
    parser.add_argument('xml_fname',
                        type=str,
                        help='EM transfer function XML file')
    args = parser.parse_args(argv[1:])

    process(args.output_mat_fname,
            args.input_iaga2002_fname,
            args.xml_fname)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    sys.exit(main())
