import logging
import sys
import os
import math
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import xml.etree.ElementTree as ET
from collections import OrderedDict, namedtuple

import numpy as np
from scipy.constants import mu_0
from scipy.io import savemat
from scipy.interpolate import CubicSpline

from ..signal.spectrum import nextpow2
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
    ax = np.mean(Bx[:lp])
    bx = np.mean(Bx[-lp-1:])
    N_range = np.arange(1, N+1)
    cx = ax * (N - N_range) / N + bx * N_range / N
    Bx0 = Bx - cx

    ay = np.mean(By[:lp])
    by = np.mean(By[-lp-1:])
    cy = ay * (N - N_range) / N + by * N_range / N
    By0 = By - cy
    # window the magnetic data time series
    wp = np.hstack((0.5 * (1 - np.cos(2 * math.pi * np.arange(0, p/2) / p)),
                    np.ones(N - p),
                    0.5 * (1 - np.cos(2 * math.pi * np.arange(p/2 - 1, -1, -1) / p))))
    Bx1 = Bx0 * wp
    By1 = By0 * wp
    # compute FFT
    Sbx = sp.fft.fft(Bx1, n=Nfft) / N
    Sby = sp.fft.fft(By1, n=Nfft) / N
    freq = np.arange(0, Nfft) / Nfft / interval
    omega = 2 * math.pi * freq

    Zxx_positive = Zxx_function(omega[1:])
    Zxx = np.hstack(([0], Zxx_positive))

    Zxx2 = np.hstack((Zxx[:(int(Nfft/2)+1)],
                      np.conj(Zxx[1:int(Nfft/2)])[::-1]))

    Zxy_positive = Zxy_function(omega[1:])
    Zxy = np.hstack(([0], Zxy_positive))

    Zxy2 = np.hstack((Zxy[:(int(Nfft/2)+1)],
                      np.conj(Zxy[1:int(Nfft/2)])[::-1]))

    Zyx_positive = Zyx_function(omega[1:])
    Zyx = np.hstack(([0], Zyx_positive))

    Zyx2 = np.hstack((Zyx[:(int(Nfft/2)+1)],
                      np.conj(Zyx[1:int(Nfft/2)])[::-1]))

    Zyy_positive = Zyy_function(omega[1:])
    Zyy = np.hstack(([0], Zyy_positive))

    Zyy2 = np.hstack((Zyy[:(int(Nfft/2)+1)],
                      np.conj(Zyy[1:int(Nfft/2)])[::-1]))

    Se_x = []
    Se_y = []
    for Zxx_i, Zxy_i, Zyx_i, Zyy_i, Sbx_i, Sby_i in zip(Zxx2, Zxy2, Zyx2, Zyy2, Sbx, Sby):
        Z = np.array([[Zxx_i, Zxy_i], [Zyx_i, Zyy_i]])
        Sb_i = np.array([[Sbx_i / mu_0], [Sby_i / mu_0]])
        Se_i = np.dot(Z, Sb_i)
        Se_x.append(Se_i[0, 0])
        Se_y.append(Se_i[1, 0])

    Ex = np.real(sp.fft.ifft(Se_x, Nfft) * N)
    Ey = np.real(sp.fft.ifft(Se_y, Nfft) * N)

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
            if value.attrib['name'] != name and value.attrib['name'] != name.upper():
                raise ValueError('name mismatch ({} != {})'.format(value.attrib['name'], name))
            values.append(complex(*map(float, value.text.split())))
        Z_array = np.array(values)
        Z_array.shape = 2, 2
        Z_map[float(period.attrib['value'])] = Z_array
    return Z_map


class XMLRecord(namedtuple('XMLRecord', ['T', 'Z', 'V', 'S', 'N'])):
    pass


def parse_xml_all(xml_fname):
    """
    Parse the E-M transfer function file *xml_fname* returning a
    record of response period (`T` in [s]), impedance (`Z`), variance
    (`V`), inverse coherent signal power (`S`), and residual
    covariance (`N`). The impedance has units [mV / km] / [nT] and
    matrices have the squared impedance units.

    This routine should, at some point, supersede :func:`parse_xml` .

    These files are available at http://ds.iris.edu/spud/emtf.
    """
    record = XMLRecord(*[[] for i in range(5)])

    tree = ET.parse(xml_fname)
    root = tree.getroot()

    data_list = root.findall('Data')
    assert len(data_list) == 1
    data = data_list[0]

    def parse_mat(tree, mat_key, is_complex=True, strict=True):
        Z_list = tree.findall(mat_key)
        assert len(Z_list) == 1
        Z = Z_list[0]
        assert len(Z) == 4
        values = []
        for value, name in zip(Z, ['Zxx', 'Zxy', 'Zyx', 'Zyy']):
            if strict and value.attrib['name'] != name and value.attrib['name'] != name.upper():
                raise ValueError('name mismatch ({} != {})'.format(value.attrib['name'], name))
            if is_complex:
                values.append(complex(*map(float, value.text.split())))
            else:
                values.append(float(value.text))
        Z_array = np.array(values)
        Z_array.shape = 2, 2
        return Z_array

    for period in data.findall('Period'):
        record.T.append(float(period.attrib['value']))
        record.Z.append(parse_mat(period, 'Z'))
        record.V.append(parse_mat(period, 'Z.VAR', is_complex=False))
        record.S.append(parse_mat(period, 'Z.INVSIGCOV', strict=False))
        record.N.append(parse_mat(period, 'Z.RESIDCOV', strict=False))
    return record


def parse_xml_header(xml_fname):
    """
    Parse the E-M transfer function file *xml_fname* returning a
    mapping of site specific information, such as location and data
    quality.
    """
    header_map = {}

    tree = ET.parse(xml_fname)
    root = tree.getroot()

    copyright_list = root.findall('Copyright')
    assert len(copyright_list) == 1
    copyright_ = copyright_list[0]

    for copyright_node in copyright_:
        if copyright_node.tag == 'Acknowledgement':
            header_map['acknowledgement'] = copyright_node.text

    site_list = root.findall('Site')
    assert len(site_list) == 1
    site = site_list[0]

    for site_node in site:
        if site_node.tag == 'Id':
            header_map['id'] = site_node.text
        elif site_node.tag == 'Location':
            for child in site_node:
                if child.tag == 'Latitude':
                    header_map['lat'] = float(child.text)
                elif child.tag == 'Longitude':
                    header_map['lon'] = float(child.text)
                elif child.tag == 'Elevation':
                    assert child.attrib['units'] == 'meters'
                    header_map['el'] = float(child.text)
                elif child.tag == 'Declination':
                    header_map['dec'] = float(child.text)
                    header_map['dec_epoch'] = float(child.attrib['epoch'])
        elif site_node.tag == 'DataQualityNotes':
            for child in site_node:
                if child.tag == 'Rating':
                    header_map['rating'] = int(child.text)
                elif child.tag == 'GoodFromPeriod':
                    header_map['good_from'] = float(child.text)
                    try:
                        header_map['good_to_mHz'] = 1/float(child.text) * 1e3
                    except ZeroDivisionError:
                        header_map['good_to_mHz'] = float('inf')
                elif child.tag == 'GoodToPeriod':
                    header_map['good_to'] = float(child.text)
                    try:
                        header_map['good_from_mHz'] = 1/float(child.text) * 1e3
                    except ZeroDivisionError:
                        header_map['good_from_mHz'] = float('inf')
                elif child.tag == 'Comments':
                    header_map['data_quality_comments'] = child.text
        elif site_node.tag == 'DataQualityWarnings':
            for child in site_node:
                if child.tag == 'Flag':
                    header_map['data_quality_flag'] = int(child.text)
                elif child.tag == 'Comments':
                    header_map['data_quality_warning_comments'] = child.text
        elif site_node.tag == 'Acknowledgment':
            header_map['acknowledgment'] = site_node.text

    return header_map


class Zw_interpolator(object):
    def __init__(self, Z_map, extrapolate0=False):
        """
        Construct a cubic-spline 3-D E-M transfer function interpolater
        using the information in *Z_map* returned from
        :func:`parse_xml` as the function samples. If *extrapolate0*,
        then 0s are inserted in the transfer function response where
        extrapolation would occur (this happens when transfer function
        response is requested at frequencies outside the range
        provided in the .XML file record).
        """
        self.Z_map = Z_map
        periods = list(Z_map.keys())
        values = list(Z_map.values())
        self.f = np.array([1/x for x in periods[::-1]])
        self.omega = 2 * math.pi * self.f
        self.Zxx_interp = CubicSpline(self.omega, [x[0, 0] for x in values[::-1]],
                                      extrapolate=False)
        self.Zxy_interp = CubicSpline(self.omega, [x[0, 1] for x in values[::-1]],
                                      extrapolate=False)
        self.Zyx_interp = CubicSpline(self.omega, [x[1, 0] for x in values[::-1]],
                                      extrapolate=False)
        self.Zyy_interp = CubicSpline(self.omega, [x[1, 1] for x in values[::-1]],
                                      extrapolate=False)
        self.key_map = {'xx': self.Zxx_interp,
                        'xy': self.Zxy_interp,
                        'yx': self.Zyx_interp,
                        'yy': self.Zyy_interp}
        self.extrapolate0 = extrapolate0


    def __call__(self, omega, key):
        """
        Return the interpolated value of the transfer function (*key* is
        xx, xy, yx, or yy) at angular frequency *omega*. The units of
        *omega* and the output are the same as *Z_map* ([rad] and
        [mV/km]/[nT] by default).
        """
        if self.extrapolate0:
            y = np.zeros_like(omega, dtype=np.complex128)
            I = np.where((omega >= self.omega[0]) & (omega <= self.omega[-1]))
            y[I] = self.key_map[key](np.asarray(omega)[I])
            return y
        return self.key_map[key](omega)


def apply_transfer_function(Bx, By, interval, xml_fname, extrapolate0=False):
    """
    Filter *Bx* and *By* (in [T]) with 3-D transfer function where
    *interval* is the sample period (in [s]). Uses the 3-D transfer
    function model given in *xml_fname*. Return Ex and Ey (in
    [V/m]). If *extrapolate0*, use 0s in the transfer function
    response at frequencies outside the range provided in *xml_fname*.
    """
    # setup surface impedance function
    Z_map = parse_xml(xml_fname)
    interp = Zw_interpolator(Z_map, extrapolate0=extrapolate0)
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
    _, data_map = parse(input_iaga2002_fname)
    interval = int((data_map.keys()[1] - data_map.keys()[0]).total_seconds())
    Bx = nan_interp([record.x * 1e-9 for record in data_map.values()])
    By = nan_interp([record.y * 1e-9 for record in data_map.values()])
    # filter with transfer function
    Ex, Ey = apply_transfer_function(Bx,
                                     By,
                                     interval,
                                     xml_fname)
    # save E field
    stn_name = os.path.basename(input_iaga2002_fname)[:3]
    j2000 = map(toJ2000, data_map.keys())
    mdict = {'Ex': Ex,
             'Ey': Ey,
             'j2000': j2000,
             'input_fname': os.path.abspath(input_iaga2002_fname),
             'xml_fname': os.path.abspath(xml_fname)}
    if save_B:
        mdict['Bx'] = Bx
        mdict['By'] = By
        mdict['Bx_raw'] = [x.x * 1e-9 for x in data_map.values()]
        mdict['By_raw'] = [x.y * 1e-9 for x in data_map.values()]
    savemat(output_mat_fname, mdict)
    return output_mat_fname


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = ArgumentParser('Compute E-field from B-field using 3-D transfer function model.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('output_mat_fname',
                        type=str,
                        help='output, containing modeled E-field, in .mat format')
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
