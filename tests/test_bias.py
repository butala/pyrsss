"""Tests for pyrsss.gnss.bias (bias removal / calibration)."""
from datetime import timedelta

import numpy as np
import pytest

from pyrsss.gnss.bias import CalibratedArc, calibrate_arcs, calibrate_dcb, estimate_receiver_bias
from pyrsss.gnss.constants import GPS_EPOCH, NS_TO_TECU, TECU_TO_NS
from pyrsss.gnss.level import LeveledArc


def make_leveled_arc(sat='G01', n=10, el=45.0, L_I=None, P_I=None):
    t0 = GPS_EPOCH + timedelta(days=1800 * 7, seconds=259200)
    arc = LeveledArc({'gps_time': [t0 + timedelta(seconds=30 * i)
                                   for i in range(n)],
                      'az': np.full(n, 180.0),
                      'el': np.full(n, el),
                      'satx': np.full(n, 1000.0),
                      'saty': np.full(n, 2000.0),
                      'satz': np.full(n, 3000.0),
                      'P_I': P_I if P_I is not None else np.full(n, 10.0),
                      'L_I': L_I if L_I is not None else np.full(n, 5.0)})
    arc.xyz = [1.0, 2.0, 3.0]
    arc.llh = [40.0, -88.0, 200.0]
    arc.stn = 'TEST'
    arc.recv_type = 'TEST'
    arc.sat = sat
    arc.L = 5.0
    arc.L_scatter = 0.1
    return arc


def test_estimate_receiver_bias():
    # L_I = model_stec + sat_bias_tecu + stn_bias_true => estimator
    # recovers stn_bias_true exactly (noise-free)
    stn_true = 3.25
    b_ns = 0.25
    sat_bias_tecu = -b_ns / TECU_TO_NS
    model = [20.0 + 0.01 * i for i in range(10)]
    L_I = np.array(model) + sat_bias_tecu + stn_true
    arc = make_leveled_arc(L_I=L_I)
    sat_biases = {'GPS': {1: (b_ns, 0.0)}}
    bias, sigma = estimate_receiver_bias([arc], [model], sat_biases)
    assert bias == pytest.approx(stn_true, rel=1e-9)
    assert sigma == pytest.approx(0.0, abs=1e-9)


def test_estimate_receiver_bias_rejects_non_gps():
    arc = make_leveled_arc(sat='R04')
    with pytest.raises(NotImplementedError):
        estimate_receiver_bias([arc], [[0.0] * 10], {'GPS': {}})


def test_calibrate_arcs():
    stn_bias = 1.5
    b_ns = 0.25
    sat_bias = -b_ns / TECU_TO_NS
    arc = make_leveled_arc(L_I=np.full(10, 5.0), P_I=np.full(10, 10.0))
    sat_biases = {'GPS': {1: (b_ns, 0.0)}}
    calibrated = calibrate_arcs([arc], sat_biases, stn_bias,
                                stn_bias_sigma=0.1)
    assert len(calibrated) == 1
    cal = calibrated[0]
    assert isinstance(cal, CalibratedArc)
    np.testing.assert_allclose(cal.sobs, 5.0 - (sat_bias + stn_bias))
    # code left uncorrected by Attila's method
    np.testing.assert_allclose(cal.sprn, 10.0)
    assert cal.sat_bias == pytest.approx(sat_bias)
    assert cal.stn_bias == pytest.approx(stn_bias)
    assert cal.stn_bias_sigma == pytest.approx(0.1)
    assert cal.sat == 'G01'
    assert cal.L == pytest.approx(5.0)


def test_calibrate_dcb():
    b_ns, s_ns = 0.25, -0.5
    arc = make_leveled_arc(sat='G01', L_I=np.full(10, 5.0), P_I=np.full(10, 10.0))
    sat_biases = {'GPS': {1: (b_ns, 0.0)}, 'GLONASS': {}}
    stn_biases = {'GPS': {'TEST': (s_ns, 0.0)}, 'GLONASS': {}}
    calibrated = calibrate_dcb([arc], sat_biases, stn_biases)
    cal = calibrated[0]
    sat_bias = b_ns * NS_TO_TECU
    stn_bias = s_ns * NS_TO_TECU
    np.testing.assert_allclose(cal.sobs, 5.0 + sat_bias + stn_bias)
    np.testing.assert_allclose(cal.sprn, 10.0 + sat_bias + stn_bias)
    assert cal.sat_bias == pytest.approx(sat_bias)
    assert cal.stn_bias == pytest.approx(stn_bias)


def test_calibrate_dcb_glonass():
    b_ns, s_ns = 0.1, 0.2
    arc = make_leveled_arc(sat='R04')
    sat_biases = {'GPS': {}, 'GLONASS': {4: (b_ns, 0.0)}}
    stn_biases = {'GPS': {}, 'GLONASS': {'TEST': (s_ns, 0.0)}}
    calibrated = calibrate_dcb([arc], sat_biases, stn_biases)
    assert calibrated[0].sat_bias == pytest.approx(b_ns * NS_TO_TECU)


def test_calibrate_dcb_unknown_satellite():
    arc = make_leveled_arc(sat='E01')
    with pytest.raises(ValueError):
        calibrate_dcb([arc], {'GPS': {}, 'GLONASS': {}},
                      {'GPS': {}, 'GLONASS': {}})
