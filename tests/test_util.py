"""Tests for pyrsss.util helpers (search, nan, date, angle, interval, lagrange)."""
from datetime import date, datetime

import numpy as np
import pytest

from pyrsss.util.angle import convert_lon, deg2tenths_of_arcminute
from pyrsss.util.date import (date2dt, fromJ2000, lt2ut, toJ2000, ut2lt)
from pyrsss.util.interval import closed, closed_open, length
from pyrsss.util.lagrange import multivariate_lagrange, poly_power_seq
from pyrsss.util.nan import nan_helper, nan_interp
from pyrsss.util.search import find_ge, find_gt, find_le, find_lt, index


A = [1, 3, 5]


def test_search_helpers():
    assert index(A, 3) == 1
    with pytest.raises(ValueError):
        index(A, 2)
    assert find_lt(A, 3) == (0, 1)
    assert find_le(A, 3) == (1, 3)
    assert find_le(A, 4) == (1, 3)
    assert find_gt(A, 3) == (2, 5)
    assert find_ge(A, 4) == (2, 5)
    assert find_ge(A, 5) == (2, 5)
    with pytest.raises(ValueError):
        find_le(A, 0)
    with pytest.raises(ValueError):
        find_gt(A, 5)


def test_nan_interp():
    y = np.array([0.0, np.nan, 2.0])
    z = nan_interp(y, silent=True)
    np.testing.assert_allclose(z, [0.0, 1.0, 2.0])
    nans, x = nan_helper(np.array([1.0, np.nan]))
    assert nans.tolist() == [False, True]
    assert x(nans).tolist() == [1]


def test_date_helpers():
    dt = datetime(2014, 6, 7, 12, 30)
    assert fromJ2000(toJ2000(dt)) == dt
    assert toJ2000(datetime(2000, 1, 1, 12)) == 0.0
    lon = 90.0
    assert ut2lt(dt, lon) == dt.replace(hour=18, minute=30)
    assert lt2ut(ut2lt(dt, lon), lon) == dt
    assert date2dt(date(2014, 6, 7)) == datetime(2014, 6, 7)


def test_angle_helpers():
    assert convert_lon(270.0) == -90.0
    assert convert_lon(10.0) == 10.0
    assert convert_lon(180.0) == 180.0
    assert deg2tenths_of_arcminute(1.0) == 600.0


def test_interval_helpers():
    i = closed_open(0.0, 2.0)
    assert 0.0 in i and 2.0 not in i
    assert length(i) == 2.0
    assert length(closed_open(1.0, float('inf'))) == float('inf')
    assert 1.0 in closed(1.0, 2.0) and 2.0 in closed(1.0, 2.0)
    # None bounds are unbounded
    assert 1e9 in closed_open(None, None)
    assert 1.0 not in closed_open(2.0, None)
    assert 5.0 in closed_open(2.0, None)


def test_lagrange_dimension_check():
    assert sorted(poly_power_seq(2, 1)) == [(0, 0), (0, 1), (1, 0)]
    # p must equal n + m choose n (here m=2, n=1 -> 3, but only 1 point)
    with pytest.raises(ValueError):
        multivariate_lagrange([(0.0, 0.0, 1.0)], 1)
