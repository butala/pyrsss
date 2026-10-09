"""
The HDZ/XYZ declination rotation that replaces geomagio's converters,
validated against closed forms and fuzzed.

geomagio's XYZAlgorithm and get_obs_from_geo/get_geo_from_obs were this
rotation and nothing more; the formulas are orthogonal, so the pins are
the identities an orthogonal rotation must satisfy plus a round-trip
fuzz over random declinations and vectors.
"""

import numpy as np
import pytest

from pyrsss.mag.hdfio import fix_sign, hdx_to_xyz, xyz_to_hdx


def test_identity_at_zero_declination():
    H, E, Z = 20000.0, -137.5, 43000.0
    X, Y, Zp = hdx_to_xyz(H, E, Z, 0.0)
    assert (X, Y, Zp) == (H, E, Z)


def test_ninety_degrees_maps_north_to_east():
    # at D = 90 deg the magnetic north axis points geographic east:
    X, Y, Zp = hdx_to_xyz(20000.0, 0.0, 5000.0, 90.0)
    assert X == pytest.approx(0.0, abs=1e-10)   # 2e4 * cos(pi/2) = 1.2e-12
    assert Y == pytest.approx(20000.0, rel=1e-12)
    assert Zp == 5000.0


def test_one_eighty_flips_the_horizontal():
    X, Y, _ = hdx_to_xyz(1000.0, 200.0, 3.0, 180.0)
    assert X == pytest.approx(-1000.0, rel=1e-12)
    assert Y == pytest.approx(-200.0, rel=1e-12)


def test_z_is_invariant_always():
    X, Y, Zp = hdx_to_xyz(1.0, 2.0, -3.0, 37.3)
    assert Zp == -3.0


def test_zero_vector_stays_zero():
    X, Y, Zp = hdx_to_xyz(0.0, 0.0, 0.0, 11.5)
    assert (X, Y, Zp) == (0.0, 0.0, 0.0)


def test_magnitude_of_the_horizontal_is_preserved():
    # the rotation is orthogonal: |(X, Y)| = |(H, E)|
    H, E = 1234.5, -678.9
    X, Y, _ = hdx_to_xyz(H, E, 0.0, 23.4)
    assert np.hypot(X, Y) == pytest.approx(np.hypot(H, E), rel=1e-14)


def test_a_known_case_by_hand():
    # D = 30 deg: X = H cos30 - E sin30, Y = H sin30 + E cos30
    X, Y, _ = hdx_to_xyz(10.0, 0.0, 0.0, 30.0)
    assert X == pytest.approx(10.0 * np.sqrt(3) / 2, rel=1e-14)
    assert Y == pytest.approx(5.0, rel=1e-14)


def test_round_trip_is_the_identity_fuzz():
    rng = np.random.default_rng(42)
    D = rng.uniform(-180.0, 180.0, 10000)
    H = rng.normal(0.0, 30000.0, 10000)
    E = rng.normal(0.0, 3000.0, 10000)
    Z = rng.normal(0.0, 40000.0, 10000)
    X, Y, Zp = hdx_to_xyz(H, E, Z, D)
    H2, E2, Z2 = xyz_to_hdx(X, Y, Zp, D)
    assert np.allclose(H2, H, rtol=0, atol=1e-10)
    assert np.allclose(E2, E, rtol=0, atol=1e-10)
    assert np.allclose(Z2, Z, rtol=0, atol=1e-12)
    assert np.allclose(np.hypot(X, Y), np.hypot(H, E), rtol=1e-14)


def test_the_inverse_is_the_rotation_by_minus_d():
    rng = np.random.default_rng(7)
    X, Y, Z = rng.normal(size=(3, 1000))
    D = rng.uniform(-180.0, 180.0, 1000)
    H, E, _ = xyz_to_hdx(X, Y, Z, D)
    X2, Y2, _ = hdx_to_xyz(H, E, Z, D)
    assert np.allclose(X2, X, rtol=0, atol=1e-12)
    assert np.allclose(Y2, Y, rtol=0, atol=1e-12)


def test_arrays_broadcast_over_declination():
    H, E, Z = np.ones(3), np.zeros(3), np.zeros(3)
    X, Y, _ = hdx_to_xyz(H, E, Z, [0.0, 90.0, 180.0])
    assert X[0] == pytest.approx(1.0, rel=1e-12)
    assert Y[1] == pytest.approx(1.0, rel=1e-12)
    assert X[2] == pytest.approx(-1.0, rel=1e-12)


def test_fix_sign_wraps_negative_declination():
    N = 360 * 60 * 10
    assert fix_sign(0.0) == 0.0
    assert fix_sign(-1.0) == pytest.approx(N - 1.0)
    assert fix_sign(N - 1.0) == N - 1.0
    with pytest.raises(AssertionError):
        fix_sign(-N - 1.0)


def test_dec_override_needs_no_igrf():
    """'decbas' in the header short-circuits the model -- pure and pinned."""
    from pyrsss.mag.hdfio import get_dec_tenths_arcminute

    assert get_dec_tenths_arcminute({'decbas': 1234.5}, '2020-01-01') == 1234.5
    assert get_dec_tenths_arcminute({'decbas': -1.0}, '2020-01-01') == \
        pytest.approx(360 * 60 * 10 - 1.0)


def test_the_converters_round_trip_through_the_rotation():
    """xy2df and he2df are the two directions of one rotation.

    The geomagio stream detour is gone; what must hold is that composing
    the two DataFrame converters is the identity -- and that a missing
    column is a message, not a KeyError.
    """
    import pandas as pd

    from pyrsss.mag.iaga2hdf import he2df, xy2df

    # IAGA records carry a DatetimeIndex; the declination lookup reads
    # the first epoch when the header has no 'decbas'
    df = pd.DataFrame({'B_H': [20000.0, 21000.0],
                       'B_E': [-137.0, 88.0],
                       'B_Z': [43000.0, 41000.0]},
                      index=pd.to_datetime(['2020-06-01T00:00:00',
                                            '2020-06-01T00:00:01']))
    header = {'decbas': 115 * 60 * 10,          # 11.5 deg in tenths of arcmin
              'Geodetic Latitude': 60.0,
              'Geodetic Longitude': 20.0}
    df2 = xy2df(df, header)
    assert {'B_X', 'B_Y'} <= set(df2.columns)
    df3 = he2df(df2, header)
    assert np.allclose(df3['B_H'], df['B_H'], rtol=0, atol=1e-10)
    assert np.allclose(df3['B_E'], df['B_E'], rtol=0, atol=1e-10)
    assert 'B_F' in df3.columns                     # synthesized once

    # and the explicit failure modes: a frame missing the other side's
    # components is a message naming what is needed, not a KeyError
    bare = pd.DataFrame({'B_Z': [1.0]},
                        index=pd.to_datetime(['2020-06-01T00:00:00']))
    with pytest.raises(ValueError, match='xy2df needs'):
        xy2df(bare, header)
    with pytest.raises(ValueError, match='he2df needs'):
        he2df(bare, header)
