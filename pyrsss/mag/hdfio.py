"""Pure HDF helpers and the HDZ/XYZ component rotation for the mag family.

Two things live here, and both are dependency-free:

* The HDF record I/O -- ``read_hdf``/``write_hdf`` (moved verbatim from
  ``iaga2hdf``, which imported geomagio for something else entirely and
  so dragged an unresolvable dependency under three innocent helpers).
* The declination rotation between the two local frames magnetometer
  data arrives in: **HDZ** (H = local magnetic north, E = local magnetic
  east, Z = down) and **XYZ** (X = geographic north, Y = geographic
  east, Z = down). This is a rotation of the horizontal *components* by
  the magnetic declination D, not a geodetic-to-geomagnetic coordinate
  transformation -- geomagio's ``XYZAlgorithm`` and
  ``get_obs_from_geo``/``get_geo_from_obs`` were this math and nothing
  more:

      X =  H cos D - E sin D        H =  X cos D + Y sin D
      Y =  H sin D + E cos D        E = -X sin D + Y cos D
      Z =  Z                        Z =  Z

  with D positive east of geographic north. The replacement is exact:
  the formulas are orthogonal rotations, |(X, Y)| = |(H, E)| and the
  round trip is the identity, which tests/test_mag_rotation.py pins
  against closed forms and over a fuzz of 10^4 random vectors.
"""

import os

import numpy as np
import pandas as pd

from ..util.angle import deg2tenths_of_arcminute


# --------------------------------------------------------------- HDF records

def write_hdf(hdf_fname, df, key, header):
    """
    Output the contents of *df* and *header* to the HDF file
    *hdf_fname* under identifier *key*.
    """
    with pd.HDFStore(hdf_fname) as store:
        store.put(key, df)
        store.get_storer(key).attrs.header = header
    return hdf_fname


def read_hdf(hdf_fname, key):
    """
    Read contents of HDF file *hdf_fname* associated with *key* and
    return a :class:`DataFrame`, header tuple.
    """
    if not os.path.isfile(hdf_fname):
        raise ValueError('file {} does not exist'.format(hdf_fname))
    with pd.HDFStore(hdf_fname) as store:
        df = store.get(key)
        try:
            header = store.get_storer(key).attrs.header
        except AttributeError:
            header = None
        return df, header


# --------------------------------------------------------- declination source

def fix_sign(x, N=360 * 60 * 10):
    """
    Convert negative tenths of arcminutes *x* to positive by checking
    bounds and taking the modulus N (360 degrees * 60 minutes per
    degree * 10 tenths per 1).
    """
    if x < 0:
        assert x > -N
        x += N
    assert x < N
    return x % N


def get_dec_tenths_arcminute(header, date):
    """
    The local magnetic declination of the sensor in *header* at *date*,
    in tenths of arcminutes (360 * 60 * 10 per turn).

    The header's own ``decbas`` override wins when present (it is the
    baseline the recording was reduced with); otherwise the declination
    is computed with ``ppigrf`` (pure-Python IGRF). The earlier
    implementation reached for a ``Point`` object that was never
    imported -- this function raised ``NameError`` before it computed
    anything -- so the source is now explicit and the failure is a
    message, not a missing symbol.
    """
    if 'decbas' in header:
        return fix_sign(float(header['decbas']))
    try:
        import ppigrf
    except ImportError as e:
        raise ImportError(
            'declination needs ppigrf (uv sync --extra mag) or a '
            "'decbas' entry in the IAGA header") from e
    date = pd.Timestamp(date).to_pydatetime()
    bx, by, _ = (float(v) for v in ppigrf.igrf(
        [date], float(header['Geodetic Latitude']),
        float(header['Geodetic Longitude']), float(header.get('Elevation', 0))))
    dec_deg = float(np.degrees(np.arctan2(by, bx)))
    import logging
    logging.getLogger('pyrsss.mag.hdfio').info(
        'declination %.4f deg from IGRF%s', dec_deg,
        f" ({header['IAGA CODE']})" if 'IAGA CODE' in header else '')
    return fix_sign(deg2tenths_of_arcminute(dec_deg))


# ------------------------------------------------------ the declination rotation

def hdx_to_xyz(H, E, Z, dec_deg):
    """
    Rotate field components from the magnetic (HDZ) to the geographic
    (XYZ) frame. *dec_deg* is the declination D, positive east of
    geographic north. Arrays broadcast; Z passes through.
    """
    D = np.radians(np.asarray(dec_deg, dtype=float))
    H = np.asarray(H, dtype=float)
    E = np.asarray(E, dtype=float)
    Z = np.asarray(Z, dtype=float)
    X = H * np.cos(D) - E * np.sin(D)
    Y = H * np.sin(D) + E * np.cos(D)
    return X, Y, Z


def xyz_to_hdx(X, Y, Z, dec_deg):
    """The inverse of :func:`hdx_to_xyz` (the rotation by ``-D``)."""
    D = np.radians(np.asarray(dec_deg, dtype=float))
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    Z = np.asarray(Z, dtype=float)
    H = X * np.cos(D) + Y * np.sin(D)
    E = -X * np.sin(D) + Y * np.cos(D)
    return H, E, Z
