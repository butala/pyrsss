"""
Broadcast ephemeris orbit propagation (GPS Keplerian and GLONASS ICD
models) and topocentric azimuth/elevation geometry. Positions are WGS84
ECEF [m].
"""
import math
from datetime import datetime, timedelta

import numpy as np

from ..util.position import Position

GPS_EPOCH = datetime(1980, 1, 6)

MU_GPS = 3.986005e14
"""GPS ICD gravitational constant [m^3/s^2]."""

OMEGA_E = 7.2921151467e-5
"""WGS84 earth rotation rate [rad/s]."""

MU_GLO = 3.9860044e14
"""GLONASS ICD gravitational constant [m^3/s^2]."""

J2_GLO = 1.0826257e-3
"""GLONASS ICD second zonal harmonic."""

AE_GLO = 6378136.0
"""GLONASS ICD earth equatorial radius [m]."""


def _solve_kepler(M, e, tol=1e-13, max_iter=20):
    """
    Solve Kepler's equation M = E - e sin E for the eccentric anomaly E.
    """
    E = M if e < 0.8 else math.pi
    for _ in range(max_iter):
        dE = (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        E -= dE
        if abs(dE) < tol:
            break
    return E


def gps_position(eph, t):
    """
    Return the ECEF [m] position of the GPS satellite with broadcast
    ephemeris *eph* (mapping of RINEX navigation parameter names, e.g.
    'sqrtA', 'M0', ...) at :class:`datetime` *t* (GPS ICD-200 Table
    20-IV).
    """
    e = eph['Eccentricity']
    A = eph['sqrtA']**2
    n = math.sqrt(MU_GPS / A**3) + eph['DeltaN']
    toe = GPS_EPOCH + timedelta(weeks=eph['GPSWeek'],
                                seconds=eph['Toe'])
    tk = (t - toe).total_seconds()
    E = _solve_kepler(eph['M0'] + n * tk, e)
    nu = math.atan2(math.sqrt(1 - e**2) * math.sin(E), math.cos(E) - e)
    phi = nu + eph['omega']
    u = phi + (eph['Cuc'] * math.cos(2 * phi) +
               eph['Cus'] * math.sin(2 * phi))
    i = (eph['Io'] + eph['IDOT'] * tk +
         eph['Cic'] * math.cos(2 * phi) + eph['Cis'] * math.sin(2 * phi))
    r = (A * (1 - e * math.cos(E)) +
         eph['Crc'] * math.cos(2 * phi) + eph['Crs'] * math.sin(2 * phi))
    Omega = (eph['Omega0'] + (eph['OmegaDot'] - OMEGA_E) * tk -
             OMEGA_E * eph['Toe'])
    x_p, y_p = r * math.cos(u), r * math.sin(u)
    x = x_p * math.cos(Omega) - y_p * math.sin(Omega) * math.cos(i)
    y = x_p * math.sin(Omega) + y_p * math.cos(Omega) * math.cos(i)
    z = y_p * math.sin(i)
    return np.array([x, y, z])


def _glonass_deriv(state, acc):
    """
    GLONASS ICD equations of motion in the Earth-fixed PZ-90 frame
    (state = [x, y, z, vx, vy, vz] all in SI, *acc* = lunisolar
    acceleration in [m/s^2]).
    """
    x, y, z, vx, vy, vz = state
    r2 = x * x + y * y + z * z
    r = math.sqrt(r2)
    r3 = r2 * r
    r5 = r3 * r2
    mu_r3 = MU_GLO / r3
    j2 = 1.5 * J2_GLO * AE_GLO**2 * MU_GLO / r5
    z2_r2 = z * z / r2
    w = OMEGA_E
    return np.array([
        vx,
        vy,
        vz,
        -mu_r3 * x - j2 * x * (1 - 5 * z2_r2) + w * w * x + 2 * w * vy + acc[0],
        -mu_r3 * y - j2 * y * (1 - 5 * z2_r2) + w * w * y - 2 * w * vx + acc[1],
        -mu_r3 * z - j2 * z * (3 - 5 * z2_r2) + acc[2],
    ])


def _rk4(state, dt, acc, nstep):
    """Propagate *state* by *dt* [s] in *nstep* Runge-Kutta 4 steps."""
    h = dt / nstep
    for _ in range(nstep):
        k1 = _glonass_deriv(state, acc)
        k2 = _glonass_deriv(state + h / 2 * k1, acc)
        k3 = _glonass_deriv(state + h / 2 * k2, acc)
        k4 = _glonass_deriv(state + h * k3, acc)
        state = state + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return state


def glonass_position(geph, t):
    """
    Return the ECEF [m] position of the GLONASS satellite with broadcast
    ephemeris *geph* (RINEX mapping: 'X', 'Y', 'Z' [km], 'VX', 'VY',
    'VZ' [km/s], 'AX', 'AY', 'AZ' [km/s^2], 'TauN'...) at :class:`datetime`
    *t*. Integrates the PZ-90 ICD equations of motion from the message
    epoch.
    """
    state = np.array([geph['X'] * 1e3, geph['Y'] * 1e3, geph['Z'] * 1e3,
                      geph['VX'] * 1e3, geph['VY'] * 1e3, geph['VZ'] * 1e3])
    acc = np.array([geph['AX'] * 1e3, geph['AY'] * 1e3, geph['AZ'] * 1e3])
    dt = (t - geph['_epoch']).total_seconds()
    if abs(dt) < 1e-9:
        return state[:3].copy()
    nstep = max(2, int(math.ceil(abs(dt) / 30.0)))
    state = _rk4(state, dt, acc, nstep)
    return state[:3].copy()


def satellite_position(nav, sat, t):
    """
    Return the ECEF [m] position of *sat* at :class:`datetime` *t* using
    the ephemeris parameters *nav* (sat -> {name -> value}).
    """
    if sat[0] == 'G':
        return gps_position(nav[sat], t)
    elif sat[0] == 'R':
        return glonass_position(nav[sat], t)
    raise ValueError('no orbit model for {}'.format(sat))


def azel(stn_llh, sat_xyz):
    """
    Return the (azimuth, elevation) [deg] of the ECEF point *sat_xyz*
    [m] observed from the geodetic station position *stn_llh*
    (lat [deg], lon [deg], height [m]). Azimuth is measured clockwise
    from North.
    """
    lat, lon, _ = stn_llh
    stn = Position(lat, lon, stn_llh[2],
                   Position.CoordinateSystem['geodetic'])
    dx, dy, dz = np.asarray(sat_xyz) - np.asarray(stn.xyz)
    phi = math.radians(lat)
    lam = math.radians(lon)
    sphi, cphi = math.sin(phi), math.cos(phi)
    slon, clon = math.sin(lam), math.cos(lam)
    e = -slon * dx + clon * dy
    n = -sphi * clon * dx - sphi * slon * dy + cphi * dz
    u = cphi * clon * dx + cphi * slon * dy + sphi * dz
    az = math.degrees(math.atan2(e, n)) % 360
    el = math.degrees(math.atan2(u, math.hypot(e, n)))
    return az, el
