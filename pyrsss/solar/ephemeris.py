"""Observer positions: where an instrument stood when it looked.

Every FOV question reduces to this: the geometry of
:mod:`pyrsss.solar.fov` places rays leaving an *apex*, and the apex is the
observer's heliocentric position at the observation time. Three backends,
matching the registry's ``ephemeris`` field:

* ``'soho-orbit-file'`` -- the NRL predictive orbit records
  (:mod:`pyrsss.solar.soho`), which is the path the IDOC pB products
  assume. The HEC columns of the record are the position directly.
* ``'earth'`` -- a ground site. The Earth is put on a *circular* 1 AU
  ecliptic orbit here -- honest about it: that is good to about 3.6 R_sun
  in absolute position (the orbit's eccentricity), which decides "does
  this coronagraph see that one's FOV" but is **not** absolute
  registration. Call :func:`horizons` for that.
* ``'horizons'`` -- JPL HORIZONS via astroquery, for any body or
  spacecraft. This is the accurate one; it is import-guarded and needs
  network.

Positions are in R_sun, heliocentric ecliptic (the HAE-like frame the SOHO
orbit files' ``hec`` columns use): x towards the vernal point, z to the
north ecliptic pole. HORIZONS results are cached on disk (see
:func:`horizons`).
"""

import json
import logging
import os
from datetime import datetime
from pathlib import Path

import numpy as np

from .soho import closest_orbit, parse_orbit_file

logger = logging.getLogger('pyrsss.solar.ephemeris')

AU_RSUN = 149_597_870.7 / 6.957e5          # km/km: 215.03 R_sun
RSUN_KM = 6.957e5                          # for orbit-file conversions
EARTH_OBLIQUITY = np.deg2rad(23.4392911)   # for ground-site offsets

# Registry observer names -> HORIZONS body names.
HORIZONS_NAMES = {
    'spacecraft:SOHO': 'SOHO',
    # HORIZONS takes the plain names; 'STEREO-AHEAD'/'STEREO-BEHIND' are
    # not valid targets there (verified against the service 2026-10)
    'spacecraft:STEREO-A': 'STEREO-A',
    'spacecraft:STEREO-B': 'STEREO-B',
    'spacecraft:SDO': 'SDO',
    'spacecraft:PSP': 'PARKER SOLAR PROBE',
    'spacecraft:PUNCH': 'PUNCH',
    'spacecraft:SolarOrbiter': 'SOLAR ORBITER',
    'spacecraft:PROBA-3': 'PROBA-3',
    'spacecraft:PROBA-2': 'PROBA-2',
    'spacecraft:Aditya-L1': 'ADITYA-L1',
    'spacecraft:ASO-S': 'ASO-S',
    'spacecraft:Hinode': 'HINODE',
    'spacecraft:GOES': 'GOES-16',
}


def earth_circular(time):
    """
    The Earth's heliocentric ecliptic position on a circular 1 AU orbit.

    Mean longitude from the J2000 mean equinox: good to ~1.7% in distance
    and ~0.01 deg in angle over decades -- plenty for a FOV survey, not
    for absolute registration (see the module docstring).
    """
    t = _to_datetime(time)
    days = (t - datetime(2000, 1, 1, 12)).total_seconds() / 86400.0
    lon = np.deg2rad(100.466 + 0.9856474 * days)   # mean solar longitude
    return AU_RSUN * np.array([np.cos(lon), np.sin(lon), 0.0])


def ground_site(time, lat_deg=19.536, lon_deg=-155.576, height_m=3397.0):
    """
    A ground observatory's heliocentric position: Earth's vector plus the
    site offset (Mauna Loa's is the default, home of MLSO).

    The site offset is at most 1 R_sun from the Earth's centre -- it moves
    the apex by that much, which matters at the 1 R_sun inner edge and is
    what the last term is for. ``height_m`` enters the same way.
    """
    earth = earth_circular(time)
    t = _to_datetime(time)
    # sidereal rotation: one turn per sidereal day, prime meridian at J2000
    days = (t - datetime(2000, 1, 1, 12)).total_seconds() / 86400.0
    lon = np.deg2rad(lon_deg + 360.9856235 * days)
    lat = np.deg2rad(lat_deg)
    r = 1.0 + height_m / (6.371e6) * (6.371e6 / 6.957e8)  # ~1 R_earth in R_sun
    offset = r * np.array([np.cos(lat) * np.cos(lon),
                           np.cos(lat) * np.sin(lon) * np.cos(EARTH_OBLIQUITY)
                           - np.sin(lat) * np.sin(EARTH_OBLIQUITY),
                           np.cos(lat) * np.sin(lon) * np.sin(EARTH_OBLIQUITY)
                           + np.sin(lat) * np.cos(EARTH_OBLIQUITY)])
    return earth + offset


def observer_position(observer, time):
    """
    The position of a registry observer name at *time*: ground sites go
    through :func:`ground_site`, everything else through :func:`horizons`.
    """
    if observer.startswith('ground:'):
        return ground_site(time)
    return horizons(observer, time)


def orbit_interpolator(observer, t0, t1, n=64):
    """
    Return ``fn(u)`` mapping loop coordinates ``u`` in ``[0, 1]`` to the
    observer's interpolated position over the window ``[t0, t1]``.

    A HORIZONS query is a network round trip; sampling the window *n*
    times and interpolating bounds the cost and makes playback smooth
    (one day of 1 AU orbit is ~1 deg of arc, so linear interpolation
    between samples is far below the geometry's accuracy). Each sample
    is logged: this is where a scene spends its startup time.
    """
    span = t1 - t0
    logger.info('sampling ephemeris for %s: %d points over %s',
                observer, n, span)
    times = [t0 + span * k / (n - 1) for k in range(n)]
    table = np.empty((n, 3))
    for k, t in enumerate(times):
        logger.debug('ephemeris %s: %d/%d at %s',
                     observer, k + 1, n, t.isoformat())
        table[k] = observer_position(observer, t)
    logger.info('sampled ephemeris for %s (%d points)', observer, n)

    def fn(u):
        x = min(max(float(u), 0.0), 1.0) * (n - 1)
        i = min(int(x), n - 2)
        frac = x - i
        return (1.0 - frac) * table[i] + frac * table[i + 1]

    return fn


def soho_orbit_file(time, dat_path):
    """
    The SOHO position from a NRL orbit ``.DAT`` (:func:`soho.parse_orbit_file`
    + the closest record). The record's HEC columns are the answer.
    """
    records = parse_orbit_file(dat_path)
    record, offset = closest_orbit(records, _to_datetime(time))
    if abs(offset.total_seconds()) > 900:
        logger.warning('closest orbit record is %s from the request', offset)
    return np.array([record.hec_x_km, record.hec_y_km, record.hec_z_km]) \
        / RSUN_KM                    # km -> R_sun


def _cache_path():
    """The ephemeris cache file (``$PYRSSS_EPHEMERIS_CACHE`` overrides)."""
    import os

    default = Path.home() / '.cache' / 'pyrsss' / 'ephemeris.json'
    return Path(os.environ.get('PYRSSS_EPHEMERIS_CACHE', default))


def _cache_load(path):
    try:
        with open(path) as fid:
            return json.load(fid)
    except FileNotFoundError:
        return {}
    except (json.JSONDecodeError, OSError) as e:
        logger.warning('ignoring unreadable ephemeris cache %s (%s)', path, e)
        return {}


def _cache_store(path, cache):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp')
    with open(tmp, 'w') as fid:
        json.dump(cache, fid, indent=0, sort_keys=True)
    os.replace(tmp, path)


def _horizons_query(name, t):
    """
    One HORIZONS round trip: the body's heliocentric ecliptic position
    in R_sun at *t* (the network call the cache exists to avoid).
    """
    try:
        from astroquery.jplhorizons import Horizons
    except ImportError as e:
        raise ImportError(
            "HORIZONS needs astroquery: `uv sync --extra solar` (or "
            "`--extra solar-viz` for the FOV scene)") from e
    # epochs wants Julian *dates* in a list -- ISO strings give
    # "BATVAR: no TLIST values found" from the API itself (verified 2026-10).
    from astropy.time import Time

    obj = Horizons(id=name, epochs=[float(Time(t).jd)],
                   location='500@10')                 # Sun-centred
    vec = obj.vectors()
    x = float(vec['x'][0]) * AU_RSUN
    y = float(vec['y'][0]) * AU_RSUN
    z = float(vec['z'][0]) * AU_RSUN
    return np.array([x, y, z])


def horizons(body, time):
    """
    The accurate position: JPL HORIZONS through astroquery, in the
    heliocentric ecliptic frame, R_sun. ``body`` is a registry observer
    name or any HORIZONS name. Needs ``astroquery`` and the network.

    Positions are cached on disk -- a HORIZONS query is deterministic,
    so the cache never expires. It lives at
    ``~/.cache/pyrsss/ephemeris.json`` (``$PYRSSS_EPHEMERIS_CACHE``
    overrides the path).
    """
    name = HORIZONS_NAMES.get(body, body)
    t = _to_datetime(time)
    key = t.isoformat()
    path = _cache_path()
    cache = _cache_load(path)
    hit = cache.get(name, {}).get(key)
    if hit is not None:
        logger.debug('ephemeris cache hit %s at %s', name, key)
        return np.array(hit)
    logger.info('querying HORIZONS for %s at %s', name, key)
    xyz = _horizons_query(name, t)
    cache.setdefault(name, {})[key] = [float(v) for v in xyz]
    _cache_store(path, cache)
    return xyz


def observer_position(observer, time, **kw):
    """
    Dispatch on a registry entry's ``observer`` / ``ephemeris`` spelling.

    ``observer`` is e.g. ``'spacecraft:SOHO'``, ``'ground:MLSO'``; *kw*
    goes to the backend (``dat_path=`` for the orbit files,
    ``lat_deg``/``lon_deg`` for ground sites).
    """
    if observer.startswith('ground:'):
        return ground_site(time,
                           lat_deg=kw.get('lat_deg', 19.536),
                           lon_deg=kw.get('lon_deg', -155.576))
    if 'dat_path' in kw:
        return soho_orbit_file(time, kw['dat_path'])
    return horizons(observer, time)


def _to_datetime(time):
    if isinstance(time, datetime):
        return time
    return datetime.fromisoformat(str(time))
