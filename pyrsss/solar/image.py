"""SkyImage: a calibrated image plus the geometry that places it.

The intercalibration step needs one thing from a FITS file -- the array,
the plate scale, where the Sun's centre is, and where the observer stood
-- and every archive spells those differently. This module normalizes:
``load_fits`` turns a file into a :class:`SkyImage` (data + a
:class:`pyrsss.solar.fov.RectPyramid` + observer position + the
coronagraph's annulus mask), and the per-archive *recipes* say which
header keys were trusted.

The first recipe is **IDOC's LASCO C2 pB family** (``kfcorona_sph``:
``pB``, ``pBs``, ``c2I``, ...), verified 2026-10 against real 2024
products: 512 x 512, ``XSUN``/``YSUN`` hold the Sun's centre in pixels,
and no ``CDELT`` at all -- the plate scale is the instrument's, so it
comes from the registry (23.8 "/px, the constant the zheyuan fetcher's
``create_pB_map`` also used). Units follow the product's own (IDOC's
are ~1e-10 B_sun; a *ratio* of two images in the same family cancels
whatever the convention is).
"""

import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np

from .ephemeris import earth_circular
from .fov import RectPyramid, from_registry
from . import registry

logger = logging.getLogger('pyrsss.solar.image')


@dataclass(frozen=True)
class SkyImage:
    """A calibrated image and everything needed to place it on the sky."""

    data: np.ndarray             # (ny, nx) the measurement
    model: RectPyramid           # the tangent-plane geometry
    position: np.ndarray         # observer, R_sun, heliocentric ecliptic
    time: datetime
    instrument_id: str
    annulus: object = None       # registry AnnulusFOV, or None for an imager
    label: str = ''              # the product, e.g. 'pB' -- panel captions

    @property
    def shape(self):
        return self.data.shape


def load_idoc(path, instrument_id='lasco_c2', position=None):
    """
    An IDOC ``kfcorona_sph`` product (LASCO C2: ``pB``, ``pBs``, ``c2I``).

    The recipe, from the files themselves: ``XSUN``/``YSUN`` is the Sun's
    centre in pixels and the plate scale is the instrument's (the registry
    value; the headers carry none). ``DATE_OBS`` + ``TIME_OBS`` give the
    epoch. Verified against 2024-11-03 products.
    """
    from astropy.io import fits

    with fits.open(path) as hdul:
        hdr = hdul[0].header
        data = np.asarray(hdul[0].data, dtype=float)
    ny, nx = data.shape
    sun_x = float(hdr['XSUN'])
    sun_y = float(hdr['YSUN'])
    inst = registry.get(instrument_id)
    model = RectPyramid.from_plate_scale(inst.fov.plate_scale_arcsec, nx, ny)
    time = datetime.strptime(
        f"{hdr['DATE_OBS']} {hdr['TIME_OBS'].split('.')[0]}", '%Y/%m/%d %H:%M:%S')
    if position is None:
        position = earth_circular(time)   # SOHO is at L1; see ephemeris
    annulus = from_registry(inst.fov) if inst.fov.kind == 'annulus' else None
    # XSUN/YSUN may sit off the array centre: the pyramid is the one whose
    # boresight passes through the Sun's recorded pixel (a half-pixel
    # convention question SphericalCT settled as 0-based centres, and a
    # shifted Sun must not silently skew every ray).
    model = _shifted_pyramid(inst.fov.plate_scale_arcsec, nx, ny, sun_x, sun_y)
    label = Path(path).stem.split('_')[-1]   # 'c2_pB.fts' -> 'pB'
    return SkyImage(data=data, model=model, position=np.asarray(position),
                    time=time, instrument_id=instrument_id, annulus=annulus,
                    label=label)


def _shifted_pyramid(scale, nx, ny, sun_x, sun_y):
    """
    The pyramid whose boresight passes through the Sun's recorded pixel
    ``(sun_x, sun_y)`` rather than the array centre: the tangent offsets
    are shifted by the centre's offset, which is what keeps an off-centre
    Sun from silently skewing every ray.
    """
    base = RectPyramid.from_plate_scale(scale, nx, ny)
    dx = (nx / 2.0 - sun_x) / (nx / 2.0) * base.tan_x
    dy = (ny / 2.0 - sun_y) / (ny / 2.0) * base.tan_y
    return ShiftedPyramid(base, dx, dy)


@dataclass(frozen=True)
class ShiftedPyramid:
    """A RectPyramid whose boresight is offset in tangent coordinates."""

    base: RectPyramid
    dx: float
    dy: float

    @property
    def tan_x(self):
        return self.base.tan_x

    @property
    def tan_y(self):
        return self.base.tan_y

    def corner_rays(self, position):
        return self.base.corner_rays(position)

    def contains(self, position, direction):
        return self.base.contains(position, direction)

    def solid_angle(self):
        return self.base.solid_angle()

    def to_tangent(self, position, direction):
        """A direction's tangent coordinates in this shifted frame."""
        from .overlap import project

        t = project(direction, position, self.base)
        return None if t is None else (t[0] - self.dx, t[1] - self.dy)


def load_fits(path, instrument_id='lasco_c2', recipe='idoc', position=None):
    """
    Load a calibrated product by recipe name ('idoc' today; the other
    archives land here as their header recipes are pinned).
    """
    if recipe == 'idoc':
        return load_idoc(path, instrument_id, position=position)
    raise ValueError(f'unknown header recipe {recipe!r}')
