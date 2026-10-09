"""Field-of-view geometry: rectangular pyramids and coronagraph annuli.

Everything here is directions, not pixels: an instrument's FOV is a set of
rays leaving its aperture, and two instruments *overlap* when their ray
sets intersect the same part of the sky. The models are pure numpy -- no
sunpy, no VTK -- so overlap and intercalibration math runs anywhere and
tests pin it against closed forms.

The frame convention is the tangent-plane one SphericalCT's D1 work settled
(see ``sphericalct.UniformPixelGrid.Config.isotropic``): for an observer at
``p`` looking at the Sun's centre, the boresight is ``b = -p/|p|``, the
local east ``e`` and north ``n`` complete the frame, and a pixel at
tangent offsets ``(tx, ty)`` lies along ``normalize(b + tx*e + ty*n)``.
A rectangular pyramid is exactly a box ``|tx| <= T_x``, ``|ty| <= T_y`` in
those coordinates -- which is what makes containment and overlap linear
algebra rather than projection code.

Distances are R_sun and angles are radians unless a name says otherwise.
"""

from dataclasses import dataclass

import numpy as np

ARCSEC = np.pi / (180.0 * 3600.0)          # 1 arcsec in radians
RSUN_KM = 6.957e8


def observer_frame(position):
    """
    The unit boresight ``b`` (towards the Sun's centre), local east ``e``
    and local north ``n`` for an observer at ``position``.

    At the ecliptic this is the fixed spherical basis, which is the identity
    that lets GOLD-era code and this module agree exactly there (the same
    argument as the isotropic ray frame). Poles are degenerate in longitude
    but the frame is still orthonormal.
    """
    p = np.asarray(position, dtype=float)
    r = np.linalg.norm(p)
    if r == 0.0:
        raise ValueError('observer at the Sun: no line of sight')
    b = -p / r
    # local north: the global z axis projected off the boresight
    z = np.array([0.0, 0.0, 1.0])
    n = z - np.dot(z, b) * b
    if np.linalg.norm(n) < 1e-12:          # observer over the pole
        n = np.array([0.0, 1.0, 0.0])
    n /= np.linalg.norm(n)
    e = np.cross(n, b)                     # right-handed: e x n = b?
    # (n, b, e) are orthonormal; fix the sign so +e is the direction of
    # increasing heliographic longitude as seen from the observer.
    if np.dot(e, np.cross(z, b)) < 0:
        e = -e
        n = np.cross(b, e)
    return b, e, n


@dataclass(frozen=True)
class RectPyramid:
    """A rectangular pyramid of rays: ``|tx| <= tan_x``, ``|ty| <= tan_y``.

    The half-extents are tangents of the half-angles -- for a plate scale
    of ``s`` arcsec and ``nx`` x ``ny`` pixels, ``tan_x = (nx/2) * s`` in
    radians to first order (exact: ``tan((nx/2) * s)``; at LASCO's 23.8"/px
    the two differ in the 8th digit, and the tangent form is what the
    pixel grid actually measures).
    """

    tan_x: float
    tan_y: float

    @classmethod
    def from_plate_scale(cls, plate_scale_arcsec, nx, ny=None):
        """The pyramid a ``nx`` x ``ny`` tangent-plane imager subtends."""
        ny = nx if ny is None else ny
        return cls(tan_x=np.tan(0.5 * nx * plate_scale_arcsec * ARCSEC),
                   tan_y=np.tan(0.5 * ny * plate_scale_arcsec * ARCSEC))

    def corner_rays(self, position):
        """
        The four corner directions of the far rectangle, and the apex.
        Returns ``(apex, corners)`` with corners in (+e+n, +e-n, -e-n, -e+n)
        order; the pyramid is the convex hull of the apex and the rectangle
        at unit boresight distance.
        """
        b, e, n = observer_frame(position)
        apex = np.asarray(position, dtype=float)
        corners = []
        for sx, sy in ((1, 1), (1, -1), (-1, -1), (-1, 1)):
            d = b + sx * self.tan_x * e + sy * self.tan_y * n
            corners.append(d / np.linalg.norm(d))
        return apex, np.array(corners)

    def contains(self, position, direction):
        """
        True when the ray from ``position`` along ``direction`` is inside
        the pyramid (on the boresight's side).
        """
        b, e, n = observer_frame(position)
        d = np.asarray(direction, dtype=float)
        d = d / np.linalg.norm(d)
        w_b = float(np.dot(d, b))
        if w_b <= 0.0:
            return False
        return (abs(float(np.dot(d, e)) / w_b) <= self.tan_x * (1 + 1e-12)
                and abs(float(np.dot(d, n)) / w_b) <= self.tan_y * (1 + 1e-12))

    def solid_angle(self):
        """The pyramid's solid angle, closed form (4 atan of the products)."""
        tx, ty = self.tan_x, self.tan_y
        return 4.0 * np.arctan(tx * ty / np.sqrt(1 + tx * tx + ty * ty))


@dataclass(frozen=True)
class AnnulusFOV:
    """A coronagraph's annulus on the sky: impact parameter in [r_in, r_out].

    The rays are every direction from the observer whose closest approach
    to the Sun's centre lands in the annulus -- which is what a coronagraph
    with a circular occulter actually selects. ``sector_deg`` limits the
    azimuth as seen in the observer's own (e, n) frame; None is the full
    circle.
    """

    r_inner_rsun: float
    r_outer_rsun: float
    sector_deg: float = None

    def impact(self, position, direction):
        """
        The ray's impact parameter (closest-approach distance to the Sun's
        centre), in the same length units as ``position``.
        """
        p = np.asarray(position, dtype=float)
        d = np.asarray(direction, dtype=float)
        d = d / np.linalg.norm(d)
        return float(np.linalg.norm(p - np.dot(p, d) * d))

    def contains(self, position, direction):
        """True when the ray's impact parameter is inside the annulus."""
        p = np.asarray(position, dtype=float)
        d = np.asarray(direction, dtype=float)
        d = d / np.linalg.norm(d)
        if float(np.dot(d, -p / np.linalg.norm(p))) <= 0.0:
            return False                      # looking away from the Sun
        b = self.impact(p, d)
        if not (self.r_inner_rsun <= b <= self.r_outer_rsun):
            return False
        if self.sector_deg is not None:
            boresight, e, n = observer_frame(p)
            half = np.deg2rad(0.5 * self.sector_deg)
            az = np.arctan2(float(np.dot(d, n)), float(np.dot(d, e)))
            if abs(az) > half:
                return False
        return True


def from_registry(fov):
    """
    Build the geometry model for a :class:`pyrsss.solar.registry.FOV`.
    """
    if fov.kind == 'rect':
        return RectPyramid.from_plate_scale(fov.plate_scale_arcsec,
                                            fov.nx, fov.ny)
    return AnnulusFOV(fov.r_inner_rsun, fov.r_outer_rsun, fov.sector_deg)


def sample_directions(model, position, n=32):
    """
    A grid of unit directions covering *model*'s FOV as seen from
    ``position`` (``n`` x ``n`` samples). For an annulus the samples cover
    the bounding box of the annulus in (e, n) tangent coordinates and are
    then filtered by ``contains``.
    """
    b, e, n_ = observer_frame(position)
    if isinstance(model, RectPyramid):
        tx = np.linspace(-model.tan_x, model.tan_x, n)
        ty = np.linspace(-model.tan_y, model.tan_y, n)
    else:
        # tangent half-width of the annulus's outer edge
        w = np.tan(np.arcsin(min(model.r_outer_rsun
                                 / np.linalg.norm(position), 1.0)))
        tx = np.linspace(-w, w, n)
        ty = np.linspace(-w, w, n)
    out = []
    for ux in tx:
        for uy in ty:
            d = b + ux * e + uy * n_
            d = d / np.linalg.norm(d)
            if model.contains(position, d):
                out.append(d)
    return np.array(out) if out else np.zeros((0, 3))
