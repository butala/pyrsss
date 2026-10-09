"""FOV overlap and intercalibration: who sees the same sky, and do they agree?

Two instruments can only be intercalibrated where they *both* look. This
module answers that in three steps, all in direction space:

1. :func:`sky_overlap` -- the common set of rays (and each FOV's fraction
   inside the other), from the geometry of :mod:`pyrsss.solar.fov`.
2. :func:`project` -- one direction's coordinates in another instrument's
   tangent plane, which is what moves a pixel from one image to the other.
3. :func:`compare` -- resample one image onto the other's pixels inside
   the overlap and report the ratio: median (the intercalibration
   constant), robust spread, and the ratio's drift with off-axis angle.

The images are arrays plus the geometry that places them -- no WCS objects,
no FITS: a caller that has read a FITS file already knows its plate scale
and observer. The math is numpy so the whole path is unit-testable against
closed forms.
"""

from dataclasses import dataclass

import numpy as np

from .fov import RectPyramid, observer_frame


@dataclass(frozen=True)
class Overlap:
    """The common sky of two FOVs: fractions, and the shared directions."""

    directions: np.ndarray      # (N, 3) unit vectors, the common rays
    fraction_a_in_b: float      # share of A's sampled rays inside B
    fraction_b_in_a: float

    @property
    def is_empty(self):
        return self.directions.shape[0] == 0


def sky_overlap(model_a, pos_a, model_b, pos_b, n=24):
    """
    Sample both FOVs and keep the directions each has in common.

    *n* x *n* is the sample grid per FOV; the fractions are over each
    FOV's own samples, so they answer "what share of A's field does B
    cover" -- the question an intercalibration study opens with.
    """
    from .fov import sample_directions

    da = sample_directions(model_a, pos_a, n=n)
    db = sample_directions(model_b, pos_b, n=n)
    common = np.array([d for d in da if model_b.contains(pos_b, d)])
    frac_a = len(common) / max(len(da), 1)
    frac_b = len(common) / max(len(db), 1)
    return Overlap(common, float(frac_a), float(frac_b))


def project(direction, position, model):
    """
    One direction's tangent coordinates ``(tx, ty)`` in an instrument's
    frame at ``position``, or None when the direction is outside the
    boresight's half-space (a coronagraph's far-side rays land there).
    """
    b, e, n = observer_frame(position)
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)
    w_b = float(np.dot(d, b))
    if w_b <= 0.0:
        return None
    return (float(np.dot(d, e)) / w_b, float(np.dot(d, n)) / w_b)


def pixel_directions(model, position, shape):
    """
    The unit ray through the centre of every pixel of a ``shape`` image,
    as an array shaped ``shape + (3,)``. ``model`` must be a
    :class:`RectPyramid` (an imager; a coronagraph's annulus mask only
    zeros pixels afterwards).
    """
    if not (hasattr(model, 'tan_x') and hasattr(model, 'tan_y')):
        raise TypeError('pixel grids need a tangent-plane model '
                        '(RectPyramid or ShiftedPyramid)')
    ny, nx = shape
    b, e, n = observer_frame(position)
    dx = getattr(model, 'dx', 0.0)      # a ShiftedPyramid's boresight offset
    dy = getattr(model, 'dy', 0.0)
    tx = (np.arange(nx) + 0.5 - nx / 2.0) / (nx / 2.0) * model.tan_x + dx
    ty = (np.arange(ny) + 0.5 - ny / 2.0) / (ny / 2.0) * model.tan_y + dy
    out = np.empty((ny, nx, 3))
    for i, uy in enumerate(ty):
        for j, ux in enumerate(tx):
            d = b + ux * e + uy * n
            out[i, j] = d / np.linalg.norm(d)
    return out


def resample(image_b, pos_b, model_b, directions, fill=np.nan):
    """
    Sample image B at arbitrary ``directions`` (bilinear in B's tangent
    plane). Out-of-view samples come back as *fill*.
    """
    image_b = np.asarray(image_b, dtype=float)
    ny, nx = image_b.shape
    out = np.full(len(directions), fill, dtype=float)
    for k, d in enumerate(directions):
        t = project(d, pos_b, model_b)
        if t is None:
            continue
        dx = getattr(model_b, 'dx', 0.0)
        dy = getattr(model_b, 'dy', 0.0)
        j = (t[0] - dx) / model_b.tan_x * (nx / 2.0) + nx / 2.0 - 0.5
        i = (t[1] - dy) / model_b.tan_y * (ny / 2.0) + ny / 2.0 - 0.5
        if not (0 <= j <= nx - 1 and 0 <= i <= ny - 1):
            continue
        i0, j0 = int(np.floor(i)), int(np.floor(j))
        i1, j1 = min(i0 + 1, ny - 1), min(j0 + 1, nx - 1)
        di, dj = i - i0, j - j0
        out[k] = (image_b[i0, j0] * (1 - di) * (1 - dj)
                  + image_b[i1, j0] * di * (1 - dj)
                  + image_b[i0, j1] * (1 - di) * dj
                  + image_b[i1, j1] * di * dj)
    return out


@dataclass(frozen=True)
class Comparison:
    """The intercalibration summary of two images over their overlap."""

    ratio_median: float          # A / B, the intercalibration constant
    ratio_spread: float          # robust half-width (1.4826 * MAD)
    n_samples: int
    offaxis_bins: np.ndarray     # bin centres (fraction of A's half-width)
    ratio_profile: np.ndarray    # median A/B in each bin

    def __str__(self):
        return (f'A/B = {self.ratio_median:.4f} +/- {self.ratio_spread:.4f} '
                f'over {self.n_samples} samples')


def compare(image_a, pos_a, model_a, image_b, pos_b, model_b,
            min_samples=25, n_profile=5):
    """
    The intercalibration constant of image A against image B: resample B
    onto A's pixels, take A/B where both are finite and positive, and
    report the median and the ratio's drift with off-axis angle.

    Raises ``ValueError`` when the overlap is too small to say anything --
    a comparison over five pixels is not a calibration.
    """
    image_a = np.asarray(image_a, dtype=float)
    dirs = pixel_directions(model_a, pos_a, image_a.shape)
    b_vec, e_vec, n_vec = observer_frame(pos_a)
    flat_dirs = dirs.reshape(-1, 3)
    sampled = resample(image_b, pos_b, model_b, flat_dirs)
    a_flat = image_a.reshape(-1)
    good = np.isfinite(sampled) & np.isfinite(a_flat) \
        & (sampled > 0) & (a_flat > 0)
    if int(good.sum()) < min_samples:
        raise ValueError(
            f'overlap too small for calibration: {int(good.sum())} samples')
    ratio = a_flat[good] / sampled[good]
    # off-axis angle as a fraction of A's half-width, for the profile
    tan_off = np.array([np.hypot(np.dot(d, e_vec), np.dot(d, n_vec))
                        / max(abs(np.dot(d, b_vec)), 1e-15)
                        for d in flat_dirs[good]])
    frac = tan_off / max(model_a.tan_x, model_a.tan_y)
    edges = np.linspace(0.0, 1.0, n_profile + 1)
    centers, profile = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (frac >= lo) & (frac < hi)
        centers.append(0.5 * (lo + hi))
        profile.append(float(np.median(ratio[sel])) if sel.any() else np.nan)
    return Comparison(
        ratio_median=float(np.median(ratio)),
        ratio_spread=float(1.4826 * np.median(np.abs(ratio
                                                     - np.median(ratio)))),
        n_samples=int(good.sum()),
        offaxis_bins=np.array(centers),
        ratio_profile=np.array(profile))
