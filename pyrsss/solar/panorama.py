"""The comparison plate: coronagraph images at one epoch, on one scale.

What an intercalibration session opens with: put the frames side by side,
Sun-centred, in **R_sun**, and see them. This module draws that plate --
one panel per image, identical spatial extent (the tangent-plane
coordinate scaled by the observer distance, which is exactly the impact
parameter to first order), the Sun's recorded pixel marked, the limb and
the annulus edges circled, and the epoch checked: "same epoch" is an
assertion, not an assumption.

Colour is deliberately *per panel* by default -- pB and total brightness
are different quantities with different unit conventions (IDOC's own
units warning), and a shared ramp across unlike products would invent a
comparison the data do not have. ``--shared-scale`` exists for like
products (pB against pBs), where one ramp is the honest picture.
"""

import logging
from pathlib import Path

import numpy as np

from .image import SkyImage

logger = logging.getLogger('pyrsss.solar.panorama')


def extent_rsun(image):
    """The half-width of an image's tangent plane, in R_sun at the observer.

    ``tan_x * |position|``: the impact parameter a ray at the frame's edge
    reaches to first order -- the axis a coronagraph actually measures in.
    """
    return image.model.tan_x * float(np.linalg.norm(image.position))


def sun_offset_rsun(image):
    """
    Where the Sun's centre sits in the R_sun plot frame: the shifted
    pyramid's boresight offset, scaled like the axes. (Zero for a
    perfectly centred product.)
    """
    r = float(np.linalg.norm(image.position))
    return (getattr(image.model, 'dx', 0.0) * r,
            getattr(image.model, 'dy', 0.0) * r)


def check_epoch(images, tol_s=60.0):
    """True when every image is within *tol_s* of the first one's epoch.

    Returns ``(ok, message)``: the plate prints it, because "same epoch"
    is the claim the comparison rests on and a two-hour-old frame beside a
    current one is a different picture of a different corona.
    """
    t0 = images[0].time
    worst = max(abs((img.time - t0).total_seconds()) for img in images)
    ok = worst <= tol_s
    return ok, (f'same epoch ({t0:%Y-%m-%d %H:%M:%S}, max offset '
                f'{worst:.0f} s)' if ok else
                f'EPOCHS DIFFER: up to {worst:.0f} s from {t0:%H:%M:%S}')


def panorama(images, path=None, shared_scale=False, title=None, dpi=130):
    """
    The plate itself: one panel per image, one spatial scale.

    *images* is a list of :class:`pyrsss.solar.image.SkyImage`. The
    spatial extent is the largest of the frames' (so nothing is cropped),
    identical in every panel. Returns the figure; saves it to *path*.
    """
    import matplotlib.pyplot as plt

    if not images:
        raise ValueError('the plate needs at least one image')
    ok, epoch_note = check_epoch(images)
    half = max(extent_rsun(img) for img in images)
    n = len(images)
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 5.8), squeeze=False)
    axes = axes[0]

    # one ramp only when asked: like products, one quantity
    levels = None
    if shared_scale:
        lo = min(float(np.nanmin(img.data)) for img in images)
        hi = max(float(np.nanmax(img.data)) for img in images)
        levels = (lo, hi)

    for ax, img in zip(axes, images):
        data = np.asarray(img.data, dtype=float)
        shown = np.ma.masked_invalid(data)
        kwargs = dict(origin='lower', cmap='gray',
                      extent=(-half, half, -half, half),
                      aspect='equal', interpolation='nearest')
        if levels is not None:
            kwargs.update(vmin=levels[0], vmax=levels[1])
        im = ax.imshow(shown, **kwargs)
        # the limb and the annulus this product's mask implies
        for radius, style, label in (
                (1.0, '-', 'limb'),
                (getattr(img.annulus, 'r_inner_rsun', None), ':', 'r_in'),
                (getattr(img.annulus, 'r_outer_rsun', None), ':', 'r_out')):
            if radius is not None and radius < 0.98 * half:
                ax.add_patch(plt.Circle((0, 0), radius, fill=False,
                                        color='white' if style == '-'
                                        else '0.7',
                                        lw=0.9, ls=style, label=label))
        sx, sy = sun_offset_rsun(img)
        ax.plot([sx], [sy], '+', color='crimson', ms=12, mew=1.4,
                label='XSUN/YSUN')
        ax.set_xlim(-half, half)
        ax.set_ylim(-half, half)
        ax.set_xlabel('R_sun (tangent plane)')
        ax.set_ylabel('R_sun')
        caption = f'{img.instrument_id} {img.label}'.strip()
        ax.set_title(f'{caption}\n{img.time:%Y-%m-%d %H:%M:%S}', fontsize=10)
        fig.colorbar(im, ax=ax, shrink=0.82,
                     label='product units' if levels is None else None)
        ax.legend(loc='upper right', fontsize=7, framealpha=0.3)

    note = epoch_note if ok else '!! ' + epoch_note
    fig.suptitle(title or f'coronagraph comparison --- {note}', fontsize=12,
                 color='black' if ok else 'crimson')
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    if path is not None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi)
        logger.info('wrote %s', path)
    return fig


def main(argv=None):
    """``pyrsss-solar-panorama``: side-by-side coronagraph frames."""
    import argparse

    import matplotlib
    matplotlib.use('Agg')

    from .image import load_fits

    parser = argparse.ArgumentParser(
        description='Show coronagraph images at one epoch on one scale.')
    parser.add_argument('fits', nargs='+', type=Path,
                        help='the calibrated FITS files, in panel order')
    parser.add_argument('--ids', default=None,
                        help='comma-separated registry ids (default: lasco_c2 '
                             'for every panel)')
    parser.add_argument('--shared-scale', action='store_true',
                        help='one colour ramp for all panels (like products)')
    parser.add_argument('--out', type=Path, default=Path('panorama.png'),
                        help='output figure (default: panorama.png)')
    args = parser.parse_args(argv)

    ids = ([s.strip() for s in args.ids.split(',')] if args.ids
           else ['lasco_c2'] * len(args.fits))
    if len(ids) != len(args.fits):
        raise SystemExit('--ids must name one registry id per FITS file')
    images = [load_fits(path, inst) for path, inst in zip(args.fits, ids)]
    panorama(images, args.out, shared_scale=args.shared_scale)
    print(f'wrote {args.out}')
    return 0


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    raise SystemExit(main())


def overlap_zoom(a, b, band=None, path=None, title=None, dpi=130):
    """
    Zoom into the two instruments' *intersecting region*: the shared
    impact-parameter band, cropped from both frames on one R_sun scale.

    Three panels: A's band, B's band (resampled onto A's grid so the two
    line up pixel for pixel), and their ratio where both are positive.
    ``band`` defaults to the shared annulus range of the two registry
    entries. The ratio is a *raw* pB comparison only when the observers
    are close (see :mod:`pyrsss.solar.survey`): Thomson scattering's
    polarization depends on the scattering angle, so distant viewpoints
    disagree for physical reasons, and the title says so when the pair
    is one of those.
    """
    import matplotlib.pyplot as plt

    from .overlap import pixel_directions, resample
    from .survey import impact_band, pair_overlap
    from . import registry

    if band is None:
        lo_a, hi_a = impact_band(registry.get(a.instrument_id))
        lo_b, hi_b = impact_band(registry.get(b.instrument_id))
        band = (max(lo_a, lo_b), min(hi_a, hi_b))
        if band[1] <= band[0]:
            raise ValueError('the two FOVs share no impact band')
    lo, hi = band
    r = float(np.linalg.norm(a.position))
    half = a.model.tan_x * r
    n = a.shape[0]
    yy, xx = np.mgrid[0:n, 0:n]
    rr = np.hypot((xx - n / 2 + 0.5) / (n / 2) * half,
                  (yy - n / 2 + 0.5) / (n / 2) * half)
    crop = (rr >= lo) & (rr <= hi)
    if not crop.any():
        raise ValueError("the two FOVs share no impact band")
    dirs = pixel_directions(a.model, a.position, a.shape).reshape(-1, 3)
    b_on_a = resample(b.data, b.position, b.model, dirs).reshape(a.shape)
    a_data = np.asarray(a.data, dtype=float)
    ratio = np.where(crop & (a_data > 0) & (b_on_a > 0),
                     a_data / np.where(b_on_a > 0, b_on_a, 1), np.nan)

    try:
        verdict = pair_overlap(registry.get(a.instrument_id),
                               registry.get(b.instrument_id), a.time)
    except KeyError:
        # product-level comparison (pB against c2I): no registry pair to
        # judge, and the two share one observer anyway
        verdict = "same observer (product-level comparison)"
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8))
    ext = (-half, half, -half, half)
    zoom = (max(lo - 0.5, -half), min(hi + 0.5, half))
    for ax, data, name in ((axes[0], a_data, a.label or a.instrument_id),
                           (axes[1], b_on_a, b.label or b.instrument_id),
                           (axes[2], ratio, f'{a.label}/{b.label} ratio')):
        shown = np.ma.masked_where(~crop, data)
        im = ax.imshow(shown, origin='lower', extent=ext, cmap='gray'
                       if name != f'{a.label}/{b.label} ratio' else 'coolwarm',
                       aspect='equal', interpolation='nearest')
        ax.set_xlim(zoom)
        ax.set_ylim(zoom)
        ax.set_title(name, fontsize=10)
        fig.colorbar(im, ax=ax, shrink=0.85)
    fig.suptitle(title or
                 f'overlap band {lo:.2f}-{hi:.2f} R_sun --- {verdict}',
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    if path is not None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi)
        logger.info('wrote %s', path)
    return fig
