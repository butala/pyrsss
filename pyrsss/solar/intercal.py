"""The intercalibration study: two calibrated images, one constant.

Given two :class:`pyrsss.solar.image.SkyImage` s that share sky, this is
the whole study: resample B onto A's pixels over the overlap (restricted
to A's own coronagraph annulus where there is one), take A/B, and report
the intercalibration constant and its drift with off-axis angle --
plus the ratio map itself, because "where do they disagree" is half the
question. :func:`plot` draws both.

The same-observer case (LASCO C2 ``pB`` against ``c2I``) is the first
thing the machinery pays for: the ratio is the polarization factor of
Thomson scattering, measured from real data, as a function of distance
from Sun centre. The cross-mission case is the same code with two
positions.
"""

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np

from .fov import observer_frame
from .image import SkyImage
from .overlap import Comparison, pixel_directions, resample

logger = logging.getLogger('pyrsss.solar.intercal')


@dataclass(frozen=True)
class IntercalReport:
    """The result of one pairwise study."""

    id_a: str
    id_b: str
    time: datetime
    comparison: Comparison
    ratio_map: np.ndarray        # A/B on A's grid, NaN outside the overlap
    annulus_mask: np.ndarray     # A's validity mask (bool), same grid

    @property
    def ratio_median(self):
        return self.comparison.ratio_median

    def __str__(self):
        return (f'{self.id_a} / {self.id_b} at {self.time:%Y-%m-%d %H:%M}: '
                f'{self.comparison}')


def _annulus_mask(image):
    """Which of A's pixels lie inside its own coronagraph annulus."""
    if image.annulus is None:
        return np.ones(image.shape, dtype=bool)
    ny, nx = image.shape
    b, e, n = observer_frame(image.position)
    dirs = pixel_directions(image.model, image.position, image.shape)
    mask = np.zeros(image.shape, dtype=bool)
    for i in range(ny):
        for j in range(nx):
            mask[i, j] = image.annulus.contains(image.position, dirs[i, j])
    return mask


def study(a, b, min_samples=25, n_profile=5):
    """
    The intercalibration of image *a* against image *b*.

    Pixels of A are used only where (1) A says they are valid (inside its
    annulus mask and finite and positive) and (2) B's resampled value is
    finite and positive. Raises ``ValueError`` under *min_samples*
    samples -- a comparison across a dozen pixels is not a calibration.
    """
    if not isinstance(a, SkyImage) or not isinstance(b, SkyImage):
        raise TypeError('study takes SkyImage pairs (see pyrsss.solar.image)')
    dirs = pixel_directions(a.model, a.position, a.shape)
    flat = dirs.reshape(-1, 3)
    sampled = resample(b.data, b.position, b.model, flat)
    a_flat = np.asarray(a.data, dtype=float).reshape(-1)
    valid = _annulus_mask(a).reshape(-1)
    good = valid & np.isfinite(sampled) & np.isfinite(a_flat) \
        & (sampled > 0) & (a_flat > 0)
    n = int(good.sum())
    if n < min_samples:
        raise ValueError(f'overlap too small for calibration: {n} samples')

    ratio = a_flat[good] / sampled[good]
    ratio_map = np.full(a_flat.shape, np.nan)
    ratio_map[good] = ratio
    b_vec, e_vec, n_vec = observer_frame(a.position)
    tan_off = np.array([np.hypot(np.dot(d, e_vec), np.dot(d, n_vec))
                        / max(abs(np.dot(d, b_vec)), 1e-15)
                        for d in flat[good]])
    frac = tan_off / max(a.model.tan_x, a.model.tan_y)
    edges = np.linspace(0.0, 1.0, n_profile + 1)
    centers, profile = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (frac >= lo) & (frac < hi)
        centers.append(0.5 * (lo + hi))
        profile.append(float(np.median(ratio[sel])) if sel.any() else np.nan)
    comparison = Comparison(
        ratio_median=float(np.median(ratio)),
        ratio_spread=float(1.4826 * np.median(np.abs(ratio
                                                     - np.median(ratio)))),
        n_samples=n,
        offaxis_bins=np.array(centers),
        ratio_profile=np.array(profile))
    logger.info('%s / %s: %s', a.instrument_id, b.instrument_id, comparison)
    return IntercalReport(id_a=a.instrument_id, id_b=b.instrument_id,
                          time=a.time, comparison=comparison,
                          ratio_map=ratio_map.reshape(a.shape),
                          annulus_mask=valid.reshape(a.shape))


def plot(report, path=None):
    """
    The ratio map and the ratio's off-axis profile, as one figure.
    Returns the figure (and saves it when *path* is given).
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    masked = np.ma.masked_invalid(report.ratio_map)
    im = axes[0].imshow(masked, origin='lower', cmap='coolwarm')
    axes[0].set_title(f'{report.id_a} / {report.id_b} '
                      f'(median {report.ratio_median:.4f})')
    fig.colorbar(im, ax=axes[0], label='A / B')
    c = report.comparison
    axes[1].plot(c.offaxis_bins, c.ratio_profile, 'o-', color='0.2')
    axes[1].axhline(c.ratio_median, color='crimson', ls='--',
                    label=f'median {c.ratio_median:.4f}')
    axes[1].set_xlabel('off-axis angle (fraction of half-width)')
    axes[1].set_ylabel('median A / B')
    axes[1].set_title(f'ratio profile ({c.n_samples} samples)')
    axes[1].legend()
    fig.tight_layout()
    if path is not None:
        fig.savefig(path, dpi=130)
    return fig


def main(argv=None):
    """``pyrsss-solar-intercal``: compare two calibrated FITS images."""
    import argparse

    import matplotlib
    matplotlib.use('Agg')

    from .image import load_fits

    parser = argparse.ArgumentParser(
        description='Intercalibrate two calibrated images over their '
                    'overlapping field of view.')
    parser.add_argument('a', type=Path, help="image A (the numerator's FITS)")
    parser.add_argument('b', type=Path, help='image B')
    parser.add_argument('--ida', default='lasco_c2', help='registry id of A')
    parser.add_argument('--idb', default='lasco_c2', help='registry id of B')
    parser.add_argument('--out', type=Path, default=Path('.'),
                        help='report directory (default: .)')
    args = parser.parse_args(argv)

    image_a = load_fits(args.a, args.ida)
    image_b = load_fits(args.b, args.idb)
    report = study(image_a, image_b)
    args.out.mkdir(parents=True, exist_ok=True)
    text_path = args.out / 'intercal_report.txt'
    text_path.write_text(str(report) + '\n')
    fig_path = args.out / 'intercal_report.png'
    plot(report, fig_path)
    print(report)
    print(f'wrote {text_path} and {fig_path}')
    return 0


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    raise SystemExit(main())
