"""The FOV overlap survey: which instruments share sky, and when.

Before any intercalibration, the question is who *can* be compared: two
instruments overlap when their fields of view cover the same plasma at
the same time. This module answers that from the registry alone -- no
network, no data -- because FOV geometry and ephemeris are knowable in
advance.

Two levels of overlap, deliberately distinguished:

* **annulus band overlap** -- the impact-parameter ranges intersect
  (2-6 R_sun for LASCO C2 against 2.5-15 for COR2 gives 2.5-6). This is
  the radial band both *can* see.
* **viewpoint separation** -- the angle between the two observers as
  seen from the Sun. This is the catch for pB: Thomson scattering's
  polarization depends on the scattering angle, so the same plasma gives
  different pB at two viewpoints, and above a few degrees the raw
  comparison means nothing without the scattering-geometry correction.
  Below ~5 deg (L1 against Earth orbit) the raw ratio is meaningful.

The survey's verdict is per pair: shared band, separation at the epoch,
and whether raw pB comparison is honest.
"""

from dataclasses import dataclass

import numpy as np

from . import registry
from .ephemeris import earth_circular, ground_site


@dataclass(frozen=True)
class PairOverlap:
    """The overlap verdict for one instrument pair at one epoch."""

    id_a: str
    id_b: str
    band: tuple            # shared impact range (R_sun), or None
    separation_deg: float  # observer-observer angle at the Sun
    raw_pb_ok: bool        # raw pB comparison honest (separation < 5 deg)

    def __str__(self):
        band = (f'{self.band[0]:.2f}-{self.band[1]:.2f} R_sun'
                if self.band else 'no shared band')
        if self.raw_pb_ok:
            verdict = 'raw pB ok'
        elif self.band is None:
            verdict = 'no shared band'
        else:
            verdict = 'check separation (ephemeris) / Thomson geometry'
        return (f'{self.id_a:8s} x {self.id_b:8s}  {band:>18s}  '
                f'sep {self.separation_deg:5.1f} deg  {verdict}')


def impact_band(inst):
    """An instrument's impact-parameter coverage (R_sun) as (lo, hi).

    Coronagraphs state it directly; an imager's tangent-plane extent
    projected at 1 AU is the practical band for a limb search.
    """
    fov = inst.fov
    if fov.kind == 'annulus':
        return (float(fov.r_inner_rsun), float(fov.r_outer_rsun))
    half = np.deg2rad(0.5 * fov.nx * fov.plate_scale_arcsec / 3600.0)
    return (0.0, float(np.tan(half) * 215.0))


def _position(inst, time):
    """The observer's position: HORIZONS when it answers, circular if not.

    The circular fallback puts *every* spacecraft at Earth's 1 AU position,
    so separations computed from it are meaningless between missions (a bug
    the first survey run showed: SOHO against STEREO-A came out 0.0 deg).
    ``separation_deg`` therefore marks the estimate it used.
    """
    if inst.observer.startswith('ground:'):
        return ground_site(time), 'ground'
    try:
        from .ephemeris import horizons
        return horizons(inst.observer, time), 'horizons'
    except Exception:
        return earth_circular(time), "circular"


def separation_deg(inst_a, inst_b, time):
    """The observer-observer angle at the Sun's centre (degrees), and how.

    Returns ``(deg, source)`` where source is 'horizons', 'ground', or
    'circular' -- the last one means both spacecraft were put at Earth's
    position and the number is only a floor.
    """
    if inst_a.observer == inst_b.observer:
        return 0.0, 'same-observer'   # one platform: no ephemeris needed
    a, src_a = _position(inst_a, time)
    b, src_b = _position(inst_b, time)
    cos = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    deg = float(np.rad2deg(np.arccos(np.clip(cos, -1.0, 1.0))))
    if src_a == 'circular' or src_b == 'circular':
        source = 'circular'
    elif src_a == 'ground' or src_b == 'ground':
        source = 'ground'
    else:
        source = 'horizons'
    return deg, source


def pair_overlap(inst_a, inst_b, time, sep_limit_deg=5.0):
    """The overlap verdict for one pair at one epoch."""
    lo_a, hi_a = impact_band(inst_a)
    lo_b, hi_b = impact_band(inst_b)
    lo, hi = max(lo_a, lo_b), min(hi_a, hi_b)
    band = (lo, hi) if hi > lo else None
    sep, source = separation_deg(inst_a, inst_b, time)
    return PairOverlap(id_a=inst_a.id, id_b=inst_b.id,
                       band=band, separation_deg=sep,
                       raw_pb_ok=band is not None and sep < sep_limit_deg
                       and source != 'circular')


def survey(time, polarized_only=True, ids=None):
    """
    Every pair of coronagraphs (optionally only the pB ones), surveyed.

    Returns the ``PairOverlap`` list, sorted with the usable pairs first
    (shared band and honest raw comparison). This is the table a
    comparison session starts from.
    """
    if ids:
        insts = [registry.get(i) for i in ids]
    else:
        insts = list(registry.coronagraphs(polarized_only=polarized_only))
    out = []
    for i in range(len(insts)):
        for j in range(i + 1, len(insts)):
            out.append(pair_overlap(insts[i], insts[j], time))
    return sorted(out,
                  key=lambda p: (not p.raw_pb_ok, p.band is None,
                                 -(p.band[1] - p.band[0]) if p.band else 0.0))


def format_survey(pairs):
    """The table as text, best pairs first."""
    lines = [str(p) for p in pairs]
    return '\n'.join(lines)
