"""PUNCH data: the full-Sun polarized-brightness mosaics.

PUNCH (four Wide Field Imagers plus one Narrow Field Imager) is the one
mission that ships *polarized* white light as its headline product: L2
polarized mosaics in Mean Solar Brightness, and L3 ``PAM``/``CAM``
(polarized / clear) low-noise mosaics on HEALPix-style projections. For
Thomson-scattering tomography the L3_PAM is the measurement.

The supported reader is **punchbowl**, the mission's own package (it knows
the trefoil mosaic layout and the MZP polarization triplet structure); this
module is the thin cache/fetch layer around it, not a re-implementation.
Import of ``punchbowl`` is deferred so the rest of ``pyrsss.solar`` works
without it.

Verified 2026-10 against the PUNCH data pages (data.nasa.gov datasets
``PUNCH_NFI_WFI_Mosaic_Level2_PTM`` etc.) and punchbowl's data docs.
"""

import logging
from pathlib import Path

logger = logging.getLogger('pyrsss.solar.punch')


def _punchbowl():
    """
    Import punchbowl or raise with the install hint.
    """
    try:
        import punchbowl  # noqa: F401
    except ImportError as e:
        raise ImportError(
            'PUNCH support needs the mission package: '
            'pip install punchbowl') from e
    return punchbowl


def levels():
    """
    The product levels this module knows, with the tomography-relevant
    one marked.
    """
    return {
        'L2_PTM': 'polarized science mosaic (MZP triplet + uncertainty)',
        'L3_PAM': 'polarized low-noise mosaic -- the tomography input',
        'L3_CAM': 'clear (unpolarized) low-noise mosaic',
    }


def fetch(day, level='L3_PAM', out_dir=Path('.'), force=False):
    """
    Fetch one PUNCH product day into *out_dir* via punchbowl.

    *day* is a ``datetime.date``/``datetime``; *level* is a key of
    ``levels()``. Returns the local paths punchbowl wrote. Network use is
    punchbowl's, and so is its layout -- this function only standardizes
    where things land and what gets logged.
    """
    punchbowl = _punchbowl()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logger.info('fetching PUNCH %s for %s into %s', level, day, out_dir)
    if hasattr(punchbowl, 'fetch'):
        return punchbowl.fetch(day, level=level, out_dir=str(out_dir),
                               force=force)
    raise RuntimeError(
        'this punchbowl version exposes no fetch(); see '
        'https://punchbowl.readthedocs.io for its current data interface')


def main(argv=None):
    """``pyrsss-solar-fetch-punch``: one day of PUNCH L3_PAM."""
    import argparse
    import datetime

    parser = argparse.ArgumentParser(
        description='Fetch PUNCH polarized mosaics (via punchbowl).')
    parser.add_argument('date', help='YYYY-MM-DD')
    parser.add_argument('--level', default='L3_PAM',
                        choices=sorted(levels()), help='(default: L3_PAM)')
    parser.add_argument('--out', type=Path, default=Path('.'),
                        help='cache directory (default: .)')
    args = parser.parse_args(argv)

    day = datetime.datetime.strptime(args.date, '%Y-%m-%d')
    paths = fetch(day, level=args.level, out_dir=args.out)
    for path in paths:
        print(path)
    return 0


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    raise SystemExit(main())
