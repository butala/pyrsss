"""
Shuttle demo: the time transport and the time strip, offline.

A small self-contained scene -- two ground-based instruments over the
default window -- so the transport (keyboard, wheel jog, drag-to-scrub,
fling coast, zoom) can be exercised without network or mission data::

    uv run --extra solar-viz pyrsss-solar-shuttle-demo

``--steps N`` renders stills instead of opening the window (for
headless checks); everything else is the interactive demo.
"""

import logging
from datetime import datetime, timedelta

import numpy as np

from ..util.cli import add_log_options, configure_logging
from . import registry
from .ephemeris import orbit_interpolator
from .shuttle import Shuttle

logger = logging.getLogger('pyrsss.solar.shuttle_demo')


def main(argv=None):
    import argparse
    from pathlib import Path

    from .viz4d import build_scene, write_stills

    parser = argparse.ArgumentParser(
        description='Demo of the shuttle time transport (offline scene).')
    parser.add_argument('--ids', default='kcor,ucomp',
                        help='registry ids (default: kcor,ucomp)')
    parser.add_argument('--t0', default=None, help='start (ISO)')
    parser.add_argument('--t1', default=None, help='end (ISO)')
    parser.add_argument('--days', type=float, default=7.0,
                        help='window length when times are given as one '
                             'epoch (default: 7)')
    parser.add_argument('--steps', type=int, default=0,
                        help='render this many stills instead of the '
                             'interactive window')
    parser.add_argument('--out', type=Path, default=Path('.'),
                        help='stills output directory')
    add_log_options(parser)
    args = parser.parse_args(argv)
    configure_logging(args)

    t1 = datetime.fromisoformat(args.t1) if args.t1 else datetime.now().replace(
        minute=0, second=0, microsecond=0)
    t0 = datetime.fromisoformat(args.t0) if args.t0 else t1 - timedelta(
        days=args.days)

    from .viz4d import PALETTE

    instruments, positions = [], {}
    for k, name in enumerate(args.ids.split(',')):
        inst = registry.get(name.strip())
        instruments.append((inst, PALETTE[k % len(PALETTE)]))
        positions[inst.id] = orbit_interpolator(inst.observer, t0, t1)

    viewer, temporals, timeline = build_scene(instruments, positions, t0, t1,
                                              offscreen=args.steps > 0)
    if args.steps > 0:
        paths = write_stills(viewer, temporals,
                             np.linspace(0.0, 1.0, args.steps), args.out)
        for path in paths:
            print(path)
        return 0

    Shuttle(viewer, t0, t1, timeline)
    logger.info('shuttle demo: drag the bottom strip, fling it, '
                'scroll the scene, +/- to zoom')
    viewer.start()
    return 0


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    raise SystemExit(main())
