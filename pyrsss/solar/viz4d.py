"""The 4-D FOV scene: source locations and rectangular pyramids over time.

What a tomographer wants to *see* before writing an inversion: where every
observer stood, and what each one's field of view covered, as the geometry
moves. This module draws that with pyviz4d -- its ``TemporalActor`` loop is
the fourth dimension, so the scene is one object animated by playback (the
slider and ``space`` are ``Viewer4D.add_playback_ui``'s), not a pile of
stills.

The scene is the Sun at the origin, one marker per observer, its orbit
track over the window, and the observer's FOV as a **rectangular pyramid
wireframe** (four apex edges plus the far rectangle) rebuilt every tick
from :mod:`pyrsss.solar.ephemeris` and :mod:`pyrsss.solar.fov`. A
coronagraph's annulus rides over its pyramid as the two circles at its
impact radii -- the same convention SphericalCT's plate draws.

Import of pyviz4d is deferred: the registry, geometry and overlap math of
this package work without it, and the tests skip the scene when it is
absent. ``uv sync --extra solar-viz`` or a sibling ``~/src/pyviz4d``
checkout (which is how SphericalCT consumes it) both work.
"""

import logging
from datetime import datetime

import numpy as np

logger = logging.getLogger('pyrsss.solar.viz4d')

SUN_COLOR = (0.95, 0.85, 0.35)
SUN_RADIUS_RSUN = 1.0


def _pyviz4d():
    try:
        import pyviz4d
        from pyviz4d import series, viz
        return pyviz4d, series, viz
    except ImportError as e:
        raise ImportError(
            'the FOV scene needs pyviz4d: uv sync --extra solar-viz '
            'or set PYTHONPATH to a pyviz4d checkout') from e


def orbit_track(position_fn, t0, t1, n=32):
    """Sampled observer positions over ``[t0, t1]`` as an (n, 3) array."""
    return np.array([position_fn(t0 + (t1 - t0) * k / (n - 1))
                     for k in range(n)])


class FOVPyramidActor:
    """
    One instrument's FOV pyramid as a moving wireframe.

    Duck-types pyviz4d's ``TemporalActor`` (``.actor`` + ``update(t)``),
    so ``Viewer4D.add_actor`` registers it and the playback loop moves it.
    ``position_fn`` maps the loop's time coordinate (float in ``[0, 1]``)
    to the observer's position; the pyramid is rebuilt from it every tick.
    """

    def __init__(self, model, position_fn, color=(0.4, 0.8, 1.0),
                 width=1.5):
        import vtk

        self.model = model
        self.position_fn = position_fn
        self.points = vtk.vtkPoints()
        self.points.SetNumberOfPoints(5)          # apex + 4 corners
        lines = vtk.vtkCellArray()
        for a, b in ((0, 1), (0, 2), (0, 3), (0, 4),   # apex edges
                     (1, 2), (2, 3), (3, 4), (4, 1)):  # far rectangle
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, a)
            line.GetPointIds().SetId(1, b)
            lines.InsertNextCell(line)
        poly = vtk.vtkPolyData()
        poly.SetPoints(self.points)
        poly.SetLines(lines)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(poly)
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(mapper)
        self.actor.GetProperty().SetColor(*color)
        self.actor.GetProperty().SetLineWidth(width)
        self.actor.GetProperty().SetLighting(0)
        self.update(0.0)

    def update(self, t):
        """Rebuild the five vertices from the observer at loop time *t*."""
        position = self.position_fn(t)
        apex, corners = self.model.corner_rays(position)
        self.points.SetPoint(0, *apex)
        for i, corner in enumerate(corners, 1):
            self.points.SetPoint(i, *corner)
        self.points.Modified()


class ObserverTrackActor:
    """The orbit track and the moving observer marker (static track, moving point)."""

    def __init__(self, position_fn, t0, t1, color=(1.0, 1.0, 1.0), n=32):
        import vtk

        self.position_fn = position_fn
        self.t0, self.t1 = t0, t1
        track = orbit_track(position_fn, t0, t1, n=n)
        self.points = vtk.vtkPoints()
        for p in track:
            self.points.InsertNextPoint(*p)
        lines = vtk.vtkCellArray()
        for i in range(len(track) - 1):
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, i)
            line.GetPointIds().SetId(1, i + 1)
            lines.InsertNextCell(line)
        poly = vtk.vtkPolyData()
        poly.SetPoints(self.points)
        poly.SetLines(lines)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(poly)
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(mapper)
        self.actor.GetProperty().SetColor(*color)
        self.actor.GetProperty().SetLineWidth(1.0)
        self.actor.GetProperty().SetOpacity(0.5)
        self.actor.GetProperty().SetLighting(0)

    def update(self, t):
        """The track is static; a temporal actor must still respond."""
        return None


def sun_actor(radius=SUN_RADIUS_RSUN, color=SUN_COLOR, opacity=0.95):
    """The Sun: an opaque sphere at the origin, the scene's one landmark."""
    import vtk

    source = vtk.vtkSphereSource()
    source.SetRadius(radius)
    source.SetThetaResolution(32)
    source.SetPhiResolution(32)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(source.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)
    actor.GetProperty().SetOpacity(opacity)
    actor.GetProperty().SetLighting(0)
    return actor


def build_scene(instruments, positions, t0, t1, size=(1400, 900),
                offscreen=True):
    """
    The 4-D scene for a set of registry instruments.

    *instruments* is a list of ``(Instrument, color)``; *positions* maps an
    instrument id to a ``position_fn(t)`` (loop time 0..1), so the caller
    chooses the ephemeris backend once per observer. Returns the
    ``Viewer4D`` and the temporal actors, ready for ``add_playback_ui`` +
    ``start()`` or for stills.
    """
    pyviz4d, series, viz = _pyviz4d()
    from . import registry
    from .fov import from_registry

    viewer = viz.Viewer4D(size=size, bg_color=(0.08, 0.08, 0.12))
    viewer.ren_win.SetOffScreenRendering(1 if offscreen else 0)
    viewer.add_actor(sun_actor())
    temporals = []
    for inst, color in instruments:
        model = from_registry(inst.fov)
        fn = positions[inst.id]
        pyramid = FOVPyramidActor(model, fn, color=color)
        viewer.add_actor(pyramid)
        viewer.add_actor(ObserverTrackActor(fn, t0, t1, color=color))
        temporals.append((inst, pyramid))
        logger.info('%s: %s FOV wired', inst.id, inst.fov.kind)
    return viewer, temporals


def write_stills(viewer, temporals, times, out_dir, prefix='fov'):
    """
    Render one PNG per entry of *times* (loop coordinates in [0, 1]),
    using pyviz4d's own render path. Returns the paths.
    """
    from pathlib import Path

    pyviz4d, series, viz = _pyviz4d()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for k, t in enumerate(times):
        for _, actor in temporals:
            actor.update(t)
        viewer.ren_win.Render()
        path = out_dir / f'{prefix}_{k:03d}.png'
        pyviz4d.primitives._render_window_to_png(viewer.ren_win, str(path))
        paths.append(path)
        logger.info('wrote %s', path)
    return paths


def main(argv=None):
    """``pyrsss-solar-fov-scene``: the 4-D FOV scene for registry instruments."""
    import argparse
    from pathlib import Path

    from . import registry
    from .ephemeris import earth_circular, ground_site, horizons

    parser = argparse.ArgumentParser(
        description='Draw instrument locations and FOV pyramids over time.')
    parser.add_argument('--ids', default='lasco_c2,kcor',
                        help='comma-separated registry ids '
                             '(default: lasco_c2,kcor)')
    parser.add_argument('--t0', default=None, help='start (ISO); default: -1d')
    parser.add_argument('--t1', default=None, help='end (ISO); default: now')
    parser.add_argument('--steps', type=int, default=5,
                        help='stills to write when not --live (default: 5)')
    parser.add_argument('--out', type=Path, default=Path('.'),
                        help='output directory (default: .)')
    parser.add_argument('--live', action='store_true',
                        help='open the interactive 4-D window instead')
    args = parser.parse_args(argv)

    t1 = datetime.fromisoformat(args.t1) if args.t1 else datetime.now()
    t0 = datetime.fromisoformat(args.t0) if args.t0 else t1.replace(
        hour=t1.hour - 24)

    def make_fn(inst):
        span = (t1 - t0).total_seconds()
        name = inst.observer

        def fn(u):
            t = datetime.fromtimestamp(t0.timestamp() + u * span)
            return ground_site(t) if name.startswith('ground:') \
                else horizons(name, t)
        return fn

    instruments, positions = [], {}
    for name in args.ids.split(','):
        inst = registry.get(name.strip())
        instruments.append((inst, (0.5, 0.8, 1.0)))
        positions[inst.id] = make_fn(inst)

    viewer, temporals = build_scene(instruments, positions, t0, t1)
    if args.live:
        viewer.add_playback_ui(1.0)
        viewer.start()
        return 0
    paths = write_stills(viewer, temporals,
                         np.linspace(0.0, 1.0, args.steps), args.out)
    for path in paths:
        print(path)
    return 0


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    raise SystemExit(main())
