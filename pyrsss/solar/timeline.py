"""
The time strip: a slim timeline along the bottom of the scene.

One rectangle, three jobs: *show* where in the window the scene is
(tick axis, instrument bands as the legend), *move* the clock (the
playhead; drag and fling policy live in :mod:`pyrsss.solar.shuttle`),
and *scale* the view (``zoom``). Drawing is plain VTK 2-D geometry in
normalized viewport coordinates -- one actor per colour, a fixed pool
of tick labels -- and the class duck-types ``TemporalActor`` (``.actor``
+ ``update(u)``) so the playhead rides the scene's clock.
"""

from datetime import datetime, timedelta

# strip rectangle in normalized viewport coordinates (x0, x1, y0, y1)
STRIP = (0.07, 0.93, 0.05, 0.135)

FRAME_COLOR = (0.45, 0.45, 0.52)
BAND_WIDTH = 3.0
PLAYHEAD_COLOR = (0.92, 0.28, 0.22)
LABEL_POOL = 10

_STEP_LADDER = [timedelta(minutes=1), timedelta(minutes=5),
                timedelta(minutes=15), timedelta(minutes=30),
                timedelta(hours=1), timedelta(hours=3), timedelta(hours=6),
                timedelta(hours=12), timedelta(days=1), timedelta(days=2),
                timedelta(days=7), timedelta(days=14), timedelta(days=30)]


def nice_tick_step(span):
    """A tick spacing between span/12 and span/5 from the step ladder."""
    for step in _STEP_LADDER:
        if span / step <= 12:
            return step
    return _STEP_LADDER[-1]


def ticks(t0, t1, view=(0.0, 1.0)):
    """Tick datetimes inside the visible view (``u`` window) of ``[t0, t1]``."""
    lo, hi = view
    start = t0 + (t1 - t0) * lo
    end = t0 + (t1 - t0) * hi
    step = nice_tick_step(end - start)
    first = start.replace(minute=0, second=0, microsecond=0)
    while first < start:
        first += step
    out, t = [], first
    while t <= end:
        out.append(t)
        t += step
    return out


def x_from_u(u, view=(0.0, 1.0)):
    """Strip x (normalized viewport) of the loop coordinate *u*."""
    lo, hi = view
    u = min(max(float(u), lo), hi)
    x0, x1 = STRIP[0], STRIP[1]
    return x0 + (x1 - x0) * (u - lo) / (hi - lo)


def u_from_x(x, view=(0.0, 1.0)):
    """Loop coordinate at strip x (clamped to the view)."""
    lo, hi = view
    x0, x1 = STRIP[0], STRIP[1]
    return min(max(lo + (hi - lo) * (x - x0) / (x1 - x0), lo), hi)


def hit(x, y):
    """True when the viewport point is inside the strip (with margin)."""
    x0, x1, y0, y1 = STRIP
    return (x0 - 0.02 <= x <= x1 + 0.02) and (y0 - 0.02 <= y <= y1 + 0.06)


class _Layer:
    """One colour's worth of line geometry."""

    def __init__(self, color, width=1.0):
        import vtk

        self.points = vtk.vtkPoints()
        self.lines = vtk.vtkCellArray()
        self.poly = vtk.vtkPolyData()
        self.poly.SetPoints(self.points)
        self.poly.SetLines(self.lines)
        coord = vtk.vtkCoordinate()
        coord.SetCoordinateSystemToNormalizedViewport()
        mapper = vtk.vtkPolyDataMapper2D()
        mapper.SetTransformCoordinate(coord)
        mapper.SetInputData(self.poly)
        self.actor = vtk.vtkActor2D()
        self.actor.SetMapper(mapper)
        self.actor.GetProperty().SetColor(*color)
        self.actor.GetProperty().SetLineWidth(width)

    def line(self, x1, y1, x2, y2):
        import vtk

        a = self.points.InsertNextPoint(x1, y1, 0.0)
        b = self.points.InsertNextPoint(x2, y2, 0.0)
        edge = vtk.vtkLine()
        edge.GetPointIds().SetId(0, a)
        edge.GetPointIds().SetId(1, b)
        self.lines.InsertNextCell(edge)

    def reset(self):
        self.points.Reset()
        self.lines.Reset()

    def done(self):
        self.points.Modified()
        self.lines.Modified()


class Timeline:
    """
    The drawn strip. *bands* is a list of ``(label, color)`` pairs, one
    per instrument; each gets a thin coloured bar (the legend *is* the
    data: every instrument covers the whole window here). The view is
    the visible ``u`` window, default everything.
    """

    def __init__(self, t0, t1, bands=(), n_frames=241):
        from pyviz4d.ui import create_text_actor

        self.t0, self.t1 = t0, t1
        self.bands = list(bands)
        self.n_frames = n_frames
        self.view = (0.0, 1.0)

        self._frame = _Layer(FRAME_COLOR, width=1.0)
        self._bands = [_Layer(color, width=BAND_WIDTH) for _, color in self.bands]
        self._playhead = _Layer(PLAYHEAD_COLOR, width=2.0)

        self.labels = []
        for _ in range(LABEL_POOL):
            label = create_text_actor('', 0.0, 0.0, font_size=11)
            label.GetTextProperty().SetColor(0.78, 0.78, 0.83)
            self.labels.append(label)

        self.update(0.0)

    @property
    def actor(self):
        """Primary VTK actor (the playhead; the rest are extras)."""
        return self._playhead.actor

    @property
    def extra_actors(self):
        """The remaining actors: frame, bands, tick labels."""
        return ([self._frame.actor]
                + [band.actor for band in self._bands]
                + self.labels)

    def _draw(self, u):
        x0, x1, y0, y1 = STRIP
        frame, bands, playhead = self._frame, self._bands, self._playhead
        for layer in [frame, playhead] + bands:
            layer.reset()

        # frame
        frame.line(x0, y0, x1, y0)
        frame.line(x0, y1, x1, y1)
        frame.line(x0, y0, x0, y1)
        frame.line(x1, y0, x1, y1)

        # instrument bands, staggered so none hides another
        for i, layer in enumerate(bands):
            y = y0 + (y1 - y0) * (0.30 + 0.18 * i)
            layer.line(x0 + 0.002, y, x1 - 0.002, y)

        # ticks
        for t in ticks(self.t0, self.t1, self.view):
            frac = ((t - self.t0).total_seconds()
                    / (self.t1 - self.t0).total_seconds())
            frame.line(x_from_u(frac, self.view), y0,
                       x_from_u(frac, self.view), y0 + 0.012)

        # playhead
        x = x_from_u(u, self.view)
        playhead.line(x, y0 - 0.012, x, y1 + 0.012)

        for layer in [frame, playhead] + bands:
            layer.done()

        # labels (fixed pool: retext and reposition, never re-add)
        ts = ticks(self.t0, self.t1, self.view)
        for k, label in enumerate(self.labels):
            if k < len(ts):
                t = ts[k]
                frac = ((t - self.t0).total_seconds()
                        / (self.t1 - self.t0).total_seconds())
                label.SetInput(t.strftime('%m-%d %H:%M'))
                label.SetPosition(x_from_u(frac, self.view) - 0.032,
                                  y0 - 0.032)
                label.SetVisibility(True)
            else:
                label.SetVisibility(False)

    def update(self, u):
        """Redraw the strip with the playhead at the loop coordinate *u*."""
        self._draw(u)

    def zoom(self, factor, center_u=None):
        """
        Zoom the visible view by *factor* (> 1 widens in time, < 1
        narrows), keeping *center_u* (default: the playhead) put.
        """
        if center_u is None:
            center_u = 0.5 * sum(self.view)
        lo, hi = self.view
        width = (hi - lo) * factor
        width = min(max(width, (self.n_frames - 1) ** -1, 1e-6), 1.0)
        frac = (center_u - lo) / (hi - lo) if hi > lo else 0.5
        lo = center_u - frac * width
        if lo < 0.0:
            lo = 0.0
        if lo + width > 1.0:
            lo = 1.0 - width
        self.view = (max(lo, 0.0), min(lo + width, 1.0))

    def u_at(self, x):
        """The loop coordinate of strip x, in the current view."""
        return u_from_x(x, self.view)
