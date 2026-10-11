"""
Shuttle controls for the 4-D scene: one small, purpose-built transport
for moving through time.

The scene is a function of a loop coordinate ``u`` in ``[0, 1]``; the
shuttle owns *when* the scene is evaluated and nothing else. It steers
the ``Viewer4D`` loop (its timer and temporal actors stay the single
driver) and reads input from three places, all feeding the same clock:

* keyboard -- ``space`` play/pause, arrows (also ``.`` ```,``) step,
  ``Home``/``End`` first/last frame, ``Up``/``Down`` speed, ``+``/``-``
  zoom the time axis
* the time strip (:mod:`pyrsss.solar.timeline`) -- drag to jog, fling
  to coast (momentum decays to a stop)
* the mouse wheel anywhere in the scene -- one notch is one frame
  (``Shift``: ten), i.e. a trackpad's two-finger scroll

The frame grid is uniform in ``u`` (``n_frames`` samples of the window)
so stepping, scrubbing and the clock never disagree about "the current
frame".
"""

import math
import time as _time
from datetime import datetime

from .timeline import hit as strip_hit

# default: cross the whole window in ~20 s at the viewer's 60 FPS timer
DEFAULT_SPEED = 1.0 / 1200.0
SPEED_STEP = 2.0 ** 0.25                    # 25% per key press
MOMENTUM_TAU = 0.7                          # coast decay time [s]
MOMENTUM_TICK = 16                          # ms
FLING_EPS = 2e-3                            # u/s: below this, stop
HINT = ("space: play/pause   arrows: step   wheel: step   "
        "drag strip: jog   fling: coast   +/-: zoom   home/end: jump")


def frame_index(u, n_frames):
    """The grid index nearest the loop coordinate *u* (clamped)."""
    u = min(max(float(u), 0.0), 1.0)
    return int(round(u * (n_frames - 1)))


def frame_u(i, n_frames):
    """The loop coordinate of grid index *i* (clamped)."""
    i = min(max(int(i), 0), n_frames - 1)
    return i / (n_frames - 1)


def frame_time(t0, t1, u):
    """The datetime at loop coordinate *u* over ``[t0, t1]`` (clamped)."""
    u = min(max(float(u), 0.0), 1.0)
    return t0 + (t1 - t0) * u


def format_epoch(t0, t1, u):
    """The scene time at *u*, as the clock shows it."""
    return frame_time(t0, t1, u).strftime('%Y-%m-%d %H:%M:%S')


def momentum_step(u, v, dt, tau=MOMENTUM_TAU):
    """
    Advance a coast of velocity *v* [u/s] by *dt* [s] with exponential
    decay *tau*: returns the new ``(u, v)`` (clamped to the window).
    """
    u = min(max(u + v * dt, 0.0), 1.0)
    v *= math.exp(-dt / tau)
    return u, v


class Shuttle:
    """
    Install the transport on a ``Viewer4D`` scene over the window
    ``[t0, t1]``. The viewer's temporal actors are driven through its
    own ``_update_scene``; nothing else touches their ``update``.
    """

    def __init__(self, viewer, t0, t1, timeline, n_frames=241,
                 speed=DEFAULT_SPEED):
        from pyviz4d.ui import create_text_actor

        self.viewer = viewer
        self.t0, self.t1 = t0, t1
        self.timeline = timeline
        self.n_frames = n_frames
        self.speed = speed

        # the loop coordinate is the viewer's own clock
        viewer._max_time = 1.0
        viewer._loop = True
        viewer.is_playing = False              # start paused: user's choice
        viewer.time_speed = speed

        self.hint_label = create_text_actor(HINT, 0.02, 0.02, font_size=12)
        viewer.ren.AddActor(self.hint_label)

        # input
        iren = viewer.iren
        iren.AddObserver('KeyPressEvent', self._on_key)
        iren.AddObserver('LeftButtonPressEvent', self._on_press)
        iren.AddObserver('MouseMoveEvent', self._on_motion)
        iren.AddObserver('LeftButtonReleaseEvent', self._on_release)
        iren.AddObserver('MouseWheelForwardEvent', self._on_wheel)
        iren.AddObserver('MouseWheelBackwardEvent', self._on_wheel)

        # drag / coast state (one timer and observer for the coast loop)
        self._dragging = False
        self._velocity = 0.0
        self._samples = []
        self._coasting = False
        self._coast_last = None
        iren.CreateRepeatingTimer(MOMENTUM_TICK)
        iren.AddObserver('TimerEvent', self._on_coast)

        # the label refresh rides the scene redraw: wrap, not replace
        self._scene_update = viewer._update_scene
        viewer._update_scene = self._update_scene
        self._update_scene()

    # ---- clock --------------------------------------------------------
    def set_time(self, u, render=True):
        u = min(max(float(u), 0.0), 1.0)
        self.viewer.current_time = u
        if render:
            self._update_scene()

    def step(self, delta):
        """Move *delta* grid frames along the window (and redraw)."""
        i = frame_index(self.viewer.current_time, self.n_frames) + delta
        self.set_time(frame_u(i, self.n_frames))

    def set_speed(self, speed):
        self.speed = max(min(float(speed), 1.0 / 30.0), 1.0 / 1e5)
        self.viewer.time_speed = self.speed

    def toggle(self):
        self.viewer.is_playing = not self.viewer.is_playing
        self._update_scene()

    def _update_scene(self):
        self._scene_update()

    # ---- mouse --------------------------------------------------------
    def _norm(self, iren):
        w, h = iren.GetSize()
        x, y = iren.GetEventPosition()
        return x / max(w - 1, 1), y / max(h - 1, 1)

    def _on_press(self, obj, _event):
        x, y = self._norm(obj)
        if not strip_hit(x, y):
            return
        self._coasting = False
        self._dragging = True
        self._samples = [(_time.monotonic(), x)]
        self.set_time(self.timeline.u_at(x))

    def _on_motion(self, obj, _event):
        if not self._dragging:
            return
        x, _ = self._norm(obj)
        now = _time.monotonic()
        self._samples.append((now, x))
        self._samples = self._samples[-4:]
        self.set_time(self.timeline.u_at(x))

    def _on_release(self, obj, _event):
        if not self._dragging:
            return
        self._dragging = False
        # fling: velocity from the last drag samples, in u/s
        if len(self._samples) >= 2:
            (t1, x1), (t2, x2) = self._samples[-2], self._samples[-1]
            dt = t2 - t1
            if dt > 1e-4:
                self._velocity = (self.timeline.u_at(x2)
                                  - self.timeline.u_at(x1)) / dt
        if abs(self._velocity) > FLING_EPS:
            self._coasting = True
            self._coast_last = _time.monotonic()
        else:
            self._velocity = 0.0

    def _on_coast(self, _obj, _event):
        if not self._coasting:
            return
        now = _time.monotonic()
        dt, self._coast_last = now - self._coast_last, now
        u, self._velocity = momentum_step(self.viewer.current_time,
                                          self._velocity, dt)
        self.set_time(u)
        if abs(self._velocity) < FLING_EPS:
            self._coasting = False

    def _on_wheel(self, obj, event):
        step = 10 if obj.GetShiftKey() else 1
        self.step(+step if event == 'MouseWheelForwardEvent' else -step)

    # ---- keyboard -----------------------------------------------------
    def _on_key(self, obj, _event):
        key = obj.GetKeySym()
        if key in ('Right', 'period'):
            self.step(+1)
        elif key in ('Left', 'comma'):
            self.step(-1)
        elif key == 'Home':
            self.set_time(0.0)
        elif key == 'End':
            self.set_time(1.0)
        elif key == 'space':
            self.toggle()
        elif key == 'Up':
            self.set_speed(self.speed * SPEED_STEP)
        elif key == 'Down':
            self.set_speed(self.speed / SPEED_STEP)
        elif key in ('plus', 'equal', 'KP_Add'):
            self.timeline.zoom(1 / 1.4)
            self._update_scene()
        elif key in ('minus', 'KP_Subtract'):
            self.timeline.zoom(1.4)
            self._update_scene()
