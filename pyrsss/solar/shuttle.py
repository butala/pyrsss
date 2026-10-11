"""
Shuttle controls for the 4-D scene: a small, purpose-built transport
for moving through time.

One idea: the scene is a function of a loop coordinate ``u`` in
``[0, 1]``; the shuttle owns *when* the scene is evaluated and nothing
else. It steers the ``Viewer4D`` loop (its timer and temporal actors
stay the single driver) and adds what a transport needs:

* ``space`` -- play / pause
* ``Right`` / ``Left`` (also ``.`` / ``,``) -- step one frame
* ``Home`` / ``End`` -- first / last frame
* ``Up`` / ``Down`` -- playback speed
* a thin scrub slider and an epoch label showing the scene time

The frame grid is uniform in ``u`` (``n_frames`` samples of the window);
stepping moves exactly one grid position, so the display and the slider
never disagree about where "the current frame" is.
"""

from datetime import datetime

# default: cross the whole window in ~20 s at the viewer's 60 FPS timer
DEFAULT_SPEED = 1.0 / 1200.0
SPEED_STEP = 2.0 ** 0.25                    # 25% per key press
HINT = ("space: play/pause    arrows: step    "
        "home/end: first/last    up/down: speed")


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
    """The scene time at *u*, as the label shows it."""
    return frame_time(t0, t1, u).strftime('%Y-%m-%d %H:%M:%S')


class Shuttle:
    """
    Install the transport on a ``Viewer4D`` scene over the window
    ``[t0, t1]``. The viewer's temporal actors are driven through its
    own ``_update_scene``; nothing else touches their ``update``.
    """

    def __init__(self, viewer, t0, t1, n_frames=241, speed=DEFAULT_SPEED):
        from pyviz4d.ui import create_slider, create_text_actor

        self.viewer = viewer
        self.t0, self.t1 = t0, t1
        self.n_frames = n_frames
        self.speed = speed

        # the loop coordinate is the viewer's own clock
        viewer._max_time = 1.0
        viewer._loop = True
        viewer.is_playing = False              # start paused: user's choice
        viewer.time_speed = speed

        self.slider = create_slider(viewer.iren, 'Time', 0.0, 1.0, 0.0,
                                    0.1, self._on_slider)
        self.hint_label = create_text_actor(HINT, 0.02, 0.02, font_size=12)
        viewer.ren.AddActor(self.hint_label)

        viewer.iren.AddObserver('KeyPressEvent', self._on_key)
        # refresh the label whenever the scene redraws: wrap, not replace
        self._scene_update = viewer._update_scene
        viewer._update_scene = self._update_scene
        self._update_scene()

    # ---- coordinate helpers -------------------------------------------
    def format(self, u):
        return format_epoch(self.t0, self.t1, u)

    def _apply(self, u, render=True):
        self.viewer.current_time = u
        self.slider.GetRepresentation().SetValue(u)
        if render:
            self._update_scene()

    def _update_scene(self):
        self._scene_update()

    def _on_slider(self, slider, _event):
        self._apply(slider.GetRepresentation().GetValue())

    # ---- transport ----------------------------------------------------
    def step(self, delta):
        """Move *delta* grid frames along the window (and redraw)."""
        i = frame_index(self.viewer.current_time, self.n_frames) + delta
        self._apply(frame_u(i, self.n_frames))

    def _on_key(self, obj, _event):
        key = obj.GetKeySym()
        if key in ('Right', 'period'):
            self.step(+1)
        elif key in ('Left', 'comma'):
            self.step(-1)
        elif key == 'Home':
            self._apply(0.0)
        elif key == 'End':
            self._apply(1.0)
        elif key == 'space':
            self.viewer.is_playing = not self.viewer.is_playing
            self._update_scene()
        elif key == 'Up':
            self.speed = min(self.speed * SPEED_STEP, 1.0 / 30.0)
            self.viewer.time_speed = self.speed
        elif key == 'Down':
            self.speed = max(self.speed / SPEED_STEP, 1.0 / 1e5)
            self.viewer.time_speed = self.speed
