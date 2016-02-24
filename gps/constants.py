import math
from datetime import datetime

import scipy.constants as const


EPOCH = datetime(1980, 1, 6)
"""The epoch for GPS time."""

F_0 = 10.23e6
"""Fundamental GPS frequency [Hz]."""

N_1 = 154
"""GPS carrier 1 frequency factor."""

N_2 = 120
"""GPS carrier 2 frequency factor."""

F_1 = N_1 * F_0
"""GPS carrier 1 frequency [Hz]."""

F_2 = N_2 * F_0
"""GPS carrier 2 frequency [Hz]."""

LAMBDA_1 = const.c / F_1
"""GPS carrier wavelength 1 [m]."""

LAMBDA_2 = const.c / F_2
"""GPS carrier 2 wavelength [m]."""

K = const.e**2 / (8 * math.pi**2 * const.epsilon_0 * const.m_e)
"""Constant used in numerious calculation [m**3 / s**2]. It is related to the
plasma frequency."""

TECU_TO_NS = K * 10**25 / const.c * (N_1**2 - N_2**2) / (F_0 * N_1 * N_2)**2
"""Conversion from [TECU] to differential delay [ns]."""

TECU_TO_KM = TECU_TO_NS * (const.c / (1e9 * 1e3))
"""Conversion factor from [TECU] to [km]."""

TECU_TO_M = TECU_TO_KM * 1e3
"""Conversion factor from [TECU] to [m]."""

M_TO_TECU = 1 / TECU_TO_M
"""Conversion factor from [m] to [TECU]."""

SHELL_HEIGHT = 450
"""Thin-shell model ionosphere height [km]."""
