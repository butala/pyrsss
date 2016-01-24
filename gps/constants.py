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
