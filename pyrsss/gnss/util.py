import math
from datetime import timedelta

from ..util.date import GPS_EPOCH


def shell_mapping(el_deg,
                  h=450,
                  R_mean=6371):
    """
    ???
    """
    el_rad = math.radians(el_deg)
    return 1 / math.sqrt(1 - (R_mean * math.cos(el_rad) / (h + R_mean))**2)


def convert_gps_week(week, seconds=0):
    """
    Convert GPS *week* and *seconds* past the beginning of the week to
    a :class:`datetime`.
    """
    return GPS_EPOCH + timedelta(seconds=week * 60 * 60 * 24 * 7) + timedelta(seconds=seconds)


def get_gps_week(dt):
    """
    Return the GPS week associated with :class:`datetime` *dt*.
    """
    return int((dt - GPS_EPOCH).total_seconds() / timedelta(weeks=1).total_seconds())
