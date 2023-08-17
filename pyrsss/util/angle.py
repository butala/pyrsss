import numpy as np


def convert_lon(lon):
    """
    Return 0 <= *lon* < 360 converted to -180 <= *lon < 180.
    """
    return lon - 360 if lon > 180 else lon


def deg2tenths_of_arcminute(deg):
    """
    Return *deg* converted to tenths of arcminutes (i.e., arcminutes *
    10).
    """
    return 10 * deg * 60


def rad2tenths_of_arcminutes(rad):
    """
    Return *rad* convert to tenths of arcminutes.
    """
    return deg2tenths_of_minute(np.degrees(rad))
