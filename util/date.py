from __future__ import division

from datetime import datetime, timedelta, time


J2000_EPOCH = datetime(2000, 1, 1, 12)
""" ??? """


def toJ2000(dt):
    """
    Convert *dt* to the number of seconds past J2000.
    """
    return (dt - J2000_EPOCH).total_seconds()


def fromJ2000(seconds):
    """
    Convert *seconds* to a :class:`datetime` object.
    """
    return J2000_EPOCH + timedelta(seconds=seconds)


def ut2lt(dt_ut, lon):
    """
    Convert :class:`datetime` expressed in UT *dt_ut* to the local
    time at longitude *lon* (in [deg]).
    """
    return dt_ut + timedelta(hours=lon * 24 / float(360))


def lt2ut(dt_lt, lon):
    """
    Convert :class:`datetime` expressed in local time *dt_lt* at
    longitude *lon* (in [deg]) to UT.
    """
    return dt_lt - timedelta(hours=lon * 24 / float(360))


def date2dt(date):
    """
    Convert :class:`date` to equivalent :class:`datetime`.
    """
    return datetime.combine(date,
                            time())
