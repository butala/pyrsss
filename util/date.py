from datetime import datetime, timedelta


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
