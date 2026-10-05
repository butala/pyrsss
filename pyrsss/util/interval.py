"""
Interval helpers (backed by the :mod:`portion` library). ``None``
bounds are interpreted as unbounded (as the retired ``intervals``
package did).
"""
import portion as P


def closed_open(lower, upper):
    """
    Return the half-bounded interval [*lower*, *upper*). ``None`` bounds
    extend to infinity.
    """
    if lower is None:
        lower = -P.inf
    if upper is None:
        upper = P.inf
    return P.closedopen(lower, upper)


def closed(lower, upper):
    """
    Return the closed interval [*lower*, *upper*]. ``None`` bounds
    extend to infinity.
    """
    if lower is None:
        lower = -P.inf
    if upper is None:
        upper = P.inf
    return P.closed(lower, upper)


def length(interval):
    """
    Return the length of *interval* (``float('inf')`` when unbounded).
    """
    try:
        return interval.upper - interval.lower
    except TypeError:
        return float('inf')
