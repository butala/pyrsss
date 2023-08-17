import numpy as np
import scipy.optimize

from ..gpstk import PyPosition, point


def ipp_from_azel(stn_pos, az, el, ht=450, tol=1e-5):
    """
    Compute the :class:`PyPosition` IPP for line-of-sight specified by
    the :class:`PyPosition` receiver position *stn_pos* and azimuth
    *az* (in degrees) and elevation *el* (in degrees) at height *ht*
    (in [km]). Use tolerance *tol* in the cost function minimization.
    """
    end_pos = point(stn_pos, az, el, 1.1 * ht * 1e3)
    stn_pos_xyz = np.array(stn_pos.xyz)
    end_pos_xyz = np.array(end_pos.xyz)
    def los_pos(s):
        los_xyz = (1 - s) * stn_pos_xyz + s * end_pos_xyz
        return PyPosition(*los_xyz)
    def J(s):
        return abs(los_pos(s).height - ht * 1e3)
    res = scipy.optimize.minimize_scalar(J, bounds=(0, 1), tol=tol)
    assert res.success
    s_star = res.x
    return los_pos(s_star)
