import math

import numpy as np
import scipy as sp


class SOS:
    def __init__(self, sos, mask):
        """
        """
        assert sos.shape[1] == 6
        assert np.allclose(sos[:, 3], 1)
        assert sos.shape == mask.shape
        self.sos = sos
        self.mask = mask

    @staticmethod
    def N_sections(Nb, Na):
        return max(math.ceil(Na / 2), math.ceil(Nb / 3))

    @staticmethod
    def get_mask(Nb, Na):
        mask = np.full((SOS.N_sections(Nb, Na), 6), False, dtype=bool)
        mask[:, :3].flat[:Nb] = True
        mask[:, 4:].flat[:Na] = True
        return mask

    @classmethod
    def I(cls, Nb, Na):
        # Need to add Nk, zi
        sos = np.zeros((SOS.N_sections(Nb, Na), 6))
        sos[:, [0, 3]] = 1
        return cls(sos, SOS.get_mask(Nb, Na))

    @classmethod
    def from_theta(cls, theta, Nb, Na):
        assert len(theta) == Na + Nb
        obj = cls.I(Nb, Na)
        obj.sos[obj.mask] = theta
        return obj

    @property
    def theta(self):
        return self.sos[self.mask]

    @property
    def Nb(self):
        # item() returns native Python type
        return np.sum(self.mask[:, :3]).item()

    @property
    def Na(self):
        # item() returns native Python type
        return np.sum(self.mask[:, 4:]).item()

    @property
    def a(self):
        a_accum = self.sos[0, 3:]
        for i in range(1, self.sos.shape[0]):
            a_accum = np.convolve(a_accum, self.sos[i, 3:])
        return a_accum

    @property
    def b(self):
        b_accum = self.sos[0, :3]
        for i in range(1, self.sos.shape[0]):
            b_accum = np.convolve(b_accum, self.sos[i, :3])
        return b_accum

    # Nk, zi
    def __call__(self, x, axis=-1, zi=None):
        return sp.signal.sosfilt(self.sos, x, axis=axis, zi=zi)

    # Nk, zi
    def residual(self, x, y, axis=-1, zi=None):
        return y - self(x, axis=axis, zi=zi)

    # Nk, zi
    def jacobian(self, u):
        M = len(u)
        y = sp.signal.sosfilt(self.sos, u)
        columns = []
        for i in range(self.mask.shape[0]):
            if self.mask.shape[0] == 1:
                # Single element cascade => H1(z) = 1, i.e., the identity system
                a2 = self.sos[0, 3:]
                Sy_b = sp.signal.lfilter(1, a2, u)
            else:
                sos1 = np.r_[self.sos[:i, :], self.sos[i+1:, :]]
                sos2 = self.sos[i, :]
                y1 = sp.signal.sosfilt(sos1, u)
                a2 = sos2[3:]
                Sy_b = sp.signal.lfilter(1, a2, y1)
            Sy_a = sp.signal.lfilter(1, a2, -y)
            for j in np.nonzero(self.mask[i, :])[0]:
                match j:
                    case 0 | 1 | 2:
                        # numerator polynomial coefficient
                        columns.append(np.pad(Sy_b[:(M-j)], (j, 0)))
                    case 4 | 5:
                        # denominator polynomial coefficient
                        k = j - 3  # k = 1 or 2
                        columns.append(np.pad(Sy_a[:(M-k)], (k, 0)))
                    case 3:
                        raise RuntimeError('The first denominator polynomial coefficient is fixed to a0=1')
                    case _:
                        raise RuntimeError('Impossible')
        return -np.c_[*columns]


def fit_nonlinear(x, y, Na, Nb, theta0=None, **kwds):
    """
    """
    if theta0 is None:
        theta0 = SOS.I(Nb, Na).theta
    result = sp.optimize.least_squares(lambda theta: SOS.from_theta(theta, Nb, Na).residual(x, y),
                                       theta0,
                                       jac=lambda theta: SOS.from_theta(theta, Nb, Na).jacobian(x),
                                       **kwds)

    if not result.success:
        raise RuntimeError(result.message)
    return SOS.from_theta(result.x, Nb, Na)
