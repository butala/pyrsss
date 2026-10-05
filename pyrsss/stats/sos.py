import math
from abc import ABC, abstractmethod

import numpy as np
import scipy as sp

from .arma import arma_sensitivity


class SOS(ABC):
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

    @classmethod
    def get_mask(cls, Nb, Na):
        mask = np.full((cls.N_sections(Nb, Na), 6), False, dtype=bool)
        mask[:, :3].flat[:Nb] = True
        mask[:, 4:].flat[:Na] = True
        return mask

    # Need to add Nk, zi
    @classmethod
    def zero(cls, Nb, Na):
        sos = np.zeros((cls.N_sections(Nb, Na), 6))
        sos[:, 3] = 1
        return cls(sos, cls.get_mask(Nb, Na))

    # Need to add Nk, zi
    @classmethod
    def I(cls, Nb, Na):
        sos = np.zeros((cls.N_sections(Nb, Na), 6))
        sos[:, [0, 3]] = 1
        return cls(sos, cls.get_mask(Nb, Na))

    @classmethod
    @abstractmethod
    def _theta0(cls, Nb, Na):
        pass

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
    @abstractmethod
    def b(self):
        pass

    @abstractmethod
    def __call__(self, x, axis=-1, zi=None):
        pass

    # Nk, zi
    def residual(self, x, y, axis=-1, zi=None):
        return y - self(x, axis=axis, zi=zi)

    # Nk, zi
    @abstractmethod
    def jacobian(self, u):
        pass

    @classmethod
    def _theta0(cls, Nb, Na):
        return cls.I(Nb, Na)

    @classmethod
    def fit_nonlinear(cls, x, y, Na, Nb, theta0=None, **kwds):
        """
        """
        if theta0 is None:
            theta0 = cls._theta0(Nb, Na).theta
        result = sp.optimize.least_squares(lambda theta: cls.from_theta(theta, Nb, Na).residual(x, y),
                                           theta0,
                                           jac=lambda theta: -cls.from_theta(theta, Nb, Na).jacobian(x),
                                           **kwds)

        if not result.success:
            raise RuntimeError(result.message)
        return cls.from_theta(result.x, Nb, Na)


# aka Cascade
class Series(SOS):
    @property
    def b(self):
        b_accum = self.sos[0, :3]
        for i in range(1, self.sos.shape[0]):
            b_accum = np.convolve(b_accum, self.sos[i, :3])
        return b_accum

    # Nk, zi
    def __call__(self, x, axis=-1, zi=None):
        return sp.signal.sosfilt(self.sos, x, axis=axis, zi=zi)

    @classmethod
    def _theta0(cls, Nb, Na):
        return cls.I(Nb, Na)

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
        return np.c_[*columns]


# aka Concurrent
class Parallel(SOS):
    @property
    def b(self):
        b_accum = np.zeros(3 * self.sos.shape[0] - (self.sos.shape[0] - 1), dtype=complex)
        for i in range(self.sos.shape[0]):
            accum_i = self.sos[i, :3]
            for k in range(self.sos.shape[0]):
                if i == k:
                    continue
                accum_i = np.convolve(accum_i, self.sos[k, 3:])
            b_accum += accum_i
        return b_accum

    # Nk, zi
    def __call__(self, x, axis=-1, zi=None):
        y = np.zeros_like(x, dtype=float if np.isrealobj(self.sos) else complex)
        for i in range(self.sos.shape[0]):
            y += sp.signal.lfilter(self.sos[i, :3], self.sos[i, 3:], x)
        return y

    @classmethod
    def _theta0(cls, Nb, Na):
        return cls.zero(Nb, Na)

    # Nk, zi
    def jacobian(self, u):
        M = len(u)
        columns = []
        for i in range(self.mask.shape[0]):
            J_ab = arma_sensitivity(self.sos[i, :3], self.sos[i, 3:], u, 0)
            # mask=True b parameters
            columns.extend(J_ab[:, 2:][:, self.mask[i, :3]].T)
            # mask=True a parameters
            columns.extend(J_ab[:, :2][:, self.mask[i, 4:]].T)
        return np.c_[*columns]
