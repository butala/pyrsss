from functools import reduce

import numpy as np
import scipy as sp

from .util import zero_pad
from .convmtx import Convmtx


class SepFilter(object):
    def __init__(self, n, h_list, mode='full'):
        """
        Construct a separable filter, i.e., the filter that operates on a
        signal with shape *n* where the $i$th component of *h_list* is
        the kernel for the $i$th dimension. The *mode* parameter can be
        one of:

        - full: standard convolution, i.e., zero-padding at the edges.

        - valid: convolution where only those portions of complete
          overlap, i.e., no zero-padding, are considered.

        - circ: circular convolution, i.e., periodic boundary
          condition at the edges.
        """
        self.n = n
        self.ndim = len(n)
        self.h_list = [np.copy(h_i) for h_i in h_list]
        self.m = tuple(map(lambda x: len(x), self.h_list))
        self.mode = mode
        self.k_full = tuple(np.array(self.n) + np.array(self.m) - 1)
        self.k_valid = tuple(np.array(self.n) - np.array(self.m) + 1)
        if self.mode == 'full':
            self.k = self.k_full
        elif self.mode == 'valid':
            self.k = self.k_valid
        elif self.mode == 'circ':
            self.k = self.k_full
        else:
            assert False
        self.H_list = list(map(lambda x: sp.fft.fft(x[1],
                                                    n=self.k_full[x[0]]),
                               enumerate(self.h_list)))

        def reducer(x, y):
            i, y_i = y
            shape = [len(y_i)] + [1]*i
            return np.reshape(y_i, shape) * x

        self.h = reduce(reducer, enumerate(reversed(self.h_list)), 1)
        self.H = reduce(reducer, enumerate(reversed(self.H_list)), 1)


    def operate(self, x):
        """
        Apply the separable filter to the signal vector *x*.
        """
        X = sp.fft.fftn(x, s=self.k_full)
        if np.isrealobj(self.h) and np.isrealobj(x):
            y = np.real(sp.fft.ifftn(self.H * X))
        else:
            y = sp.fft.ifftn(self.H * X)

        if self.mode == 'full' or self.mode == 'circ':
            return y
        elif self.mode == 'valid':
            slice_list = []
            for i in range(self.ndim):
                if self.m[i]-1 == 0:
                    slice_list.append(slice(None, None, None))
                else:
                    slice_list.append(slice(self.m[i]-1, -(self.m[i]-1), None))
            return y[*slice_list]
        else:
            assert False


    def __mul__(self, x):
        """
        Apply the separable filter (via the multiplication operator) to
        the signal vector *x*.
        """
        return self.operate(x)


    def asmatrix(self):
        """
        Return the sparse matrix representation of the separable filter.
        """
        h_matrix = np.array([1])
        for i in range(self.ndim):
            if self.mode == 'circ':
                h_i = Convmtx([self.k[i]], self.h_list[i], mode=self.mode)
            else:
                h_i = Convmtx([self.n[i]], self.h_list[i], mode=self.mode)
            h_matrix = sp.sparse.kron(h_matrix, h_i)
        return h_matrix


def random_validation(N,
                      ndim_max=4,
                      n_max=10,
                      m_max=10,
                      modes=['full', 'valid', 'circ']):
    """
    Validate the :class:`SepFilter` implementation by comparing *N*
    direct circular and full / valid convolutions with those computed
    using :class:`SepFilter`. Limit the kernel and signal vector
    dimension to *n_dim_max*, and the length per dimension to *n_max*
    for the signal and *m_max* for the convolution kernel. Only
    consider convolution types given in *modes*.
    """
    for i in range(1, N+1):
        ndim = np.random.randint(1, ndim_max)
        n = np.random.randint(1, n_max, ndim)
        m = np.random.randint(1, m_max, ndim)
        mode = modes[np.random.randint(0, 2)]

        for j, n_j in enumerate(n):
            if n_j - m[j] + 1 <= 0:
                n[j] = m[j]

        if mode == 'full' or mode == 'circ':
            k = m + n - 1
        elif mode == 'valid':
            k = n - m + 1
        else:
            assert False

        assert (np.array(k) > 0).all()

        print('{} of {} (n={} k={} mode={})'.format(i, N, n, k, mode))

        h_list = [np.random.randn(m_i) for m_i in m]

        x = np.random.randn(*n)

        H = SepFilter(n, h_list, mode=mode)

        if mode == 'circ':
            y_true = sp.signal.convolve(x, H.h, mode='full')
        else:
            y_true = sp.signal.convolve(x, H.h, mode=mode)
        y_sep = H * x
        assert np.allclose(y_sep, y_true)

        if mode == 'circ':
            H_matrix_true = Convmtx(k, H.h, mode=mode)
            x = zero_pad(x, k)
        else:
            H_matrix_true = Convmtx(n, H.h, mode=mode)

        H_matrix_sep = H.asmatrix()

        y_matrix_true = H_matrix_true * x.flat
        y_matrix_sep = H_matrix_sep * x.flat
        assert np.allclose(y_matrix_sep, y_matrix_true)

    return True


if __name__ == '__main__':
    random_validation(100)
