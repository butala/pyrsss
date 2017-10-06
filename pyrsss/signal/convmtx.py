from __future__ import division

import numpy as NP
import scipy.linalg
import scipy.signal
import scipy.sparse as sparse

from util import zero_pad


class Convmtx(sparse.coo_matrix):
    def __new__(cls, n, H, mode='full'):
        """
        Construct sparse convolution matrix to operate on vector of
        dimension *n* with the kernel *H*. The *mode* parameter can be
        one of:

        - full: standard convolution, i.e., zero-padding at the edges.

        - valid: convolution where only those portions of complete
          overlap, i.e., no zero-padding, are considered.

        - circ: circular convolution, i.e., periodic boundary
          condition at the edges.
        """
        def toeplitz_mapper_full(h):
            if (h == 0).all():
                return sparse.coo_matrix((k[-1], n[-1]))
            else:
                c = h
                r = NP.array([c[0]] + [0]*(n[-1]-1))
                return sparse.coo_matrix(scipy.linalg.toeplitz(c, r))

        def toeplitz_mapper_valid(h):
            if (h == 0).all():
                return sparse.coo_matrix((k[-1], n[-1]))
            else:
                r = NP.zeros(n[-1])
                r[:len(h)] = h
                c = NP.zeros(k[-1])
                c[0] = r[0]
                return sparse.coo_matrix(scipy.linalg.toeplitz(c, r))

        def toeplitz_mapper_circ(h):
            if (h == 0).all():
                return sparse.coo_matrix((k[-1], n[-1]))
            else:
                c = h
                r = NP.zeros(n[-1])
                r[0] = c[0]
                r[1:] = h[:0:-1]
                return sparse.coo_matrix(scipy.linalg.toeplitz(c, r))

        def block_mapper_full(n, k, blocks):
            c = [blocks[i] for i in range(k)]
            r = [c[0]] + [None]*(n-1)
            return sparse.bmat(scipy.linalg.toeplitz(c, r).tolist(), format='coo')

        def block_mapper_valid(n, k, blocks):
            r = []
            for i in range(n):
                if (n - k - i < 0):
                    r.append(None)
                else:
                    r.append(blocks[n - k - i])

            c = []
            for i in range(n-k, n):
                c.append(blocks[i])

            return sparse.bmat(scipy.linalg.toeplitz(c, r).tolist(), format='coo')

        def block_mapper_circ(n, k, blocks):
            c = [blocks[i] for i in range(k)]
            r = []
            r.append(blocks[0])
            r.extend(blocks[:0:-1])
            return sparse.bmat(scipy.linalg.toeplitz(c, r).tolist(), format='coo')

        m = H.shape

        if mode == 'full':
            k = tuple(NP.array(n) + NP.array(m) - 1)
            toeplitz_mapper = toeplitz_mapper_full
            block_mapper = block_mapper_full

            H_zp = zero_pad(H, k)
            c_list = NP.split(H_zp.flatten(), NP.prod(k[:-1]))
        elif mode == 'valid':
            k = tuple(NP.array(n) - NP.array(m) + 1)
            toeplitz_mapper = toeplitz_mapper_valid
            block_mapper = block_mapper_valid

            H_zp = zero_pad(H[...,::-1], n)
            c_list = NP.split(H_zp.flatten(), NP.prod(n[:-1]))
        elif mode == 'circ':
            assert((NP.array(m) <= NP.array(n)).all())
            k = n
            toeplitz_mapper = toeplitz_mapper_circ
            block_mapper = block_mapper_circ

            H_zp = zero_pad(H, k)
            c_list = NP.split(H_zp.flatten(), NP.prod(k[:-1]))
        else:
            raise ValueError('Unknown mode {0}'.format(mode))

        blocks = map(toeplitz_mapper, c_list)

        for n_i, k_i in zip(n[-2::-1], k[-2::-1]):
            if mode == 'full' or mode == 'circ':
                blocks = map(lambda x: block_mapper(n_i, k_i, x),
                             NP.split(NP.array(blocks), len(blocks)/k_i))
            elif mode =='valid':
                blocks = map(lambda x: block_mapper(n_i, k_i, x),
                             NP.split(NP.array(blocks), len(blocks)/n_i))
            else:
                raise ValueError('Unknown mode {0}'.format(mode))

        return blocks[0]


def ndcircconv(x, h):
    """
    Compute the circular convolution of real, n-dimensional vectors
    *x* and *h*.
    """
    n = x.shape
    m = h.shape
    k = NP.array(n) + NP.array(m) - 1
    return NP.real(NP.fft.ifftn(NP.fft.fftn(h, s=k) * NP.fft.fftn(x, s=k))).flat


def random_validation(N,
                      n_dim_max=4,
                      n_max=10,
                      m_max=10):
    """
    Validate the :class:`Convmtx` implementation by comparing *N*
    direct circular and full / valid convolutions with those computed
    using :class:`Convmtx`. Limit the kernel and signal vector
    dimension to *n_dim_max*, and the length per dimension to *n_max*
    for the signal and *m_max* for the convolution kernel.
    """
    print('Testing circ mode')
    for i in range(1, N + 1):
        n_dim = NP.random.random_integers(1, n_dim_max)
        n = NP.random.random_integers(1, n_max, n_dim)
        m = NP.random.random_integers(1, m_max, n_dim)
        x = NP.random.randn(*n)
        h = NP.arange(NP.prod(m)) + 1
        h.shape = m
        k = NP.array(n) + NP.array(m) - 1
        print('{} of {} (n={} k={} mode=circ)'.format(i, N, n, k))
        H = Convmtx(k, h, mode='circ')
        y_true = ndcircconv(x, h)
        x_zp = zero_pad(x, k)
        y_mtx = H * x_zp.flat
        assert(NP.allclose(y_true, y_mtx))
    print('')
    print('Testing full and valid modes')
    mode_list = ['full', 'valid']
    for i in range(1, N + 1):
        n_dim = NP.random.random_integers(1, n_dim_max)
        n = NP.random.random_integers(1, n_max, n_dim)
        m = NP.random.random_integers(1, m_max, n_dim)

        mode = mode_list[NP.random.random_integers(0, 1)]
        if mode == 'full':
            k = tuple(NP.array(n) + NP.array(m) - 1)
        elif mode == 'valid':
            k = tuple(NP.array(n) - NP.array(m) + 1)
        else:
            assert(False)

        x = NP.random.randn(*n)
        h = NP.arange(NP.prod(m)) + 1
        h.shape = m

        if (NP.array(k) <= 0).any():
            assert(mode == 'valid')
            n, m = m, n
            x, h = h, x
            k = tuple(NP.array(n) - NP.array(m) + 1)

        print('{} of {} (n={} k={} mode={})'.format(i, N, n, k, mode))

        if (NP.array(k) > 0).all():
            H = Convmtx(n, h, mode=mode)
            y_true = scipy.signal.convolve(x, h, mode=mode)
            y_mtx = H * x.flat
            y_mtx.shape = k

            assert(NP.allclose(y_true, y_mtx))
        else:
            try:
                y_true = scipy.signal.convolve(x, h, mode=mode)
            except ValueError:
                pass
            else:
                assert(NP.prod(y_true.shape) == 0)
    return True


if __name__ == '__main__':
    random_validation(100)
