import numpy as np
import scipy.signal

from .sepfilter import SepFilter


def differentiator(n, Hz=1):
    """
    Return linear phase impulse response for a length *n* filter that
    approximates the differential operator. The sampling frequency is
    *Hz*.

    The filter length *n* must be even. The remez function returns
    type 3 (for *n* odd) and 4 (for *n* even) linear phase
    filters. However, type 3 linear phase filters are 0 at $\omega =
    0$ and $\omega = \pi$.
    """
    if n % 2 == 1:
        raise ValueError('the filter length n must be even')
    return scipy.signal.remez(n,
                              [0, Hz / 2],
                              [1],
                              Hz=Hz,
                              type='differentiator') * Hz * 2 * np.pi


class GradientFilter(SepFilter):
    def __init__(self, n, axis, order=2, mode='valid'):
        """
        Construct gradient operator for signal of dimension *n* for
        dimension *axis*. Use a filter kernel of length *order* (must
        be even). Use convolution type *mode*.
        """
        self.n = n
        self.ndim = len(self.n)
        self.axis = axis
        if axis < 0 or axis >= self.ndim:
            raise ValueError('0 <= axis (= {0}) < ndim = {1}'.format(axis, self.ndim))

        self.d = differentiator(order)

        h_list = []
        for i in reversed(range(self.ndim)):
            if i == axis:
                h_list.append(self.d)
            else:
                h_list.append(np.array([1]))

        super(GradientFilter, self).__init__(n, h_list, mode=mode)


class Gradient2Filter(SepFilter):
    def __init__(self, n, axis, order=3, mode='valid'):
        """
        Construct a second-order gradient operator for signal of dimension
        *n* for dimension *axis*. Use a filter kernel of length
        *order* (must be odd). Use convolution type *mode*.
        """
        # assert that the filter length is odd
        assert(order % 2 == 1)
        self.n = n
        self.ndim = len(self.n)
        self.axis = axis
        if axis < 0 or axis >= self.ndim:
            raise ValueError('0 <= axis (= {0}) < ndim = {1}'.format(axis, self.ndim))

        self.d = differentiator(int(order/2) + 1)
        self.d2 = np.convolve(self.d, self.d)

        self.mode = mode

        h_list = []
        m = []
        for i in reversed(range(self.ndim)):
            if i == axis:
                h_list.append(self.d2)
            else:
                h_list.append(np.array([1]))
            m.append(len(h_list[-1]))
        self.m = m

        if mode == 'circ':
            n_prime = array(n) - m + 1
            super(Gradient2Filter, self).__init__(n_prime, h_list, mode=mode)
        else:
            super(Gradient2Filter, self).__init__(n, h_list, mode=mode)


def plot_diff_spectra(l=[2, 4, 10]):
    """Plot the magnitude response of differentiators of length *l* and
    include the response of the ideal differentiator. Each length in
    *l* must be even.
    """
    from collections import OrderedDict

    import numpy as np
    import matplotlib.pyplot as plt

    from pyrsss.signal.spectrum import spectrum

    h_map = OrderedDict()
    for l_i in l:
        h_map[l_i] = differentiator(l_i)

    oversample = 8
    H_map = OrderedDict()
    for l_i, h_i in h_map.items():
        H_map[l_i] = spectrum(h_i, oversample=oversample)

    f_ideal = np.linspace(0, 0.5, 129)
    H_ideal = np.abs(f_ideal) * 2 * np.pi

    plt.figure(figsize=(8, 4))
    plt.plot(f_ideal,
             H_ideal,
             c='C0',
             zorder=10,
             label='l=$\infty$')
    for i, (l_i, (f_i, H_i)) in enumerate(H_map.items(), 1):
        plt.plot(f_i,
                 np.abs(H_i),
                 c='C{}'.format(i),
                 label='l={}'.format(l_i))
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.title('Comparison of Differentiator Magnitude Responses')

    plt.show()


if __name__ == '__main__':
    plot_diff_spectra()
