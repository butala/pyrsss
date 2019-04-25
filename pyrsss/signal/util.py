import numpy as NP


def zero_pad(x, n):
    """
    Return *x* zero padded to the dimensions specified in *n*.
    """
    assert len(n) == x.ndim
    return NP.pad(x,
                  [(0, n_i - s_i) for n_i, s_i in zip(n, x.shape)],
                  'constant')
