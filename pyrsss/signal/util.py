import numpy as NP


def zero_pad(x, n):
    """
    Return *x* zero padded to the dimensions specified in *n*.
    """
    y = NP.zeros(n)
    for index in NP.ndindex(x.shape):
        y[index] = x[index]
    return y
