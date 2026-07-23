import numpy as np
import scipy as sp


def get_a_and_zi(theta, Na):
    """Unpack the parameters *theta* where *Na* is the number of denominator polynomial terms and
    the parameters are packed as [a[1:], z_i]. Return the tuple ([1, a[1:]], z_i).
    """
    assert len(theta) == 2*Na
    return np.insert(theta[:Na], 0, 1), theta[Na:]


def ar_residual(x_hat, Na, x, y):
    """Return the residual between the AR filter *x_hat* (the tuple of *Na* AR and *Na* initial
    conditions) for the input *x* and output *y*.
    """
    assert(len(x) == len(y))
    a_hat, zi_hat = get_a_and_zi(x_hat, Na)
    y_hat = sp.signal.lfilter([0.], a_hat, x, zi=zi_hat)[0]
    return y - y_hat


def ar_sensitivity(a, x, zi):
    """
    Return the matrix \\partial y_i(theta) / \\partial theta_j
    where theta = [a, zi] and

    y_i(theta) = lfilter(0, a, x, zi=zi)[i].
    """
    assert np.isclose(a[0], 1)
    assert len(zi) == len(a) - 1

    # SWITCH TO SOS FILTER!
    w = sp.signal.lfilter(1, a, x)
    z = sp.signal.lfilter(-1, a, w)

    w2 = sp.signal.lfilter(1, a, np.pad(zi, (0, len(x)-len(zi))))
    v = z + sp.signal.lfilter(-1, a, w2)

    Da_r = np.zeros(len(a) - 1)
    Da_c = np.pad(v[:-1], (1, 0))
    Da = sp.linalg.toeplitz(Da_c, r=Da_r)

    delta = np.zeros(len(x))
    delta[0] = 1
    Dzi_c = sp.signal.lfilter(1, a, delta)
    Dzi_r = np.zeros(len(zi))
    Dzi = sp.linalg.toeplitz(Dzi_c, r=Dzi_r)

    return -np.c_[Da, Dzi]


def ar_fit_nonlinear(x, y, Na, x_hat0=None, **kwds):
    """Calculate the nonlinear fit between the *Na* (number of poles) AR model with *Na* initial
    condition parameters and the input *x* and output *y*. The optimization starts at *x_hat0* (a
    vector with all 0s when `None`). The output is the tuple of the AR and initial conditions.
    Additional keyword arguments **kwds are passed to `sp.optimize.least_squares`.
    """
    if x_hat0 is None:
        x_hat0 = np.zeros(2*Na)
    assert len(x_hat0) == 2*Na

    def jac_wrapper(x_hat, Na, x, y, *args, **kwds):
        a, zi = get_a_and_zi(x_hat, Na)
        return ar_sensitivity(a, x, zi)

    result = sp.optimize.least_squares(ar_residual,
                                       x_hat0,
                                       jac=jac_wrapper,
                                       args=(Na, x, y),
                                       **kwds)

    if not result.success:
        raise RuntimeError(result.message)
    x_hat = result.x
    return get_a_and_zi(x_hat, Na)
