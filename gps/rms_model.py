import math
from collections import OrderedDict

import numpy as NP


SNR_8 = OrderedDict([#( 5, 46.0),
                     (10, 35.5),
                     (15, 26.5),
                     (20, 19.0),
                     (25, 13.5),
                     (30,  9.5),
                     (35,  7.0),
                     (40,  5.0),
                     (45,  4.5),
                     (50,  4.0),
                     (55,  3.0),
                     (60,  3.0),
                     (65,  3.0),
                     (70,  2.0),
                     (75,  2.0),
                     (80,  2.0),
                     (85,  3.0),
                     (90,  2.0)])


def fit_power_law(x, y):
    """
    """
    ln_x = NP.log(x)
    ln_y = NP.log(y)
    # least squares solution
    A = NP.empty((len(x), 2))
    A[:, 0] = 1
    A[:, 1] = ln_x
    #b_ls = NP.linalg.lstsq(A, ln_y)[0]
    # total least-squares solution
    X = NP.empty((len(x), 3))
    X[:, :2] = A
    X[:, 2] = ln_y
    U, S, V = NP.linalg.svd(X, 1)
    b_tls = (V[-1, :] / -V[-1, -1])[:2]
    alpha = math.exp(b_tls[0])
    beta = b_tls[1]
    return alpha, beta


class RMSModel(OrderedDict):
    def __init__(self, model_spec=SNR_8):
        """
        ???
        """
        super(RMSModel, self).__init__()
        self.el = NP.array(model_spec.keys())
        self.rms = NP.array(model_spec.values())
        self.alpha, self.beta = fit_power_law(self.el,
                                              self.rms)

    def __call__(self, el):
        """
        ???

        note that el is in [deg]
        """
        if el < 0:
            raise RuntimeError('el = {} < 0'.format(el))
        elif el > 90:
            raise RuntimeError('el = {} > 90'.format(el))
        if 0 <= el < min(self.el):
            el = min(self.el)
        return self.alpha * el**self.beta


if __name__ == '__main__':
    rms_model = RMSModel()

    print(rms_model(0))
    print(rms_model(22.2))

    el = NP.linspace(5, 89.5, 129)

    import pylab as PL

    PL.plot(el,
            [rms_model(el_i) for el_i in el],
            marker='o',
            c='b',
            ls='-',
            mec='None')
    PL.xlabel('Elevation angle [deg]')
    PL.ylabel('Modeled RMS [TECU]')
    PL.title('Power Law Fit a * el^b, a = {:.1f}, b = {:.3f}'.format(rms_model.alpha,
                                                                     rms_model.beta))
    PL.xlim(10, 90)
    PL.ylim(0, 50)

    PL.show()
