from __future__ import division

import numpy as NP


def chapman(z, Nm, Hm, H_O):
    """
    Return Chapman function electron density at height *z*, maximum
    electron density at the F-peak *Nm*, height at the maximum *Hm*,
    and the scale height of atomic oxygen *H_O*.
    """
    return Nm * NP.exp((1 - ((z - Hm) / H_O) - NP.exp(-(z - Hm) / H_O)) / 2)
