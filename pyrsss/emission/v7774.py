"""
The following (ALPHA1 through K3) are from Table 2.3 from Makela's
dissertation.
"""

ALPHA1 = 7.8e-13
"""[cm^3 / s]"""


BETA_7774 = 0.42
""""""


K1 = 1.3e-15
"""[cm^3 / s]"""


K2 = 1.5e-7
"""[cm^3 / s]"""


K3 = 1.4e-10
"""[cm^3 / s]"""


def emission_v7774_rr(ne, o):
    """
    """
    return ALPHA1 * ne**2


def emission_v7774_ii(ne, o, o_plus):
    """
    """
    return (BETA_7774 * K1 * K2 * o * o_plus * ne) / (K2 * o_plus + K3 * o)


def emission_v7774(ne, o, o_plus):
    """
    """
    return emission_v7774_rr(ne, o) + emission_v7774_ii(ne, o, o_plus)
