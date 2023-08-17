import logging
from itertools import product

import numpy as np
from scipy.misc import comb
import sympy as sym

logger = logging.getLogger('pyrsss.util.lagrange')


def poly_power_seq(m, n):
    """
    Return the list of length *m* tuples with sum less than or equal
    to *n*.
    """
    return filter(lambda x: sum(x) <= n,
                  product(*[range(n, -1, -1) for i in range(m)]))


def lagrange_interpolator(points, delta, delta_i_list):
    """
    Construct the multidimensional Lagrange interpolating polynomial
    defined by *delta* and *delta_i_list* at *points*.

    This is (7) in the reference

    Kamron Saniee, A Simple Expression for Multivariate Lagrange
    Interpolation, SIAM, 2007.
    """
    f = 0
    for i, delta_i in enumerate(delta_i_list):
        f += points[i][-1] * delta_i / delta
    # fix for corner case when polynomial reduces to a single float
    if isinstance(f, float):
        return sym.Float(f)
    return f


def multivariate_lagrange(points, n):
    """
    Given the list of values *points*, construct a multidimensional
    (with dimension equal to the length of each of the given points)
    polynomial of degree *n* using the method developed in the
    following reference:

    Kamron Saniee, A Simple Expression for Multivariate Lagrange
    Interpolation, SIAM, 2007.

    The method has a restriction: if p is the number of *point*s and m
    is one less than the problem domain (i.e., m = `len(points[0]) -
    1`) then p = *n* + m choose *n*.
    """

    # check dimensions
    p = len(points)
    m = len(points[0]) - 1

    if not comb(n + m, n, exact=True) == p:
        raise ValueError('dimension mismatch')

    # setup symbols and build expression
    z_vec = sym.Matrix([points[i][-1] for i in range(p)])
    coef = sym.symbols(' '.join(['a{}'.format(i) for i in range(1, p + 1)]))
    var = sym.symbols(' '.join(['x{}'.format(i) for i in range(1, m + 1)]))
    # handle trivial case, n = 0
    if n == 0:
        return sym.Float(points[0][-1]), var, 1., [1.]
    z_terms = [sym.Poly.from_dict({x: 1}, var) for x in poly_power_seq(m, n)]
    z = sym.Poly.from_dict(dict(zip(poly_power_seq(m, n), coef)), var)
    # build M matrix
    M_rows = []
    for point_i in points:
        z_i = z(*point_i[:-1])
        M_rows.append([float(z_i.coeff(x)) for x in coef])
    M = np.array(M_rows)
    delta = np.linalg.det(M)
    # compute delta_i
    delta_i_list = []
    B = np.ones(p - 1)
    C = np.array(z_terms[:-1])
    D = M[-1, -1]
    for i in range(p):
        # using block matrix property of determinants which assumes
        # that the matrix A below is non-singular (the typical case)
        # --- the code fails safe to calculating the symbolic
        # determinant when A is singular (which will be slow for large
        # problems)
        A = np.vstack((M[ :i,  :-1],
                       M[i+1:, :-1]))
        try:
            b = np.linalg.solve(A, B)
            # using row interchange property
            delta_i_list.append((-1)**(p-1-i) * np.linalg.det(A) * (D - np.dot(C, b)))
        except np.linalg.linalg.LinAlgError:
            logger.warning('singular matrix encountered (this will happen when, e.g., points contain many zeros) --- resorting to symbolic determinant calculation which will be slow for large problems')
            row_i = sym.Matrix([z_terms])
            M_i = sym.Matrix(M_rows)
            M_i.row_del(i)
            M_i = M_i.row_insert(p, row_i)
            # using row interchange property
            delta_i_list.append(M_i.det() * (-1)**(p-1-i))
    f = lagrange_interpolator(points, delta, delta_i_list)
    return sym.Poly(f.simplify()), var, delta, delta_i_list


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    # check corner case
    # answer = 1
    points0 = [(1, 1, 1)]
    n0 = 0
    poly0, _, _, _ = multivariate_lagrange(points0, n0)
    print(poly0)
    print('')

    # Example 1 from the paper
    # answer = x1 + x2 + 1
    points1 = [(0, 0, 1),
               (0, 1, 2),
               (1, 1, 3)]
    n1 = 1
    poly1, _, _, _ = multivariate_lagrange(points1, n1)
    print(poly1)
    print('')

    # Example 1 from the paper
    # answer = 2 * x1**2 + x1 * x2 - 4 * x2 - 3
    points2 = [( 0,  1,  -7),
               ( 2,  1,   3),
               ( 1,  3, -10),
               (-2, -1,  11),
               (-3,  2,   1),
               (-1,  2, -11)]
    n2 = 2
    poly2, _, _, _ = multivariate_lagrange(points2, n2)
    print(poly2)
    print('')

    # Random large test
    # p_large = 36
    # n_large = 7

    p_large = 66
    n_large = 10

    points_large = [np.random.randn(3) for i in range(p_large)]

    poly_large, _, _, _ = multivariate_lagrange(points_large, n_large)
    print(poly_large)
