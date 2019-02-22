import logging
from itertools import combinations
from collections import namedtuple

import numpy as NP
import scipy.stats

logger = logging.getLogger('pyrsss.stats.linregress')


"""Linear regression assessment metrics""",
Assessment = namedtuple('Assessment',
                        'std_error ' \
                        't_stat '\
                        'p_value '\
                        'RSE '\
                        'R2 '\
                        'F_stat '\
                        'Cp '\
                        'AIC '\
                        'BIC '\
                        'adjusted_R2')


def linregress_assessment(y, X, beta_hat):
    """
    Return an ss:`Assessment` of the linear regression *beta_hat* for
    the model specified by *y* and *X* (see :func:`linregress`)
    composed of several metrics.
    """
    n, p_plus_1 = X.shape
    p = p_plus_1 - 1
    assert len(y) == n
    RSS = NP.sum((y - NP.dot(X, beta_hat))**2)
    RSE = NP.sqrt(RSS / (n - p - 1))
    std_error = NP.sqrt(NP.diag(RSE**2 * NP.linalg.inv(NP.dot(X.T, X))))
    t_stat = beta_hat / std_error
    p_value = [2*scipy.stats.t.sf(NP.abs(t_stat_i), n - 2) for t_stat_i in t_stat]
    TSS = NP.sum((y - NP.mean(y))**2)
    R2 = 1 - RSS / TSS
    F_stat = (TSS - RSS) / p / (RSS / (n - p - 1))
    sigma_hat = RSE
    d = p
    Cp = (RSS + 2 * d * sigma_hat**2) / n
    AIC = (RSS + 2 * d * sigma_hat**2) / (n * sigma_hat**2)
    BIC = (RSS + NP.log(n) * d * sigma_hat**2) / n
    adjusted_R2 = 1 - (RSS / (n - d - 1)) / (TSS / (n - 1))
    return Assessment(std_error,
                      t_stat,
                      p_value,
                      RSE,
                      R2,
                      F_stat,
                      Cp,
                      AIC,
                      BIC,
                      adjusted_R2)


def linregress(df, y_column):
    """
    Compute the linear regression given the model

    .. math::
    \mathbf{y} &= \mathbf{X}\,\boldsymbol{\beta} + \boldsymbol{\epsilon} \\
    \boldsymbol{\epsilon} &\overset{\text{i.i.d.}}{\sim} \mathcal{N}\bigl(\mathbf{0},\, \sigma^2\,\mathbf{I}\bigr)

    where :math:`\boldsymbol{\beta}` are the :math:`p` parameters,
    :math:`y` is the length :math:`n` column in *df* with identifier
    *y_column*, and the :math:`n \times p` matrix
    :math:`\boldsymbol{X}` is constructed from the remaining columns
    in *df*.

    Return the tuple with the estimated parameters in the first
    element and an :class:`Assessment` (see
    :func:`linregress_assessment`) as the second.
    """
    y = df[y_column].values
    n = len(y)
    cols = [NP.ones(n)]
    for col in df.columns:
        if col == y_column:
            continue
        cols.append(df[col].values)
    X = NP.column_stack(cols)
    p = X.shape[1] - 1
    # compute least-squares fit
    beta_hat, delme, _, _ = NP.linalg.lstsq(X, y)
    return beta_hat, linregress_assessment(y, X, beta_hat)


def best_subset_selection(df, y_column):
    """
    Compute the exhaustive multiple linear regression for the
    combination of all parameters of the model *df* with respect to
    *y_column* (see :func:`linregress`). Return the tuple of lists
    with 1) the best parameter first, 2) the :class:`Assessment` for
    each regression, and 3) the *df* columns for the best fit.
    """
    predictor_cols = list(df.columns)
    predictor_cols.remove(y_column)
    beta_hat_max = []
    assessment_max = []
    R2_max_col = []
    for k in range(1, len(predictor_cols) + 1):
        logger.info('processing {}/{} predictors'.format(k, len(predictor_cols)))
        beta_hat_max.append(None)
        assessment_max.append(None)
        R2_max_col.append(None)
        for col_combo in combinations(predictor_cols, k):
            df_fit = df[list(col_combo) + [y_column]]
            beta_hat, assessment = linregress(df_fit, y_column)
            if (assessment_max[-1] is None) or \
               (assessment.R2 > assessment_max[-1].R2):
                beta_hat_max[-1] = beta_hat
                assessment_max[-1] = assessment
                R2_max_col[-1] = list(col_combo)
    return beta_hat_max, assessment_max, R2_max_col
