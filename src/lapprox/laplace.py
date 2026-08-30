"""Core Laplace-approximation routines.

This module provides a numerical Hessian estimator and the Laplace
approximation of integrals that can be written as ``Z = int exp(f(x)) dx``.
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)


def numerical_second_partial(func, x0, dim1, dim2, **kwargs):
    """Numerically estimate a second partial derivative of ``func``.

    Args:
        func (callable): function to differentiate, taking an ``N``-dimensional
            point as its first argument.
        x0 (array): ``N``-dimensional point at which to evaluate the derivative.
        dim1 (int): dimension of the first partial derivative.
        dim2 (int): dimension of the second partial derivative.

    Keyword Args:
        eps (float): finite-difference step scale (default ``1e-5``).
        priors (dict): optional ``{name: (kind, a, b)}`` mapping used to scale the
            step per parameter (``kind`` in ``"Gaussian"``, ``"Uniform"``,
            ``"Jeffreys"``). Any other keyword arguments are forwarded to ``func``.

    Returns:
        float: the mixed second partial derivative with respect to ``dim1``
        and ``dim2``.
    """
    eps = kwargs.get("eps", 1e-5)
    priors = kwargs.get("priors")

    x_upup = x0.copy()
    x_updown = x0.copy()
    x_downup = x0.copy()
    x_downdown = x0.copy()

    # Scale the step by the parameter's prior width; this fails if it is zero.
    if priors is None:
        scale1 = 1.0
        scale2 = 1.0
    else:
        key = list(priors.keys())
        scale1 = _prior_scale(priors[key[dim1]])
        scale2 = _prior_scale(priors[key[dim2]])

    if scale1 == 0:
        scale1 = 1
    if scale2 == 0:
        scale2 = 1
    dim1adj = eps / 2 * scale1
    dim2adj = eps / 2 * scale2

    x_upup[dim1] += dim1adj
    x_upup[dim2] += dim2adj

    x_updown[dim1] += dim1adj
    x_updown[dim2] -= dim2adj

    x_downup[dim1] -= dim1adj
    x_downup[dim2] += dim2adj

    x_downdown[dim1] -= dim1adj
    x_downdown[dim2] -= dim2adj

    func_upup = func(x_upup, **kwargs)
    func_updown = func(x_updown, **kwargs)
    func_downup = func(x_downup, **kwargs)
    func_downdown = func(x_downdown, **kwargs)

    return (func_upup - func_updown - func_downup + func_downdown) / (4 * dim1adj * dim2adj)


def _prior_scale(prior):
    """Return the characteristic width of a single prior tuple."""
    kind = prior[0]
    if kind == "Gaussian":
        return prior[1]  # the mean
    if kind in ("Uniform", "Jeffreys"):
        return (prior[2] - prior[1]) / 2
    return 1.0


def calculate_hessian(func, vals, **kwargs):
    """Estimate the Hessian matrix of a generic function.

    Args:
        func (callable): ``N``-dimensional function to estimate the Hessian of.
        vals (array): length-``N`` point at which to estimate the matrix of
            second derivatives.
        **kwargs: forwarded to :func:`numerical_second_partial`.

    Returns:
        numpy.ndarray: the ``N x N`` Hessian of ``func`` at ``vals``.
    """
    hessian = np.zeros([len(vals), len(vals)])
    for i in range(len(vals)):
        for j in range(len(vals)):
            hessian[i][j] += numerical_second_partial(func, vals, i, j, **kwargs)

    return np.array(hessian)


def laplace_approximation(func, x0, **kwargs):
    r"""Calculate the Laplace approximation of an integral of a specific form.

    A challenging integral that can be written in terms of an exponent

    .. math::
        Z = \int \exp(f(x))\, dx

    can be estimated as approximately

    .. math::
        \left[\frac{(2\pi)^{2}}{\det|H(x_{0})|}\right]^{\frac{1}{2}} \exp(f(x_{0}))

    where ``H`` is the function's Hessian matrix and ``x0`` is a region of high
    probability.

    Args:
        func (callable): the function ``f(x)`` in the exponent of the term to
            estimate -- not the full integrand, just the exponent.
        x0 (array): length-``N`` local maximum around which to compute the
            approximation. In practice one should optimise ``func`` first and
            pass the optimum here.
        **kwargs: forwarded to ``func`` and :func:`calculate_hessian`.

    Returns:
        tuple: ``(logA, logB)`` where ``logA = f(x0)`` and ``logB`` is the
        logarithm of the Hessian term
        ``[(2*pi)^2 / det|H(x0)|]^(1/2)``. Logarithms are returned because these
        terms routinely overflow floating-point precision.
    """
    logA = func(x0, **kwargs)

    H = calculate_hessian(func, x0, **kwargs)
    logB = np.log(((np.pi * 2) ** 2 / (np.abs(np.linalg.det(H)))) ** (1 / 2))
    logger.debug("Determinant: %s", np.linalg.det(H))

    return logA, logB
