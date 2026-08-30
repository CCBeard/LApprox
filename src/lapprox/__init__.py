"""LApprox: the Laplace approximation for Bayesian model comparison.

The public API lives here:

>>> from lapprox import laplace_approximation, calculate_hessian
"""

from lapprox.laplace import (
    calculate_hessian,
    laplace_approximation,
    numerical_second_partial,
)

__version__ = "0.1.0"

__all__ = [
    "laplace_approximation",
    "calculate_hessian",
    "numerical_second_partial",
    "__version__",
]
