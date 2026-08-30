"""Tests for the core Laplace-approximation routines.

The two-dimensional Gaussian ``f(x, y) = -3x^2 - y^2`` has a known integral,

    integral of exp(f) over R^2 = sqrt(pi/3) * sqrt(pi) = pi / sqrt(3),

and the Laplace approximation is exact for a Gaussian, so these are tight checks.
"""

import numpy as np
import pytest

from lapprox import calculate_hessian, laplace_approximation

EXACT_INTEGRAL = np.pi / np.sqrt(3)


def gaussian_exponent(x0, **kwargs):
    return -3 * x0[0] ** 2 - x0[1] ** 2


def test_laplace_approximation_matches_exact_integral():
    log_a, log_b = laplace_approximation(gaussian_exponent, [0.0, 0.0])
    assert np.exp(log_a + log_b) == pytest.approx(EXACT_INTEGRAL, rel=1e-3)


def test_calculate_hessian_of_quadratic():
    hessian = calculate_hessian(gaussian_exponent, [0.0, 0.0])
    expected = np.array([[-6.0, 0.0], [0.0, -2.0]])
    np.testing.assert_allclose(hessian, expected, atol=1e-4)


@pytest.mark.parametrize("n", [1, 2, 3, 5])
def test_hessian_is_square_with_input_dimension(n):
    def paraboloid(x0, **kwargs):
        return -np.sum(np.asarray(x0) ** 2)

    hessian = calculate_hessian(paraboloid, np.zeros(n))
    assert hessian.shape == (n, n)
