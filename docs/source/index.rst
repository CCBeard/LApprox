Welcome to LApprox's documentation!
===================================

LApprox is a lightweight Python package for computing the Laplace approximation
on a variety of models, and for comparing Bayesian models via their evidence.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api

What is the Laplace approximation?
----------------------------------

The Laplace approximation is a fast, computationally inexpensive way to estimate
the value of an integral of a particular form.

A challenging integral, when it can be written in terms of an exponential

.. math::
   Z = \int \exp(f(x))\, dx

can be estimated as approximately

.. math::
   \left[\frac{(2\pi)^{2}}{\det|H(x_{0})|}\right]^{\frac{1}{2}} \exp(f(x_{0}))

where :math:`H` is the function's Hessian matrix and :math:`x_{0}` is a region of
high probability.

Many integrals of interest across the sciences cannot be evaluated analytically
and must be approximated numerically. Often even the numerical evaluation is
intractable, especially when it must be repeated many times. The Laplace
approximation is a convenient workaround when an exact answer is not required --
for example when calculating the Bayesian evidence of a model.

When is the approximation accurate?
-----------------------------------

The Laplace approximation is most accurate when the integrand has a single
dominant mode :math:`x_{0}` and that mode is far from the bounds of integration.

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
