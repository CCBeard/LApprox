.. LApprox documentation master file, created by
   sphinx-quickstart on Fri Sep 22 09:35:15 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to LApprox's documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   LApprox.rst
   Likelihood_Functions.rst

LApprox is a lightweight python module intended to help when calculating the Laplace Approximation on a variety of models.

What is the Laplace approximation?

The Laplace approximation is a fast, computationally inexpensive way to calculate the value of an integral of a particular form.

A challenging integral, when possible to write in terms of an exponential:
    .. math::
        Z = \int(exp(f(x))dx)

Can be estimated as approximately:

    .. math::
        [\frac{(2\pi)^{2}}{det|H(x_{0})|}]^{\frac{1}{2}} * exp(f(x_{0}))

where H is the function's Hessian matrix, and x0 is a region of high probability.


Many integrals of interest in a wide variety of scientific fields are impossible to calculate analytically, and must
be approximated numerically. Often, even calculating them numerically is computationally intractable, especially if the
operation must be performed many times. The Laplace Approximation can be a convenient workaround if an exact answer is not necessary.

When is the approximation accurate?

The Laplace approximation is most accurate when the function being integrated has a single dominant mode x0, and when that mode is far from the bounds of integration.

A common situation where the Laplace approximation might be useful is the calculation of the Bayesian evidence of a model, as this is often a very complicated integral.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
