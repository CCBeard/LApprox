# Examples

Run these from this directory after installing the package (`pip install -e ".[rv]"`
from the repo root).

## `multivariate_gaussian.ipynb`

A fully self-contained walkthrough: define a 2-D Gaussian exponent, estimate its
Hessian, and check the Laplace approximation against the exact integral
`π / √3 ≈ 1.8138`. Needs only `numpy`, `scipy` and `matplotlib`.

## `rv_likelihood.ipynb`

Computes the Bayesian evidence of a radial-velocity model for Kepler-21 using the
helpers in `lapprox.likelihoods` and the data in `data/Kepler21_rv.csv`.

> **This notebook is not reproducible with the released `radvel`.**
> It runs with `GP=True` and `Kernel="KJ1"`, which call
> `radvel.likelihood.ChromaticLikelihood` and the `"Chromatic_1"` kernel — both
> from a custom `radvel` fork, not mainline `radvel`. To run it you need that
> fork installed. For a version that works against released `radvel`, set
> `gp = False` and remove the GP hyper-parameters and priors.
