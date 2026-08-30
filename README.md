# LApprox

[![Documentation Status](https://readthedocs.org/projects/lapprox/badge/?version=latest)](https://lapprox.readthedocs.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](pyproject.toml)

<img src="docs/source/_static/logo.png" alt="LApprox logo" width="180"/>

A lightweight Python package for computing the **Laplace approximation** and
using it to compare Bayesian models by their evidence.

## What it does

A challenging integral that can be written in exponential form,

```
Z = ∫ exp(f(x)) dx
```

can be estimated as

```
Z ≈ [ (2π)² / |det H(x₀)| ]^(1/2) · exp(f(x₀))
```

where `H` is the Hessian of `f` and `x₀` is a mode of the integrand. This is a
fast, cheap alternative to numerical integration when an exact answer is not
required — most usefully for the Bayesian evidence of a model. It is most
accurate when the integrand has a single dominant mode far from the integration
bounds.

The package provides a numerical Hessian estimator, the Laplace approximation
itself, and optional helpers for radial-velocity (RV) likelihoods built on
[`radvel`](https://radvel.readthedocs.io).

## Installation

```bash
pip install "git+https://github.com/CCBeard/LApprox.git"
```

Optional extras:

```bash
pip install "lapprox[rv]"    # radial-velocity likelihood helpers (needs radvel)
pip install -e ".[dev]"      # tests, linting and docs, for development
```

> **Note:** the chromatic-kernel RV example (`examples/rv_likelihood.ipynb`)
> depends on a custom `radvel` fork — see [`examples/README.md`](examples/README.md).

## Quickstart

```python
import numpy as np
from lapprox import laplace_approximation

# f(x, y) = -3x² - y²  →  ∫ exp(f) dx dy = π / √3 ≈ 1.8138
def f(x0, **kwargs):
    return -3 * x0[0] ** 2 - x0[1] ** 2

log_a, log_b = laplace_approximation(f, [0.0, 0.0])
print(np.exp(log_a + log_b))   # 1.8137993642342178
```

## Documentation

Full API reference: <https://lapprox.readthedocs.io>

Worked examples live in [`examples/`](examples/):

- `multivariate_gaussian.ipynb` — a self-contained example with a known answer.
- `rv_likelihood.ipynb` — an RV evidence calculation (requires the `radvel` fork).

## Citing

If you use LApprox in published work, please cite this repository. A Zenodo DOI
will be added here. <!-- TODO: mint DOI -->

## License

[MIT](LICENSE) © Corey Beard
