# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-30

First packaged release. No change to numerical behaviour — this pass makes the
project installable and maintainable.

### Added

- Installable `lapprox` package with a `src/` layout and a real PEP 621
  `pyproject.toml` (declared dependencies, optional `rv` / `docs` / `dev`
  extras).
- MIT `LICENSE`.
- Test suite (`tests/`) covering the Hessian estimator and the Laplace
  approximation against a known integral.
- GitHub Actions CI: linting, tests on Python 3.9–3.12, and a docs build.
- `.gitignore`; build artifacts and OS cruft removed from version control.

### Changed

- Public API renamed to PEP 8 (`snake_case`):
  - `Laplace_Approximation` → `laplace_approximation`
  - `Calculate_Hessian` → `calculate_hessian`
  - `NDeriv_2` → `numerical_second_partial`
  - `Calculate_Likelihood_RadVel` → `calculate_likelihood_radvel`
  - `LikelihoodXPrior` → `log_likelihood_times_prior`
  - `Prior_Components` → `prior_components`
- Modules renamed: `LApprox.py` → `lapprox/laplace.py`,
  `Likelihood_Functions.py` → `lapprox/likelihoods.py`.
- `radvel` is now an optional dependency, imported lazily.
- Stray `print` calls replaced with module loggers.
- Sphinx docs build cleanly with warnings treated as errors; switched to the
  Read the Docs theme.
- Example notebooks moved to `examples/` and updated to the new import paths.

### Removed

- Dead / commented-out code, duplicate imports, unused `scipy` imports.
- Committed `docs/build/` output, `.DS_Store` files, `__pycache__`.
