# LApprox
Lightweight package for utilizing the Laplace Approximation to compare Bayesian models.

![LA](LA.png)


Project documentation can be found here: [LApprox Documentation](https://lapprox.readthedocs.io/en/latest/index.html)

LApprox.py is the only file required to run this package, though we include a Likelihood_Functions.py file for convenience. This project was designed for the purpose of calculating the Laplace Approximation for RV likelihood functions, and these are non-trivial to work with. Hence, we share our implementation.

An example of a simple, exact computation done for a multivariate Gaussian distribution is available in the Example_Multivariate_Gaussian.ipynb file.

A more complicated example, utilizing our Radial Velocity (RV) likelihood functions, is visible in Example_RV_Likelihood.ipynb.
