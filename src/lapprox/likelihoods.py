"""Radial-velocity likelihood helpers for use with :mod:`lapprox.laplace`.

These functions build a `radvel <https://radvel.readthedocs.io>`_ model from a
dictionary of parameter values and return its log-likelihood, optionally
multiplied by a set of priors. ``radvel`` is an optional dependency; install it
with ``pip install "lapprox[rv]"``.

.. note::
   The chromatic-kernel code paths (``kernel="KJ1"`` / ``"KJ2"``) rely on a
   custom ``radvel`` fork and will not run against the released package.
"""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def calculate_likelihood_radvel(vals, **kwargs):
    """Compute the RV log-likelihood of a model given parameter values and data.

    Args:
        vals (array): current value of each varying parameter, ordered to match
            ``priors``.

    Keyword Args:
        priors (dict): ``{name: (kind, a, b)}`` for the varying parameters.
        NPlanets (int): number of planets in the system (default ``1``).
        filename (str): path to the RV data (used if ``dataframe`` is absent).
        dataframe (pandas.DataFrame): RV data with ``time``, ``mnvel``, ``errvel``
            (and optionally ``tel``) columns.
        GP (bool): whether to use a Gaussian-process likelihood.
        hparam_dict (dict): GP hyper-parameters.
        kernel (str): GP kernel name.
        val_dict (dict): full parameter dictionary (varying and fixed).

    Returns:
        float: the composite log-probability of the model.
    """
    priors = kwargs.get("priors")
    if priors is None:
        logger.warning("No priors detected")
    NPlanets = kwargs.get("NPlanets", 1)
    filename = kwargs.get("filename")
    GP = kwargs.get("GP", False)
    hparam_dict = kwargs.get("hparam_dict")
    kernel = kwargs.get("kernel")
    val_dict = kwargs.get("val_dict")
    data = kwargs.get("dataframe")

    radvel = _import_radvel()

    # Update the value dictionary with the latest values of the varying params.
    for counter, key in enumerate(priors.keys()):
        val_dict[key] = vals[counter]

    if data is None:
        if filename is None:
            logger.warning("Pass a filename or dataframe to the likelihood")
        data = pd.read_csv(filename, sep=" ")

    t = np.array(data.time)
    vel = np.array(data.mnvel)
    errvel = np.array(data.errvel)

    time_base = np.min(t)  # reference time for the linear/curvature terms
    try:
        tel_arr = np.array(data.tel)
        telgrps = data.groupby("tel").groups
        instnames = telgrps.keys()
    except AttributeError:
        tel_arr = np.repeat("test", len(t))
        instnames = ["test"]

    telnames = np.array([inst for inst in instnames])

    vary_dict = {key: (key in priors.keys()) for key in val_dict.keys()}

    if hparam_dict is not None:
        hnames = np.array(list(hparam_dict.keys()))

    params = radvel.Parameters(NPlanets, basis="per tc e w k")

    for i in range(NPlanets):
        n = str(i + 1)
        params["per" + n] = radvel.Parameter(value=val_dict["per" + n], vary=vary_dict["per" + n])
        params["tc" + n] = radvel.Parameter(value=val_dict["tc" + n], vary=vary_dict["tc" + n])
        params["e" + n] = radvel.Parameter(value=val_dict["e" + n], vary=vary_dict["e" + n])
        params["w" + n] = radvel.Parameter(value=val_dict["w" + n], vary=vary_dict["w" + n])
        params["k" + n] = radvel.Parameter(value=val_dict["k" + n], vary=vary_dict["k" + n])

    # dvdt is the linear trend term and curv the curvature term.
    params["dvdt"] = radvel.Parameter(value=val_dict["dvdt"], vary=vary_dict["dvdt"])
    params["curv"] = radvel.Parameter(value=val_dict["curv"], vary=vary_dict["curv"])

    if GP:
        for key in hparam_dict.keys():
            params[key] = radvel.Parameter(value=val_dict[key], vary=vary_dict[key])

    for tel_suffix in instnames:
        params["gamma_" + tel_suffix] = radvel.Parameter(
            value=val_dict["gamma_" + tel_suffix], vary=vary_dict["gamma_" + tel_suffix]
        )
        params["jit_" + tel_suffix] = radvel.Parameter(
            value=val_dict["jit_" + tel_suffix], vary=vary_dict["jit_" + tel_suffix]
        )

    # Combine all parameters into a Keplerian model.
    model = radvel.model.RVModel(params, time_base=time_base)

    likes = []  # one likelihood per instrument

    def initialize(tel_suffix):
        # A separate likelihood object per instrument, sharing one RVModel.
        try:
            indices = telgrps[tel_suffix]
        except (KeyError, NameError):
            logger.warning("telgrps did not initialize properly")
            indices = np.ones(len(t), dtype=bool)
        if GP:
            like = radvel.likelihood.GPLikelihood(
                model, t[indices], vel[indices], errvel[indices],
                hnames[tel], suffix="_" + tel_suffix, kernel="QuasiPer",
            )
        else:
            like = radvel.likelihood.RVLikelihood(
                model, t[indices], vel[indices], errvel[indices], suffix="_" + tel_suffix,
            )
        like.params["gamma_" + tel_suffix] = radvel.Parameter(
            value=val_dict["gamma_" + tel_suffix], vary=True
        )
        like.params["jit_" + tel_suffix] = radvel.Parameter(
            value=val_dict["jit_" + tel_suffix], vary=True
        )
        likes.append(like)

    def initialize_chromatic(kernel_name):
        return radvel.likelihood.ChromaticLikelihood(
            model=model, t=t, vel=vel, errvel=errvel, suffix=telnames, hnames=hnames,
            kernel_name=kernel_name, tel=tel_arr, telnames=telnames,
        )

    if GP is False:
        for tel in instnames:
            initialize(tel)
    elif kernel == "KJ2":
        likes.append(initialize_chromatic("Chromatic_2"))
    elif kernel == "KJ1":
        likes.append(initialize_chromatic("Chromatic_1"))

    # Merge into a composite likelihood for the final calculation.
    like = radvel.likelihood.CompositeLikelihood(likes)

    return like.logprob()


def log_likelihood_times_prior(x0, **kwargs):
    """Return the RV log-likelihood plus the log-prior contribution.

    Args:
        x0 (array): current value of each varying parameter.

    Keyword Args:
        optimizing (bool): when ``True`` the sign is flipped so the result can be
            minimised. All keyword arguments are forwarded to
            :func:`calculate_likelihood_radvel`, including ``priors``.

    Returns:
        float: ``log L(x0) + sum(log prior_i(x0))``, negated if ``optimizing``.
    """
    loglikelihood = calculate_likelihood_radvel(x0, **kwargs)

    priors = kwargs.get("priors")
    if priors is None:
        logger.warning("You need a prior")
    optimizing = kwargs.get("optimizing", False)

    for counter, key in enumerate(priors.keys()):
        con, noncon = prior_components(priors[key], x0[counter])
        loglikelihood += np.log(con * noncon)

    return loglikelihood * -1 if optimizing else loglikelihood


def prior_components(prior, value):
    """Split a prior into its constant and non-constant parts.

    The non-constant part is multiplied by the likelihood before the Laplace
    approximation is applied; the constant part factors out of the integral and
    its logarithm is added to the evidence afterwards.

    Supported priors (``prior[0]`` selects the form):

    ===============  =========================================================
    ``Uniform``      ``U(a, b) = 1/(b-a)`` for ``a < x < b`` else ``0``
    ``Gaussian``     ``G(a, b) = 1/(b*sqrt(2*pi)) * exp(-0.5*(x-a)^2/b^2)``
    ``Jeffreys``     ``J(a, b) = 1/x * 1/log(b/a)``
    ===============  =========================================================

    Args:
        prior (tuple): ``(kind, a, b)``.
        value (float): the parameter value at which to evaluate the prior.

    Returns:
        tuple: ``(constant, nonconstant)``.
    """
    kind = prior[0]

    if kind == "Uniform":
        constant = 1 / (prior[2] - prior[1])
        nonconstant = 1.0 if (value > prior[1] and value < prior[2]) else 0.0
    elif kind == "Gaussian":
        constant = 1 / (prior[2] * np.sqrt(np.pi * 2))
        nonconstant = np.exp(-0.5 * (value - prior[1]) ** 2 / (prior[2] ** 2))
    elif kind == "Jeffreys":
        constant = 1 / np.log(prior[2] / prior[1])
        nonconstant = 1 / value
    else:
        logger.warning(
            "Unrecognised prior %r; add its functional form to prior_components", kind
        )
        constant = 1.0
        nonconstant = 1.0

    return constant, nonconstant


def _import_radvel():
    """Import ``radvel`` lazily so it stays an optional dependency."""
    try:
        import radvel
    except ImportError as exc:  # pragma: no cover - depends on optional install
        raise ImportError(
            "The RV likelihood helpers require radvel. Install it with "
            '`pip install "lapprox[rv]"`. The chromatic-kernel paths additionally '
            "need the custom radvel fork noted in examples/README.md."
        ) from exc
    return radvel
