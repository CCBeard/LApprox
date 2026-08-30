# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Prefer the installed package; fall back to the source tree for local builds.
sys.path.insert(0, os.path.abspath("../../src"))

try:
    from lapprox import __version__ as _version
except ImportError:
    _version = "0.1.0"

# -- Project information ------------------------------------------------------

project = "LApprox"
copyright = "2023-2026, Corey Beard"
author = "Corey Beard"
release = _version
version = _version

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

templates_path = ["_templates"]
exclude_patterns = []

autodoc_mock_imports = ["radvel", "pandas"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_logo = "_static/logo.png"
html_theme_options = {
    "logo_only": True,
}
