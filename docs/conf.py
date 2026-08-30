"""Sphinx configuration.

Build with ``tox -e docs`` or::

    sphinx-build -b html docs docs/_build/html
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

project = "exlink"
author = "Simone Coniglio"
copyright = "2015-2026, Simone Coniglio"  # noqa: A001
release = "1.0.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx_autodoc_typehints",
    "myst_parser",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
}

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
}
napoleon_google_docstring = True
napoleon_numpy_docstring = False

templates_path: list[str] = []
exclude_patterns = ["_build", "figures"]
html_theme = "furo"
html_static_path: list[str] = []
