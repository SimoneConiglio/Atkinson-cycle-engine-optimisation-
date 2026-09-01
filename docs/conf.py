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
    "sphinxcontrib.bibtex",
]

# One bibliography file for the whole site.  ``author_year`` is the citation
# style expected in the engineering-optimization literature this work sits in:
# the text names the authors, so footnote-style numbering would hide exactly
# the information the sentence depends on.
#
# sphinxcontrib-bibtex ships no bibliography style whose *labels* are author
# and year, so an inline "Duran and Grossmann [1986]" would point at an entry
# labelled "[4]".  The style below supplies matching labels; it is the recipe
# from the sphinxcontrib-bibtex documentation, narrowed to what is needed here.
bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"


def _setup_author_year_labels() -> str:
    from pybtex.plugin import register_plugin
    from pybtex.style.formatting.unsrt import Style as UnsrtStyle
    from pybtex.style.labels import BaseLabelStyle

    class AuthorYearLabelStyle(BaseLabelStyle):
        def format_labels(self, sorted_entries):
            for entry in sorted_entries:
                persons = entry.persons.get("author") or entry.persons.get("editor") or []
                names = [
                    (" ".join(str(part) for part in person.last_names)
                     if person.last_names else str(person))
                    .replace("{", "")
                    .replace("}", "")
                    for person in persons
                ]
                if not names:
                    who = entry.fields.get("organization", entry.key)
                elif len(names) == 1:
                    who = names[0]
                elif len(names) == 2:
                    who = f"{names[0]} and {names[1]}"
                else:
                    who = f"{names[0]} et al."
                yield f"{who}, {entry.fields.get('year', 'n.d.')}"

    class AuthorYearStyle(UnsrtStyle):
        default_sorting_style = "author_year_title"
        default_label_style = AuthorYearLabelStyle

    register_plugin("pybtex.style.formatting", "authoryear", AuthorYearStyle)
    return "authoryear"


bibtex_default_style = _setup_author_year_labels()

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

# MyST parses no math by default, so `$...$` and `$$...$$` would render as
# literal dollar signs.  ``dollarmath`` turns both into MathJax; ``amsmath``
# additionally accepts bare LaTeX environments such as ``\begin{align}``.
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "colon_fence",
    "deflist",
    "smartquotes",
]
# Allow a display equation to follow text on the same line, and permit digits
# immediately after a closing ``$`` so that "$10^{-3}$ target" parses.
myst_dmath_double_inline = True
myst_dmath_allow_digits = True
myst_heading_anchors = 3

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
