"""The documentation is built from these docstrings, so they are checked too.

Two failures are possible here that no other test would catch, because both
produce a *valid* string that renders wrongly rather than an error.

The first is a LaTeX command in a docstring that Python reads as an escape.
``:math:`\\beta``` in a non-raw string is a backspace followed by ``eta``, and
``\\times`` is a tab followed by ``imes`` -- the backslash is consumed at
compile time, so Sphinx never sees it and the rendered page shows a control
character where a symbol should be.  This is how ``\\beta = 3.00`` reached the
published API page as ``<BS>eta = 3.00``.

The second is a ``\\tag`` in a display equation that also contains ``\\\\``.
Sphinx's MathJax writer wraps any such equation in ``\\begin{split}``, where
``\\tag`` is illegal, and MathJax refuses the whole equation.  The fix is the
``math`` directive with ``:nowrap:`` and an explicit ``equation`` environment;
this test is what stops a plain ``$$...$$`` from creeping back.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

SOURCE_ROOTS = ("src", "tests", "examples")
DOCS = Path(__file__).resolve().parents[1] / "docs"

# The escapes Python interprets that also begin a LaTeX command someone would
# plausibly write: \a \b \f \n \r \t \v.  Their compiled forms are control
# characters, and a newline is the only one a docstring has any reason to hold.
# Spelled as ordinals rather than as escapes so that this module does not trip
# its own check.
INTERPRETED = frozenset(chr(c) for c in (0x07, 0x08, 0x0C, 0x0D, 0x09, 0x0B))

DISPLAY_MATH = re.compile(r"\$\$(.*?)\$\$", re.S)


def _python_files() -> list[Path]:
    root = Path(__file__).resolve().parents[1]
    return sorted(p for d in SOURCE_ROOTS for p in (root / d).rglob("*.py"))


@pytest.mark.parametrize("path", _python_files(), ids=lambda p: p.name)
def test_no_latex_command_is_eaten_by_a_python_escape(path: Path) -> None:
    """No string literal holds a control character a backslash should have kept."""
    offenders = []
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            found = sorted(INTERPRETED.intersection(node.value))
            if found:
                offenders.append((node.lineno, [hex(ord(c)) for c in found]))
    assert not offenders, (
        f"{path.name} has string literals holding control characters: {offenders}. "
        "A LaTeX command was written in a non-raw string -- mark the literal r'' "
        "or double the backslash."
    )


@pytest.mark.parametrize("path", sorted(DOCS.glob("*.md")), ids=lambda p: p.name)
def test_no_tagged_equation_would_be_wrapped_in_split(path: Path) -> None:
    """A tagged multi-line equation must use the nowrap directive, not ``$$``."""
    offenders = [
        re.search(r"\\tag\{([^}]*)\}", m.group(1)).group(1)  # type: ignore[union-attr]
        for m in DISPLAY_MATH.finditer(path.read_text())
        if "\\tag" in m.group(1) and ("\\\\" in m.group(1) or "\n\n" in m.group(1))
    ]
    assert not offenders, (
        f"{path.name} tags equations that Sphinx will wrap in \\begin{{split}}, "
        f"where \\tag is illegal: {offenders}. Use the math directive with "
        ":nowrap: and an explicit equation environment."
    )
