"""The two standard MDO pictures of the problem, drawn by GEMSEO.

Section 3 argues the formulation in prose and section 4 maps it onto modules.
Neither shows the *shape* of the problem, and two conventional diagrams do:

* an **N2 chart** -- the disciplines on the diagonal and their couplings off
  it.  Everything above the diagonal is feed-forward and everything below it is
  feedback, so a glance says whether an MDA is needed at all.  Here exactly one
  entry sits below the diagonal -- :class:`~exlink.disciplines.StructureDiscipline`
  returning the section diameters to
  :class:`~exlink.disciplines.DynamicsDiscipline` -- and that pair is the fixed
  point of section 3.3.
* an **XDSM** -- the same graph with the *process* on it: who calls whom, in
  what order, and how many times.  It is what distinguishes MDF from IDF at a
  glance, and it shows the MDA sitting inside the optimizer's loop, which is
  the whole content of the choice section 3.6 makes.

Both come from GEMSEO rather than being drawn by hand, so they cannot drift
from the scenario the study actually solves: they are generated from the very
object :func:`exlink.scenarios.build_range_scenario` hands to the optimizer.

What each output needs
----------------------
The N2 chart is matplotlib and writes a PNG anywhere.  The XDSM has two
back-ends: an interactive HTML one that works everywhere, and pyXDSM, which
emits TikZ and needs a LaTeX toolchain to compile.  The pyXDSM output is the
one worth publishing -- the HTML renderer overlaps its node captions with the
variable lists when the names are as long as this problem's -- so
:func:`write_xdsm` prefers it and falls back to HTML with a warning when
``pdflatex`` is absent.  Converting the PDF to a PNG for the documentation
needs ``pdftoppm`` from Poppler; without it the PDF is left in place.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Any

from .constants import DEFAULT_SPEC, DEFAULT_TARGETS, DesignTargets, EngineSpec
from .design import Design
from .reference import RELIABLE_DESIGN

if TYPE_CHECKING:  # pragma: no cover - typing only
    from gemseo.scenarios.base_scenario import BaseScenario

LOGGER = logging.getLogger(__name__)

PINNED_MODULE = 0.8
"""Gear module the published scenario is pinned to [mm]."""

PINNED_TEETH = 48
"""Teeth on the small gear, pinning ``I`` and leaving ten continuous variables."""

RENDER_DPI = 200
"""Resolution the XDSM PDF is rasterised at, for the documentation."""


def build_published_scenario(
    design: Design = RELIABLE_DESIGN,
    spec: EngineSpec = DEFAULT_SPEC,
    targets: DesignTargets = DEFAULT_TARGETS,
    **kwargs: Any,
) -> BaseScenario:
    """The scenario the diagrams are drawn from.

    The range problem of section 3.10 at the design section 6.4 arrives at,
    with the gear pair pinned as every run in the study pins it.  Kept in one
    place so the N2 chart and the XDSM cannot be drawn from different problems.

    Args:
        design: Starting design; the study's result by default.
        spec: Fixed engine data.
        targets: Constraint right-hand sides.
        **kwargs: Forwarded to
            :func:`exlink.scenarios.build_range_scenario`.

    Returns:
        The scenario.
    """
    from .scenarios import build_range_scenario

    return build_range_scenario(
        initial=design,
        module=PINNED_MODULE,
        teeth=PINNED_TEETH,
        spec=spec,
        targets=targets,
        **kwargs,
    )


def write_n2(
    path: str | Path = "n2.png",
    scenario: BaseScenario | None = None,
    fig_size: tuple[float, float] = (11.0, 8.0),
) -> Path:
    """Write the N2 chart of the coupled problem.

    Args:
        path: Destination; the suffix chooses the format.
        scenario: The scenario to read the disciplines from; the published one
            if omitted.
        fig_size: Figure size [in].

    Returns:
        The path written.
    """
    from gemseo import generate_n2_plot

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    chosen = scenario if scenario is not None else build_published_scenario()
    generate_n2_plot(chosen.disciplines, target, fig_size=fig_size, save=True, show=False)
    # GEMSEO drops an interactive companion beside the image whether or not one
    # was asked for; the documentation ships the static chart only.
    target.with_suffix(".html").unlink(missing_ok=True)
    return target


def write_xdsm(
    directory: str | Path = ".",
    scenario: BaseScenario | None = None,
    file_name: str = "xdsm",
    dpi: int = RENDER_DPI,
) -> Path:
    """Write the XDSM of the coupled problem, as a PNG where the tools allow.

    Prefers the pyXDSM back-end, which needs ``pdflatex``, and rasterises its
    PDF with ``pdftoppm``.  Falls back a step at a time: to the PDF if Poppler
    is missing, and to the interactive HTML if LaTeX is missing.  The return
    value says which was produced.

    Args:
        directory: Directory to write into.
        scenario: The scenario; the published one if omitted.
        file_name: Stem of the written file.
        dpi: Resolution for the PNG.

    Returns:
        The path written -- ``.png``, ``.pdf`` or ``.html``.
    """
    from gemseo import generate_xdsm

    target = Path(directory)
    target.mkdir(parents=True, exist_ok=True)
    chosen = scenario if scenario is not None else build_published_scenario()

    if shutil.which("pdflatex") is None:
        LOGGER.warning(
            "pdflatex not found: writing the interactive HTML XDSM instead. "
            "Its node captions overlap the variable lists on this problem."
        )
        generate_xdsm(
            chosen,
            directory_path=target,
            file_name=file_name,
            save_html=True,
            show_html=False,
        )
        return target / f"{file_name}.html"

    generate_xdsm(
        chosen,
        directory_path=target,
        file_name=file_name,
        save_html=False,
        save_pdf=True,
        pdf_build=True,
        pdf_cleanup=True,
        show_html=False,
    )
    pdf = target / f"{file_name}.pdf"
    if shutil.which("pdftoppm") is None:
        LOGGER.warning("pdftoppm not found: leaving the XDSM as a PDF.")
        return pdf

    subprocess.run(
        [
            "pdftoppm",
            "-png",
            "-r",
            str(int(dpi)),
            "-singlefile",
            str(pdf),
            str(target / file_name),
        ],
        check=True,
    )
    pdf.unlink(missing_ok=True)
    # pyXDSM leaves its TikZ source behind even under ``pdf_cleanup``.
    for suffix in (".tex", ".tikz"):
        (target / f"{file_name}{suffix}").unlink(missing_ok=True)
    png = target / f"{file_name}.png"
    _trim(png)
    return png


def _trim(path: Path, pad: int = 16) -> None:
    """Crop a rasterised diagram to its ink, leaving a small margin.

    ``pdftoppm`` pads to the PDF's page box, which pyXDSM sizes generously; the
    result is a figure that reads small in a document because most of it is
    white.  Args: ``path`` the image, ``pad`` the margin to leave [px].
    """
    from PIL import Image

    image = Image.open(path).convert("RGB")
    mask = image.convert("L").point(lambda value: 255 if value < 250 else 0)
    box = mask.getbbox()
    if box is None:
        return
    left, top, right, bottom = box
    image.crop(
        (
            max(left - pad, 0),
            max(top - pad, 0),
            min(right + pad, image.width),
            min(bottom + pad, image.height),
        )
    ).save(path, optimize=True)
