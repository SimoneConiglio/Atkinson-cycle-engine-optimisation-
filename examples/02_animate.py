"""Write the figures and the animated GIF into ``figures/``.

python examples/02_animate.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from exlink import analyse
from exlink.animation import animate_dashboard, save
from exlink.plots import (
    plot_cycle,
    plot_mechanism,
    plot_motion,
    plot_overview,
    plot_torque,
)
from exlink.reference import REFINED_DESIGN

OUTDIR = Path("figures")


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    design = REFINED_DESIGN
    analysis = analyse(design, samples=720)

    for name, figure in {
        "motion": plot_motion(design, analysis),
        "cycle": plot_cycle(design, analysis),
        "torque": plot_torque(design, analysis),
        "mechanism": plot_mechanism(design, analysis=analysis),
        "overview": plot_overview(design, analysis),
    }.items():
        path = OUTDIR / f"{name}.png"
        figure.savefig(path, dpi=140, bbox_inches="tight")
        print(f"written: {path}")

    animation = animate_dashboard(design, frames=120)
    print(f"written: {save(animation, OUTDIR / 'exlink.gif', fps=25, dpi=90)}")


if __name__ == "__main__":
    main()
