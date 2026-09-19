"""The design box cut into boxes, and what that is worth here.

Runs the geometric efficiency problem three ways over the same 2 x 2 x 3 grid,
which is the comparison §5.11 reports:

1. the box subdivision on its **default** box start, the center of the box,
   which finds nothing at all -- the center is on the penalty plateau, where the
   objective and five of the seven constraints are flat;
2. the box subdivision on the **restored** box start of
   :mod:`exlink.subdivision`, which finds designs;
3. **every box enumerated** with the same start policy and the same sub-problem
   solver, which is what the master has to be measured against: it is the only
   way to tell a master that stopped early from a subdivision that had nothing
   left to find.

Takes about three minutes.

    python examples/13_box_subdivision.py
"""

from __future__ import annotations

import time

from gemseo_box_subdivision import BoxSubdivisionSettings

from exlink.cli import configure_logging
from exlink.design import GLOBAL_BOUNDS
from exlink.reference import REFINED_DESIGN
from exlink.scenarios import _best_design, analyse, build_scenario, is_feasible
from exlink.subdivision import (
    enumerate_boxes,
    format_enumeration,
    format_subdivision,
    maximise_efficiency_by_subdivision,
)

GRID = {"q_1": 2, "q_2": 2, "theta_f": 3}
SAMPLES = 180


def baseline() -> tuple[float, int]:
    """What §4.2 already reports: one SLSQP solve from the refined reference."""
    scenario = build_scenario(
        "neg_efficiency",
        bounds=GLOBAL_BOUNDS,
        initial=REFINED_DESIGN,
        samples=SAMPLES,
        relax_equalities=True,
    )
    scenario.execute(algo_name="SLSQP", max_iter=60)
    design = _best_design(scenario, objective="neg_efficiency", fallback=REFINED_DESIGN)
    analysis = analyse(design, samples=SAMPLES)
    evaluations = len(scenario.formulation.optimization_problem.database)
    if not is_feasible(analysis):
        return float("nan"), evaluations
    return float(analysis.metrics.efficiency), evaluations


def main() -> None:
    configure_logging(verbose=False)

    efficiency, evaluations = baseline()
    print(f"single SLSQP solve from the reference: eta = {efficiency:.4f} "
          f"in {evaluations} analyses\n")

    settings = BoxSubdivisionSettings(
        convexity_margin=0.13,  # the spread of the objective over the boxes
        trust_region_radius=2,
        max_iter=40,
        sub_problem_max_iter=60,
    )

    cold = maximise_efficiency_by_subdivision(
        GRID, samples=SAMPLES, restoring=False, settings=settings
    )
    print(format_subdivision(cold, "box subdivision, box-centre starts"))
    print()

    warm = maximise_efficiency_by_subdivision(
        GRID, samples=SAMPLES, seed=0, probes=256, restarts=2, settings=settings
    )
    print(format_subdivision(warm, "box subdivision, restored starts"))
    print()

    began = time.perf_counter()
    boxes = enumerate_boxes(GRID, samples=SAMPLES, seed=0, probes=256, restarts=2)
    print(format_enumeration(boxes, "every box, restored starts"))
    print(f"  {time.perf_counter() - began:.0f} seconds")


if __name__ == "__main__":
    main()
