"""The gear choice as a real mixed-integer problem: bi-level outer approximation.

The discrete part of this design -- which standard module, how many teeth --
was previously handled by solving the two or three nearest buildable lattice
points and keeping the best.  That is enumeration with an arbitrary budget: no
bound, no stopping criterion, and no way to know whether an unvisited point
would have won.

This module states it as the mixed-integer nonlinear program it is and hands it
to the ``gemseo-bilevel-outer-approximation`` plugin, which implements the
Duran-Grossmann decomposition directly.

The formulation
---------------
.. code-block:: text

    main       gear_choice, a one-hot selection over the lattice   (categorical)
               solved by BILEVEL_MASTER_OUTER_APPROXIMATION
                    |                                    ^
                    | I, gear_module, gear_teeth         | linearisations and
                    v      (catalogue interpolation)     | feasibility
    sub        the ten remaining linkage variables       (continuous)
               MDF over the coupled disciplines, solved by SLSQP

The ``Benders`` formulation performs the split itself.  Categorical variables
go to the main problem; everything continuous goes to a sub-scenario which it
wraps in an ``MDOScenarioAdapterBenders``, so the main problem optimises the
sub-problem's *optimum*.  The sub-problem keeps the MDF formulation and the
exact Jacobians this package already supplies, and SLSQP solves it.

The lattice as a catalogue
---------------------------
The gear pair reaches the physics through three numbers: the centre distance
``I = 1.5 m z``, which is dominant and geometric, and the module and tooth count
themselves, which set the blank mass and the face width.  A single
``CatalogueDesignSpace`` categorical variable drives all three through one
interpolation chain, so the main problem picks a lattice point and the three
follow.  At unit SIMP penalty the interpolation is exactly ``I = sum_j y_j I_j``:
linear in the selection, analytically differentiable, and therefore something
the outer-approximation master can linearise.

This is why :class:`~exlink.disciplines.RangeDiscipline` takes the gear pair as
*inputs* rather than construction data -- under this formulation they are set by
the main problem.

Infeasible sub-problems
-----------------------
At several lattice points the continuous problem has no feasible solution:
pinning ``I`` throws the design off the equalities ``STE = 74`` and
``epsilon = 16``, which ``I`` was one of the variables used to satisfy.

Outer approximation handles that by construction, and it is the reason the
constraints here are attached with ``main_level=True``.  That flag tells the
formulation the constraint depends on the categorical choice and can render the
sub-problem infeasible, so the main problem carries an ``is_feasible``
condition alongside the linearisations.  A lattice point whose sub-problem
cannot be solved is then excluded on evidence rather than silently returning a
number, which is exactly the information the previous enumeration discarded.

The honest caveat
-----------------
Finite convergence of outer approximation to the global optimum requires the
sub-problem convex in the continuous variables for each fixed categorical
choice.  This problem is emphatically nonconvex, so the main problem's bound is
a bound under an assumption this problem violates.  Because the lattice is
small, :func:`exhaustive` solves every candidate, so what the decomposition
returns can be checked against the true best over the lattice.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from gemseo import create_scenario
from gemseo.settings.formulations import MDF_Settings
from gemseo.settings.opt import SLSQP_Settings
from gemseo_bilevel_outer_approximation.algos.design_space.catalogue_design_space import (
    CatalogueDesignSpace,
)

from .constants import DEFAULT_SPEC, DEFAULT_TARGETS, DesignTargets, EngineSpec
from .design import VARIABLE_NAMES, Bounds, Design
from .dynamics import DEFAULT_SPEED_RPM
from .gears import lattice_inter_axle
from .materials import DEFAULT_MATERIAL, DEFAULT_SAFETY, Material, SafetyFactors

RANGE_VIOLATION_OUTPUTS: tuple[str, ...] = ("runs_violation", "gear_violation")
"""The range constraints in violation form, which the adapter can address."""

GEAR_CHOICE = "gear_choice"
"""The categorical variable: which lattice point, as a one-hot selection."""

MASTER_ALGORITHM = "BILEVEL_MASTER_OUTER_APPROXIMATION"
"""The main-problem solver supplied by the plugin."""


@dataclass(frozen=True)
class LatticePoint:
    """One buildable gear pair, and the centre distance it forces."""

    module: float
    teeth: int

    @property
    def inter_axle(self) -> float:
        """``I = 1.5 m z`` [mm]."""
        return lattice_inter_axle(self.module, self.teeth)

    def __str__(self) -> str:  # pragma: no cover - presentation only
        return f"m={self.module:g}, z={self.teeth} (I={self.inter_axle:.3f})"


@dataclass
class MinlpResult:
    """Outcome of the mixed-integer solve."""

    design: Design | None
    point: LatticePoint | None
    km_per_litre: float
    feasible: bool
    iterations: int
    """Main-problem iterations, i.e. sub-problems solved."""

    seconds: float
    candidates: list[LatticePoint] = field(default_factory=list)
    scenario: Any = None
    reason: str = ""

    @property
    def solves(self) -> int:
        """Sub-problems solved -- the cost that matters."""
        return self.iterations


def candidates_from_design(
    design: Design,
    speed_rpm: float = DEFAULT_SPEED_RPM,
    count: int = 3,
    samples: int = 360,
    spec: EngineSpec = DEFAULT_SPEC,
    limit: int = 6,
) -> list[LatticePoint]:
    """Buildable lattice points around a design, nearest workable first.

    Ranked by whether the pair can carry the tooth load and only then by
    distance -- see :func:`exlink.gears.buildable_neighbours` for why the
    obvious ranking is the wrong one.

    Args:
        design: The design whose centre distance to bracket.
        speed_rpm: Crankshaft speed [rev/min], for the tooth load.
        count: Tooth counts to consider on each side, per module.
        samples: Crank angles per revolution.
        spec: Fixed engine data.
        limit: Most candidates to return.

    Returns:
        The candidate lattice points.
    """
    from .coupled import solve_for_design
    from .gears import buildable_neighbours

    sized = solve_for_design(design, speed_rpm=speed_rpm, samples=samples, spec=spec)
    tangential = float(np.max(np.abs(sized.loads.gear_force))) * math.cos(spec.pressure_angle)
    ranked = buildable_neighbours(design.I, tangential, count=count)
    return [LatticePoint(module, teeth) for module, teeth, _value, _pair in ranked][:limit]


def build_minlp_scenario(
    candidates: list[LatticePoint],
    initial: Design,
    bounds: Bounds | None = None,
    speed_rpm: float = DEFAULT_SPEED_RPM,
    samples: int | None = None,
    sub_max_iter: int = 30,
    spec: EngineSpec = DEFAULT_SPEC,
    targets: DesignTargets = DEFAULT_TARGETS,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
) -> Any:
    """Assemble the bi-level mixed-integer scenario.

    Args:
        candidates: The lattice points the main problem may choose between.
        initial: Starting design; its ``I`` is supplied by the catalogue instead.
        bounds: The continuous design box.
        speed_rpm: Crankshaft speed [rev/min].
        samples: Crank angles per revolution.
        sub_max_iter: SLSQP budget for each sub-problem.
        spec: Fixed engine data.
        targets: Constraint right-hand sides.
        material: The material.
        safety: The design factors.

    Returns:
        A scenario to execute with :data:`MASTER_ALGORITHM`.

    Raises:
        ValueError: If no candidates are given.
    """
    from gemseo.core.mdo_functions.mdo_function import MDOFunction

    from .disciplines import (
        COUPLED_SAMPLES,
        DynamicsDiscipline,
        ExlinkDiscipline,
        RangeDiscipline,
        StructureDiscipline,
    )
    from .scenarios import (
        COUPLED_INEQUALITY_OUTPUTS,
        DEFAULT_MDA,
        DEFAULT_MDA_SETTINGS,
        EQUALITY_OUTPUTS,
        INEQUALITY_OUTPUTS,
        BearingMarginDiscipline,
    )

    if not candidates:
        msg = "the main problem needs at least one lattice point to choose between"
        raise ValueError(msg)

    n_samples = COUPLED_SAMPLES if samples is None else samples
    box = Bounds.around(initial, relative=0.35) if bounds is None else bounds

    # -- design space: continuous linkage variables, plus the categorical choice --
    space = CatalogueDesignSpace()
    start = np.clip(initial.to_array(), box.lower, box.upper)
    for index, name in enumerate(VARIABLE_NAMES):
        if name == "I":
            continue  # supplied by the catalogue interpolation
        space.add_variable(
            name,
            lower_bound=float(box.lower[index]),
            upper_bound=float(box.upper[index]),
            value=float(start[index]),
            type_=CatalogueDesignSpace.DesignVariableType.FLOAT,
        )
    nearest = int(np.argmin([abs(item.inter_axle - initial.I) for item in candidates]))
    space.add_categorical_variable(
        name=GEAR_CHOICE,
        value=[nearest],
        catalogue=list(range(len(candidates))),
    )

    # -- one interpolation chain, three catalogue outputs -------------------------
    gears = space.get_catalogue_interpolation_discipline_from_dict(
        variable=GEAR_CHOICE,
        dictionary={
            "I": {
                "catalogue": np.array([item.inter_axle for item in candidates]),
                "penalty": 1.0,
            },
            "gear_module": {
                "catalogue": np.array([item.module for item in candidates]),
                "penalty": 1.0,
            },
            "gear_teeth": {
                "catalogue": np.array([float(item.teeth) for item in candidates]),
                "penalty": 1.0,
            },
        },
    )

    disciplines = [
        gears,
        ExlinkDiscipline(samples=n_samples, spec=spec, targets=targets),
        DynamicsDiscipline(
            speed_rpm=speed_rpm,
            samples=n_samples,
            material=material,
            safety=safety,
            spec=spec,
        ),
        StructureDiscipline(samples=n_samples, material=material, safety=safety, spec=spec),
        BearingMarginDiscipline(limit=targets.max_bearing_load),
        RangeDiscipline(
            speed_rpm=speed_rpm,
            samples=n_samples,
            material=material,
            safety=safety,
            spec=spec,
        ),
    ]

    scenario = create_scenario(
        disciplines=disciplines,
        formulation_name="Benders",
        objective_name="neg_range",
        design_space=space,
        maximize_objective=False,
        sub_problem_formulation_settings=MDF_Settings(
            main_mda_name=DEFAULT_MDA, main_mda_settings=dict(DEFAULT_MDA_SETTINGS)
        ),
        sub_problem_algo_settings=SLSQP_Settings(max_iter=sub_max_iter),
    )

    # Every constraint here depends on the categorical choice through ``I``, and
    # pinning ``I`` to a lattice point is precisely what can make the sub-problem
    # infeasible.  ``main_level=True`` is how that is declared: the sub-problem
    # still enforces the constraint, and the main problem additionally carries
    # the feasibility condition, so an unsolvable lattice point is excluded on
    # evidence instead of quietly returning a number.
    inequality = MDOFunction.ConstraintType.INEQ
    for name in INEQUALITY_OUTPUTS:
        scenario.add_constraint(name, constraint_type=inequality, main_level=True)
    # The two equalities go in as equalities.  The plugin's source notes that
    # the master does not linearise equality constraints and suggests relaxing
    # them into inequality bands first -- but doing that here is what makes the
    # run exhaust memory: four extra sub-problem constraints, each carrying a
    # post-optimal sensitivity through an MDF sub-problem whose coupling
    # variables are the 45 367 load-history entries
    # :func:`exlink.formulations.coupling_dimension` counts.  With equalities
    # the formulation adds a single ``is_feasible`` condition instead and the
    # run fits in 0.5 GB.  The cost is that these two constraints inform the
    # master only through feasibility, not through a linearisation.
    for name in EQUALITY_OUTPUTS:
        scenario.add_constraint(
            name, constraint_type=MDOFunction.ConstraintType.EQ, main_level=True
        )
    for name in COUPLED_INEQUALITY_OUTPUTS:
        scenario.add_constraint(name, constraint_type=inequality, main_level=True)
    # The violation form, not the margin form: see RANGE_OUTPUTS for why
    # ``positive=True`` cannot be used under this formulation.
    for name in RANGE_VIOLATION_OUTPUTS:
        scenario.add_constraint(name, constraint_type=inequality, main_level=True)
    return scenario


def solve(
    initial: Design,
    candidates: list[LatticePoint] | None = None,
    bounds: Bounds | None = None,
    speed_rpm: float = DEFAULT_SPEED_RPM,
    max_iter: int = 20,
    sub_max_iter: int = 25,
    samples: int | None = None,
    posa: float = 2.0,
    adapt: bool = True,
    min_dfk: float = 0.0,
    max_step: float = 1000.0,
    **kwargs: Any,
) -> MinlpResult:
    """Solve the gear choice and the linkage together, by bi-level outer approximation.

    Args:
        initial: Starting design.
        candidates: Lattice points to choose between; derived from ``initial``
            by :func:`candidates_from_design` if omitted.
        bounds: The continuous design box.
        speed_rpm: Crankshaft speed [rev/min].
        max_iter: Main-problem iteration budget.
        sub_max_iter: SLSQP budget per sub-problem.
        samples: Crank angles per revolution.
        posa: Post-optimal sensitivity amplification.  Multiplies the cut
            slopes, making each linearisation steeper and so less likely to cut
            off a lattice point that is actually better.  ``1.0`` is the raw
            linearisation, which on a nonconvex problem terminates early: with
            it this problem stops after two sub-solves, 0.6 % short of the best
            lattice point.
        adapt: Adaptive convexification.  Corrects the slopes by a secant
            method over the visited history so the cuts behave like valid
            supports even where the true value function is not convex.  This is
            the plugin's own answer to premature convergence.
        min_dfk: Convexity margin used by the adaptive correction.
        max_step: Main-problem trust radius.
        **kwargs: Forwarded to :func:`build_minlp_scenario`.

    Returns:
        The result, with the chosen lattice point and the design that goes
        with it.
    """
    from .performance import evaluate

    began = time.perf_counter()
    points = (
        candidates_from_design(initial, speed_rpm=speed_rpm)
        if candidates is None
        else candidates
    )
    scenario = build_minlp_scenario(
        points,
        initial,
        bounds=bounds,
        speed_rpm=speed_rpm,
        samples=samples,
        sub_max_iter=sub_max_iter,
        **kwargs,
    )
    reason = ""
    try:
        scenario.execute(
            algo_name=MASTER_ALGORITHM,
            max_iter=max_iter,
            normalize_design_space=False,
            posa=posa,
            adapt=adapt,
            min_dfk=min_dfk,
            max_step=max_step,
        )
    except Exception as error:  # a failed run is a result, not a crash
        reason = f"{type(error).__name__}: {error}"

    design, point = _extract(scenario, points, initial)
    outcome = (
        evaluate(
            design,
            speed_rpm=speed_rpm,
            module=point.module,
            teeth=point.teeth,
            samples=samples or 360,
        )
        if design is not None and point is not None
        else None
    )
    problem = scenario.formulation.optimization_problem
    return MinlpResult(
        design=design,
        point=point,
        km_per_litre=outcome.km_per_litre if outcome else 0.0,
        feasible=bool(outcome and outcome.feasible),
        iterations=len(problem.database),
        seconds=time.perf_counter() - began,
        candidates=points,
        scenario=scenario,
        reason=reason,
    )


def _extract(
    scenario: Any, candidates: list[LatticePoint], fallback: Design
) -> tuple[Design | None, LatticePoint | None]:
    """Recover the chosen lattice point and the sub-problem's design.

    The main problem's solution vector is the one-hot selection; the continuous
    variables live in the sub-scenario the adapter wrapped, so both halves have
    to be read back and reassembled into an eleven-variable design.
    """
    problem = scenario.formulation.optimization_problem
    solution = problem.solution
    if solution is None or solution.x_opt is None:
        return None, None

    selection = np.asarray(solution.x_opt, dtype=float).ravel()
    if selection.size < len(candidates):
        return None, None
    point = candidates[int(np.argmax(selection[: len(candidates)]))]

    adapter = scenario.formulation.sub_problem_scenario_adapter
    sub = adapter.scenario.formulation.optimization_problem
    sub_solution = sub.solution
    if sub_solution is None or sub_solution.x_opt is None:
        return fallback.replace(I=point.inter_axle), point

    values = dict(
        zip(
            sub.design_space.variable_names,
            np.asarray(sub_solution.x_opt, dtype=float).ravel(),
            strict=False,
        )
    )
    base = fallback.to_array()
    vector = [
        point.inter_axle if name == "I" else values.get(name, base[index])
        for index, name in enumerate(VARIABLE_NAMES)
    ]
    return Design.from_array(vector), point


def exhaustive(
    initial: Design,
    candidates: list[LatticePoint] | None = None,
    bounds: Bounds | None = None,
    speed_rpm: float = DEFAULT_SPEED_RPM,
    sub_max_iter: int = 25,
    samples: int | None = None,
    **kwargs: Any,
) -> list[tuple[LatticePoint, float, bool]]:
    """Solve every candidate separately, to check what the decomposition returned.

    One full continuous solve per lattice point, which is what the previous
    enumeration did and what the decomposition is measured against.

    Args:
        initial: Starting design.
        candidates: Lattice points; derived from ``initial`` if omitted.
        bounds: The continuous design box.
        speed_rpm: Crankshaft speed [rev/min].
        sub_max_iter: SLSQP budget per point.
        samples: Crank angles per revolution.
        **kwargs: Forwarded to the scenario builder.

    Returns:
        ``(point, km/L, feasible)`` per candidate; ``km/L`` is zero when the
        continuous problem had no feasible solution there.
    """
    from .performance import evaluate
    from .scenarios import _best_design, build_range_scenario

    points = (
        candidates_from_design(initial, speed_rpm=speed_rpm)
        if candidates is None
        else candidates
    )
    box = Bounds.around(initial, relative=0.35) if bounds is None else bounds
    n_samples = samples or 360
    rows: list[tuple[LatticePoint, float, bool]] = []
    for point in points:
        pinned = initial.replace(I=point.inter_axle)
        try:
            scenario = build_range_scenario(
                bounds=box,
                initial=pinned,
                speed_rpm=speed_rpm,
                module=point.module,
                teeth=point.teeth,
                samples=n_samples,
                **kwargs,
            )
            scenario.execute(algo_name="SLSQP", max_iter=sub_max_iter)
        except Exception:
            rows.append((point, 0.0, False))
            continue
        design = _best_design(scenario, objective="neg_range", fallback=pinned)
        outcome = evaluate(
            design,
            speed_rpm=speed_rpm,
            module=point.module,
            teeth=point.teeth,
            samples=n_samples,
        )
        rows.append((point, outcome.km_per_litre, outcome.feasible))
    return rows


def format_result(result: MinlpResult, title: str = "bi-level outer approximation") -> str:
    """Render a :class:`MinlpResult` as an aligned table."""
    lines = [title, "=" * len(title), ""]
    lines.append(f"  candidates on the lattice: {len(result.candidates)}")
    for item in result.candidates:
        mark = "  <- chosen" if item == result.point else ""
        lines.append(f"    {item}{mark}")
    lines.append("")
    if result.point is not None:
        lines.append(f"  chosen      {result.point}")
        lines.append(f"  range       {result.km_per_litre:.1f} km/L")
        lines.append(f"  feasible    {result.feasible}")
    lines.append(f"  iterations  {result.iterations}")
    lines.append(f"  seconds     {result.seconds:.0f}")
    if result.reason:
        lines.append(f"  stopped     {result.reason}")
    return "\n".join(lines)
