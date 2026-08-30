"""GEMSEO problem formulations for the EX-link sizing problem.

The report's final formulation is

.. math::
    \\begin{cases}
    l_b \\le X \\le u_b \\\\
    \\min_X f(X) = (-\\eta, H, B)^T \\\\
    c(X) \\le 0 \\\\
    c_{eq}(X) = 0
    \\end{cases}

with ``c = (mra - 10, W - 0.985, g - 0.01, 10 - d, \\gamma - 0.02)^T`` and
``c_eq = (STE - 74, \\epsilon - 16)^T``.

Four entry points cover the workflow the report went through, in the order it
went through it:

:func:`maximise_efficiency`
    Single objective, all constraints, global search.  The report's "big
    population GA" step; :data:`DEFAULT_GLOBAL_ALGORITHM` is differential
    evolution, which -- unlike NSGA-II -- handles the equality constraints
    directly.

:func:`refine`
    Augmented Lagrangian polish from a starting design.  The report's last step,
    taken because the external penalty function is accurate only for small ``r``
    and badly conditioned when ``r`` is small.

:func:`pareto_front`
    Multi-objective search over all three objectives.  The report's MOEA step.

:func:`local_pareto`
    The trick that finally worked for the report: run the MOEA in a box shrunk
    around a design already known to be good, rather than over the whole space.

:func:`sweep_moving_limits`
    The report's *first* Pareto method: keep efficiency as the only objective
    and walk a moving upper limit on ``H`` and ``B`` down, re-solving each time.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np
from gemseo import create_scenario
from gemseo.algos.design_space import DesignSpace
from gemseo.algos.opt.factory import OptimizationLibraryFactory
from gemseo.algos.pareto.pareto_front import ParetoFront
from gemseo.scenarios.base_scenario import BaseScenario

from .constants import DEFAULT_SPEC, DEFAULT_TARGETS, DesignTargets, EngineSpec
from .design import GLOBAL_BOUNDS, VARIABLE_NAMES, Bounds, Design
from .disciplines import ExlinkDiscipline
from .kinematics import DEFAULT_SAMPLES
from .model import Analysis, analyse
from .reference import PUBLISHED_DESIGN

LOGGER = logging.getLogger(__name__)

DEFAULT_GLOBAL_ALGORITHM = "DIFFERENTIAL_EVOLUTION"
"""Global search used in place of the report's genetic algorithm."""

DEFAULT_LOCAL_ALGORITHM = "Augmented_Lagrangian_order_0"
"""Local refinement; the report's own final step."""

DEFAULT_MOEA = "PYMOO_NSGA2"
"""Multi-objective evolutionary algorithm, the report's MOEA."""

OBJECTIVE_NAMES: tuple[str, str, str] = ("neg_efficiency", "height", "width")
"""``f(X) = (-eta, H, B)^T``."""

INEQUALITY_OUTPUTS: tuple[str, ...] = (
    "rod_angle_margin",
    "compatibility_margin",
    "tdc_gap_margin",
    "clearance_margin",
    "side_load_margin",
)
"""Discipline outputs that must stay ``<= 0``."""

EQUALITY_OUTPUTS: tuple[str, ...] = ("stroke_error", "compression_ratio_error")
"""Discipline outputs that must equal 0."""

DEFAULT_EQUALITY_TOLERANCE: dict[str, float] = {
    "stroke_error": 0.05,
    "compression_ratio_error": 0.05,
}
"""Half-widths used when relaxing the equalities into pairs of inequalities.

An algorithm such as NSGA-II cannot take equality constraints.  Relaxing
``STE = 74`` to ``|STE - 74| <= 0.05 mm`` and ``epsilon = 16`` to
``|epsilon - 16| <= 0.05`` keeps the design within manufacturing tolerance of
the specification -- the same pragmatism the report applies to ``g``, which it
constrains to a tolerance rather than to zero.
"""

MOEA_EQUALITY_TOLERANCE: dict[str, float] = {
    "stroke_error": 1.0,
    "compression_ratio_error": 0.5,
}
"""Looser half-widths for the equalities during a multi-objective run.

See :func:`moea_targets` for why the multi-objective stage needs a relaxed
problem at all.
"""

MOEA_TDC_GAP: float = 1.0
"""Relaxed bound on ``g`` during a multi-objective run [mm]."""


def moea_targets(
    targets: DesignTargets = DEFAULT_TARGETS, tdc_gap: float = MOEA_TDC_GAP
) -> DesignTargets:
    """Relax ``g`` for a population-based run.

    Why this is necessary is worth stating plainly, because it is the same wall
    the report hit ("the results were still disappointing ... even with big
    population I still got solutions that were even worse than the ones I got
    with the gradient based methods").

    Sampling 2000 designs uniformly from a box around a known-good solution and
    testing each constraint on its own gives:

    ==================  ============================
    constraint          satisfied by
    ==================  ============================
    ``d >= 10``         96 % of analysable designs
    ``W <= 0.985``      76 %
    ``mra <= 10``       28 %
    ``|eps - 16|``      8 %
    ``gamma <= 0.02``   6 %
    ``|STE - 74|``      4 %
    ``g <= 0.01``       **0.1 %**
    ==================  ============================

    ``g``, the gap between the two top dead centres, is nominally an inequality
    but at 0.01 mm it is a *third equality in disguise*.  With three equalities
    in eleven variables the feasible set is a thin sheet, and the joint hit rate
    for a random population is of order 1e-7.  NSGA-II duly returns a "front" of
    one point, because one point is all it ever found that was feasible.

    So the multi-objective stage runs on a relaxed problem: it is a scouting
    device that maps the *shape* of the efficiency-versus-size trade-off.  The
    design finally quoted comes from :func:`refine`, which drives ``g``, ``STE``
    and ``epsilon`` back onto their true targets -- exactly the two-stage
    workflow the report ends up with::

        front = pareto_front()
        final = refine(front.design)

    On this problem :func:`sweep_moving_limits` is the more reliable way to trace
    the trade-off, since every one of its steps is a warm-started local solve on
    the *unrelaxed* problem.

    Args:
        targets: The targets to relax.
        tdc_gap: The relaxed bound on ``g`` [mm].

    Returns:
        A copy of ``targets`` with ``max_tdc_gap`` widened.
    """
    return replace(targets, max_tdc_gap=tdc_gap)


@dataclass
class Outcome:
    """Result of one optimization run."""

    design: Design
    """The best design found (for multi-objective runs, the highest-efficiency
    point of the front)."""

    analysis: Analysis
    """Its full analysis."""

    algorithm: str
    scenario: BaseScenario = field(repr=False)
    front: list[Design] = field(default_factory=list, repr=False)
    """Pareto-optimal designs; empty for single-objective runs."""

    @property
    def feasible(self) -> bool:
        """Whether the returned design satisfies every constraint."""
        return is_feasible(self.analysis)

    @property
    def n_evaluations(self) -> int:
        """Design points recorded in the top-level problem database.

        For the augmented Lagrangian this counts *outer* iterations only: each
        one solves a sub-problem with its own database, so the true number of
        analyses is larger.
        """
        return len(self.scenario.formulation.optimization_problem.database)

    def summary(self) -> str:
        """A short human-readable report of the outcome."""
        return format_analysis(self.analysis, title=f"{self.algorithm} result")


INEQUALITY_TOLERANCE = 1.0e-4
"""Slack allowed on the inequality constraints when judging feasibility.

GEMSEO's own default (``ineq_tolerance``), and it matters here: an augmented
Lagrangian drives the active constraints *onto* their bounds, so a converged
design routinely lands a few parts in 1e5 outside one of them -- ``g``, in
practice, which the optimizer pushes hard against 0.01 mm.  Demanding ``<= 0``
exactly would label designs infeasible that every solver in the chain considers
converged.
"""

EQUALITY_TOLERANCE = 0.05
"""Half-width allowed on ``STE - 74`` [mm] and ``epsilon - 16``.

Looser than GEMSEO's 0.01 default because both are physical dimensions: 0.05 mm
on a 74 mm stroke is well inside what the linkage could be manufactured to.
"""


def is_feasible(
    analysis: Analysis,
    targets: DesignTargets = DEFAULT_TARGETS,
    tolerance: float = EQUALITY_TOLERANCE,
    inequality_tolerance: float = INEQUALITY_TOLERANCE,
) -> bool:
    """Whether an analysed design satisfies every constraint.

    Args:
        analysis: The analysis to check.
        targets: Constraint right-hand sides.
        tolerance: Half-width allowed on the two equality constraints.
        inequality_tolerance: Slack allowed on the inequality constraints; see
            :data:`INEQUALITY_TOLERANCE`.

    Returns:
        ``True`` if the design was analysable and meets every constraint within
        those tolerances.
    """
    from .model import equality_constraints, inequality_constraints

    if not analysis.valid:
        return False
    if np.any(inequality_constraints(analysis, targets) > inequality_tolerance):
        return False
    return bool(np.all(np.abs(equality_constraints(analysis, targets)) <= tolerance))


def build_design_space(
    bounds: Bounds = GLOBAL_BOUNDS, initial: Design | None = None
) -> DesignSpace:
    """Build the GEMSEO design space for the eleven variables.

    Args:
        bounds: The box ``l_b <= X <= u_b``.
        initial: Starting point; defaults to the centre of the box.

    Returns:
        A design space with one scalar variable per entry of
        :data:`exlink.design.VARIABLE_NAMES`.
    """
    if initial is None:
        start = 0.5 * (bounds.lower + bounds.upper)
    else:
        start = np.clip(initial.to_array(), bounds.lower, bounds.upper)

    space = DesignSpace()
    for index, name in enumerate(VARIABLE_NAMES):
        space.add_variable(
            name,
            lower_bound=float(bounds.lower[index]),
            upper_bound=float(bounds.upper[index]),
            value=float(start[index]),
        )
    return space


def _attach_constraints(
    scenario: BaseScenario,
    relax_equalities: bool,
    equality_tolerance: dict[str, float] | None,
    max_height: float,
    max_width: float,
) -> None:
    """Attach the report's constraints to a scenario."""
    for name in INEQUALITY_OUTPUTS:
        scenario.add_constraint(name, constraint_type="ineq")

    tolerance = dict(DEFAULT_EQUALITY_TOLERANCE)
    if equality_tolerance:
        tolerance.update(equality_tolerance)

    for name in EQUALITY_OUTPUTS:
        if relax_equalities:
            half_width = tolerance[name]
            # |residual| <= half_width, written as two one-sided inequalities.
            scenario.add_constraint(
                name,
                constraint_type="ineq",
                value=half_width,
                constraint_name=f"{name}_upper",
            )
            scenario.add_constraint(
                name,
                constraint_type="ineq",
                value=-half_width,
                positive=True,
                constraint_name=f"{name}_lower",
            )
        else:
            scenario.add_constraint(name, constraint_type="eq")

    if np.isfinite(max_height):
        scenario.add_constraint(
            "height", constraint_type="ineq", value=max_height, constraint_name="height_limit"
        )
    if np.isfinite(max_width):
        scenario.add_constraint(
            "width", constraint_type="ineq", value=max_width, constraint_name="width_limit"
        )


def build_scenario(
    objective: str | Sequence[str] = "neg_efficiency",
    bounds: Bounds = GLOBAL_BOUNDS,
    initial: Design | None = None,
    samples: int = DEFAULT_SAMPLES,
    spec: EngineSpec = DEFAULT_SPEC,
    targets: DesignTargets = DEFAULT_TARGETS,
    relax_equalities: bool = False,
    equality_tolerance: dict[str, float] | None = None,
    max_height: float = float("inf"),
    max_width: float = float("inf"),
) -> BaseScenario:
    """Assemble a constrained GEMSEO scenario.

    Args:
        objective: Output name, or several for a multi-objective run.
        bounds: The design box.
        initial: Starting design.
        samples: Crank angles per revolution.
        spec: Fixed engine data.
        targets: Constraint right-hand sides.
        relax_equalities: Express ``STE = 74`` and ``epsilon = 16`` as pairs of
            inequalities.  Required for algorithms that cannot take equality
            constraints, such as NSGA-II.
        equality_tolerance: Overrides for :data:`DEFAULT_EQUALITY_TOLERANCE`.
        max_height: Moving upper limit on ``H`` [mm].
        max_width: Moving upper limit on ``B`` [mm].

    Returns:
        A scenario ready to :meth:`~gemseo.scenarios.base_scenario.BaseScenario.execute`.
    """
    discipline = ExlinkDiscipline(samples=samples, spec=spec, targets=targets)
    space = build_design_space(bounds, initial)
    scenario = create_scenario(
        [discipline],
        objective,
        space,
        formulation_name="DisciplinaryOpt",
    )
    _attach_constraints(scenario, relax_equalities, equality_tolerance, max_height, max_width)
    return scenario


def _best_design(scenario: BaseScenario, objective: str = "neg_efficiency") -> Design:
    """Extract the best design from a solved single-objective scenario."""
    problem = scenario.formulation.optimization_problem
    solution = problem.solution
    if solution is not None and solution.x_opt is not None:
        return Design.from_array(solution.x_opt)
    return Design.from_array(problem.database.get_x_vect(-1))


def _run(
    scenario: BaseScenario,
    algorithm: str,
    settings: dict[str, Any],
    samples: int,
    spec: EngineSpec,
    targets: DesignTargets,
) -> Outcome:
    """Execute a scenario and package the outcome."""
    scenario.execute(algo_name=algorithm, **settings)
    design = _best_design(scenario)
    return Outcome(
        design=design,
        analysis=analyse(design, samples=samples, spec=spec, targets=targets),
        algorithm=algorithm,
        scenario=scenario,
    )


def maximise_efficiency(
    algorithm: str = DEFAULT_GLOBAL_ALGORITHM,
    bounds: Bounds = GLOBAL_BOUNDS,
    initial: Design | None = None,
    max_iter: int = 4000,
    samples: int = DEFAULT_SAMPLES,
    spec: EngineSpec = DEFAULT_SPEC,
    targets: DesignTargets = DEFAULT_TARGETS,
    max_height: float = float("inf"),
    max_width: float = float("inf"),
    relax_equalities: bool = False,
    seed: int | None = 1,
    **settings: Any,
) -> Outcome:
    """Maximise ``eta`` subject to every constraint.

    This is the report's single-objective step, with the size objectives either
    dropped or turned into the moving limits ``max_height`` / ``max_width``.

    Args:
        algorithm: Any GEMSEO optimizer; the default is differential evolution.
        bounds: The design box.
        initial: Starting design.
        max_iter: Evaluation budget.
        samples: Crank angles per revolution.
        spec: Fixed engine data.
        targets: Constraint right-hand sides.
        max_height: Moving upper limit on ``H`` [mm].
        max_width: Moving upper limit on ``B`` [mm].
        relax_equalities: See :func:`build_scenario`.
        seed: Random seed, for the stochastic algorithms that accept one.
        **settings: Extra algorithm settings passed through to GEMSEO.

    Returns:
        The outcome, whose :attr:`Outcome.feasible` flag should be checked.
    """
    scenario = build_scenario(
        "neg_efficiency",
        bounds=bounds,
        initial=initial,
        samples=samples,
        spec=spec,
        targets=targets,
        relax_equalities=relax_equalities,
        max_height=max_height,
        max_width=max_width,
    )
    options: dict[str, Any] = {"max_iter": max_iter}
    if seed is not None and algorithm in {"DIFFERENTIAL_EVOLUTION", "DUAL_ANNEALING"}:
        options["seed"] = seed
    if algorithm.startswith("Augmented_Lagrangian"):
        options["sub_algorithm_name"] = "NELDER-MEAD"
        options["sub_algorithm_settings"] = {"max_iter": 300}
    options.update(settings)
    return _run(scenario, algorithm, options, samples, spec, targets)


def refine(
    design: Design = PUBLISHED_DESIGN,
    algorithm: str = DEFAULT_LOCAL_ALGORITHM,
    bounds: Bounds | None = None,
    relative: float = 0.25,
    max_iter: int = 400,
    samples: int = DEFAULT_SAMPLES,
    spec: EngineSpec = DEFAULT_SPEC,
    targets: DesignTargets = DEFAULT_TARGETS,
    sub_algorithm_name: str = "NELDER-MEAD",
    sub_algorithm_settings: dict[str, Any] | None = None,
    **settings: Any,
) -> Outcome:
    """Polish a design with the augmented Lagrangian, as the report's last step.

    The report reaches its published solution by taking a point off the Pareto
    front and running an augmented Lagrangian from it, precisely because the
    external penalty function cannot be both accurate and well conditioned.

    Args:
        design: The starting design.
        algorithm: Local algorithm; the default is the augmented Lagrangian.
        bounds: The box to search in; defaults to a window around ``design``.
        relative: Half-width of that default window, as a fraction of ``|X_0|``.
        max_iter: Outer-iteration budget.
        samples: Crank angles per revolution.
        spec: Fixed engine data.
        targets: Constraint right-hand sides.
        sub_algorithm_name: Solver for the augmented Lagrangian sub-problems.
        sub_algorithm_settings: Settings for that solver.
        **settings: Extra algorithm settings passed through to GEMSEO.

    Returns:
        The refined outcome.
    """
    box = Bounds.around(design, relative=relative) if bounds is None else bounds
    scenario = build_scenario(
        "neg_efficiency",
        bounds=box,
        initial=design,
        samples=samples,
        spec=spec,
        targets=targets,
    )
    options: dict[str, Any] = {"max_iter": max_iter}
    if algorithm.startswith("Augmented_Lagrangian"):
        options["sub_algorithm_name"] = sub_algorithm_name
        options["sub_algorithm_settings"] = sub_algorithm_settings or {"max_iter": 300}
    options.update(settings)
    return _run(scenario, algorithm, options, samples, spec, targets)


def _supported_settings(algorithm: str, candidates: dict[str, Any]) -> dict[str, Any]:
    """Keep only the settings the given GEMSEO algorithm actually declares.

    Lets one set of "run the full budget" defaults be offered to every
    multi-objective algorithm without breaking those that do not have, say, a
    hypervolume criterion.

    Args:
        algorithm: The GEMSEO algorithm name.
        candidates: Settings to filter.

    Returns:
        The subset of ``candidates`` accepted by that algorithm.
    """
    factory = OptimizationLibraryFactory()
    library = factory.get_class(factory.algo_names_to_libraries[algorithm])
    fields = set(library.ALGORITHM_INFOS[algorithm].Settings.model_fields)
    return {name: value for name, value in candidates.items() if name in fields}


def _extract_front(scenario: BaseScenario) -> list[Design]:
    """Collect the Pareto-optimal designs of a solved multi-objective scenario.

    GEMSEO builds the front from the *feasible* points in the database, so an
    empty result means the run never found a feasible design -- see
    :data:`MOEA_EQUALITY_TOLERANCE`.
    """
    problem = scenario.formulation.optimization_problem
    pareto = ParetoFront.from_optimization_problem(problem)
    inputs = np.atleast_2d(np.asarray(pareto.x_optima, dtype=float))
    if inputs.size == 0:
        return []
    return [Design.from_array(row) for row in inputs]


def pareto_front(
    algorithm: str = DEFAULT_MOEA,
    bounds: Bounds = GLOBAL_BOUNDS,
    initial: Design | None = None,
    pop_size: int = 200,
    max_gen: int = 60,
    samples: int = 360,
    spec: EngineSpec = DEFAULT_SPEC,
    targets: DesignTargets | None = None,
    equality_tolerance: dict[str, float] | None = None,
    seed: int = 1,
    **settings: Any,
) -> Outcome:
    """Approximate the Pareto front of ``(-eta, H, B)``.

    Args:
        algorithm: A multi-objective algorithm; NSGA-II by default.
        bounds: The design box.
        initial: Starting design, used only by algorithms that take one.
        pop_size: Population size.  The report needed at least 550 individuals
            over the full space and grew to 2492; start smaller and grow.
        max_gen: Generation budget.  On this problem it matters more than
            ``pop_size``: seeded in a box around a good design, 20 generations
            still return a single point whatever the population, because the
            run has not yet worked its way into the thin feasible region, while
            35 generations give fronts of several dozen.  Spend extra budget
            here before spending it on population.
        samples: Crank angles per revolution.  Lowered by default, because a
            MOEA spends tens of thousands of evaluations here.
        spec: Fixed engine data.
        targets: Constraint right-hand sides; defaults to :func:`moea_targets`,
            which relaxes ``g``.  Read its docstring before overriding this.
        equality_tolerance: Half-widths for the relaxed equalities; defaults to
            :data:`MOEA_EQUALITY_TOLERANCE`.
        seed: Random seed.
        **settings: Extra algorithm settings passed through to GEMSEO.

    Returns:
        The outcome, with :attr:`Outcome.front` populated and
        :attr:`Outcome.design` set to the most efficient point of the front.
        Its designs satisfy the *relaxed* problem, so
        :attr:`Outcome.feasible` will usually be ``False`` until the chosen
        point has been through :func:`refine`.
    """
    targets = moea_targets() if targets is None else targets
    scenario = build_scenario(
        list(OBJECTIVE_NAMES),
        bounds=bounds,
        initial=initial,
        samples=samples,
        spec=spec,
        targets=targets,
        relax_equalities=True,
        equality_tolerance=equality_tolerance or MOEA_EQUALITY_TOLERANCE,
    )
    # GEMSEO's default convergence tests stop NSGA-II after about five
    # generations on this problem -- long before a front has formed -- because
    # the objectives barely move while the population is still hunting for
    # feasible designs.  Disable them so the generation budget is the budget
    # that actually applies, keeping only the ones this algorithm declares.
    options: dict[str, Any] = _supported_settings(
        algorithm,
        {
            "ftol_rel": 0.0,
            "ftol_abs": 0.0,
            "xtol_rel": 0.0,
            "xtol_abs": 0.0,
            "hv_tol_rel": 0.0,
            "hv_tol_abs": 0.0,
            "stop_crit_n_hv": max_gen + 1,
            "pop_size": pop_size,
            "max_gen": max_gen,
            "seed": seed,
        },
    )
    options["max_iter"] = pop_size * (max_gen + 1)
    options.update(settings)
    scenario.execute(algo_name=algorithm, **options)

    front = _extract_front(scenario)
    if not front:
        LOGGER.warning(
            "%s found no feasible design in %d evaluations; widen "
            "equality_tolerance, enlarge the population, or seed the run with "
            "local_pareto()",
            algorithm,
            len(scenario.formulation.optimization_problem.database),
        )
    analyses = [analyse(d, samples=samples, spec=spec, targets=targets) for d in front]
    if analyses:
        best_index = int(np.argmax([a.metrics.efficiency for a in analyses]))
        best_design, best_analysis = front[best_index], analyses[best_index]
    else:  # pragma: no cover - only if the algorithm returns nothing
        best_design = _best_design(scenario)
        best_analysis = analyse(best_design, samples=samples, spec=spec, targets=targets)

    return Outcome(
        design=best_design,
        analysis=best_analysis,
        algorithm=algorithm,
        scenario=scenario,
        front=front,
    )


def local_pareto(
    design: Design,
    relative: float = 0.1,
    absolute_angle: float = 20.0,
    **kwargs: Any,
) -> Outcome:
    """Run the MOEA in a box shrunk around a known-good design.

    This is the move that made the report's multi-objective step work: over the
    full eleven-dimensional space the MOEA returned solutions worse than the
    gradient-based ones, but seeded near an existing optimum it produced a
    usable local front.  Combining several such fronts into one starting
    population, then re-running over the full box, is what eventually gave the
    published front.

    Args:
        design: The centre of the shrunken box.
        relative: Half-width as a fraction of ``|X_0|`` for the lengths.
        absolute_angle: Half-width in degrees for ``theta_f`` and ``theta_r``.
        **kwargs: Forwarded to :func:`pareto_front`.

    Returns:
        The outcome of the local multi-objective run.
    """
    box = Bounds.around(design, relative=relative, absolute_angle=absolute_angle)
    return pareto_front(bounds=box, initial=design, **kwargs)


def sweep_moving_limits(
    limits: Iterable[float],
    algorithm: str = DEFAULT_LOCAL_ALGORITHM,
    on_height: bool = True,
    **kwargs: Any,
) -> list[Outcome]:
    """Trace a Pareto front by walking a moving size limit downwards.

    The report's first approach to the multi-objective problem: treat ``H`` (or
    ``B``) as a constraint rather than an objective, solve for maximum
    efficiency, tighten the limit, and repeat.  Each solve is an ordinary
    single-objective problem, and the sequence of solutions traces the
    efficiency-versus-size trade-off.

    Args:
        limits: The successive upper limits [mm], typically decreasing.
        algorithm: The single-objective algorithm to use at each step.  The
            default is the augmented Lagrangian: each step has to satisfy the
            two equalities *and* the relaxed-nothing inequalities exactly, and a
            derivative-free method such as COBYLA does not get there reliably
            within a sensible budget.
        on_height: Apply the limit to ``H``; otherwise to ``B``.
        **kwargs: Forwarded to :func:`maximise_efficiency`.

    Returns:
        One outcome per limit, in the order given.
    """
    outcomes: list[Outcome] = []
    previous: Design | None = kwargs.pop("initial", None)
    for limit in limits:
        limit_settings: dict[str, Any] = (
            {"max_height": limit} if on_height else {"max_width": limit}
        )
        outcome = maximise_efficiency(
            algorithm=algorithm, initial=previous, **limit_settings, **kwargs
        )
        outcomes.append(outcome)
        if outcome.feasible:
            # Warm-start the next, tighter problem from this solution.
            previous = outcome.design
        LOGGER.info(
            "moving limit %s <= %.1f mm -> eta = %.4f (%s)",
            "H" if on_height else "B",
            limit,
            outcome.analysis.metrics.efficiency,
            "feasible" if outcome.feasible else "infeasible",
        )
    return outcomes


def format_analysis(analysis: Analysis, title: str = "analysis") -> str:
    """Render an analysis as an aligned text table."""
    m = analysis.metrics
    t = DEFAULT_TARGETS
    lines = [title, "=" * len(title), "", str(analysis.design), ""]
    if not m.valid:
        lines.append(f"PENALISED: {m.reason}")
        return "\n".join(lines)
    rows = [
        ("eta   efficiency", f"{100 * m.efficiency:8.3f} %", "maximise"),
        ("H     height", f"{m.height:8.2f} mm", "minimise"),
        ("B     width", f"{m.width:8.2f} mm", "minimise"),
        ("STE   expansion stroke", f"{m.expansion_stroke:8.3f} mm", f"= {t.expansion_stroke}"),
        ("STC   compression stroke", f"{m.compression_stroke:8.3f} mm", ""),
        ("eps   compression ratio", f"{m.compression_ratio:8.3f}", f"= {t.compression_ratio}"),
        ("mra   piston rod angle", f"{m.rod_angle:8.3f} deg", f"<= {t.max_rod_angle}"),
        ("W     compatibility", f"{m.compatibility:8.4f}", f"<= {t.max_transmission}"),
        ("g     TDC gap", f"{m.tdc_gap:8.4f} mm", f"<= {t.max_tdc_gap}"),
        ("d     clearance", f"{m.clearance:8.2f} mm", f">= {t.min_clearance}"),
        ("gamma side load ratio", f"{m.side_load_ratio:8.4f}", f"<= {t.max_side_load}"),
        ("      mean torque", f"{m.mean_torque:8.1f} N.mm", ""),
    ]
    for label, value, target in rows:
        lines.append(f"  {label:<26} {value:>14}   {target}")
    lines.append("")
    lines.append(f"  feasible: {is_feasible(analysis)}")
    return "\n".join(lines)
