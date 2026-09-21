"""The design box cut into boxes: applying ``gemseo-box-subdivision`` here.

What the method is
------------------
The `box-subdivision outer approximation
<https://simoneconiglio.github.io/gemseo-box-subdivision/>`_ cuts a design space
into a Cartesian grid of boxes and makes the choice of a box a categorical
variable.  A mixed-integer master decides which box to open next from the cuts
of the boxes already solved; a local solver does the rest inside it.
Exploration and local exploitation stay in two distinct levels, which is exactly
the split §4.6's multistart improvises by hand -- perturb, project, re-solve --
without a rule for where to go next or a bound saying when to stop.

On its own benchmark the method reaches the optimum of Rastrigin in five
dimensions from every starting point.  This module asks whether it transfers to
this engine, and the answer turns on one property of this design space that
Rastrigin does not have.

Why a box cannot simply be handed to a solver here
--------------------------------------------------
:func:`exlink.model.analyse` gates a design twice, and over 20 000 uniform
draws from :data:`~exlink.design.GLOBAL_BOUNDS`, **91.2 %** fail the first
gate -- the crankshaft cannot turn through a full revolution, ``W >= 1`` -- and
a further **2.8 %** fail the second, the piston motion not having four monotone
phases.  Six per cent are analysable and **none** is feasible.  What the
disciplines report on the other ninety-four is
:class:`~exlink.constants.PenaltyValues` -- ``eta = 0``, ``H = B = 1000`` -- and
the penalty is a *constant*, so on 94 % of the design box the objective and
five of the seven constraints are a flat plateau with an exactly zero
gradient.

The consequence for a box subdivision is immediate and fatal if left alone.
The default policy starts each sub-problem at the center of its box, which is
the right policy on a problem that has a value everywhere: it is inside the box
and independent of the order the boxes are visited.  Here the center of a box is
almost surely on the plateau, SLSQP terminates after **one** evaluation, and the
master builds that box's cut from a point at which nothing was computed.  Run
that way over a 2 x 2 x 3 grid, the method spends 65 analyses and reports no
design at all.

So the method does not fail on this problem for the reason a method usually
fails on a problem.  It fails because it assumes a box is somewhere a solver can
start, and here it is not.

The box start this module supplies
----------------------------------
``gemseo-box-subdivision`` hands that policy to the caller through
``scenario_adapter_cls``.  :func:`create_restoring_adapter_class` fills it with
three stages, run inside the box the master has just chosen and nowhere else:

1. **Probe.** A uniform sample of the box, kept if
   :attr:`~exlink.model.Analysis.valid`.  An analysis costs 0.29 ms at 180
   crank angles, so a few hundred probes of a box cost a fraction of one
   evaluation of the range problem, which costs 10.4 s.
2. **Descend to analysability**, when no probe lands there.  ``W`` is smooth on
   the plateau and its gradient is exact, so minimising it inside the box is an
   ordinary bound-constrained solve -- the one quantity of the penalised branch
   that still carries information.
3. **Restore.** From the best analysable probe, maximise the *smallest* scaled
   constraint margin, in the epigraph form :mod:`exlink.restoration` uses, so
   the start is not merely feasible but interior.  With one addition that
   matters more than it sounds: the restoration is free to step back onto the
   plateau, where it is stuck for good, so the solve keeps the best
   **analysable** iterate it ever saw rather than the point SLSQP stops at.
   Without that bookkeeping, five restorations in nine ended on the plateau,
   worse than the probe they started from.

Only then does the box's sub-problem run, from a point inside its own box that
the disciplines can evaluate.

What this does not repair
-------------------------
The feasible set is a codimension-two manifold -- ``STE = 74`` and
``epsilon = 16``, relaxed to bands of 0.05 -- and restoration heads for the
nearest interior point, not the best one.  A box whose feasible sliver
restoration does not find is reported as an infeasible box, and the master cuts
it away on that evidence.  That is a true statement about this run, not about
the box, and the difference is not small: on the box that holds the best design
found anywhere, **one** of the sixty iterates the sub-problem passes through
lands inside the band, and whether it does is not stable across seeds.

§5.11 reports what the method is worth here once all of that is paid for.  The
short version is that it runs, it finds designs, and on this problem it does
not reach anything the single SLSQP solve of §4.2 had not already reached.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from .constants import DEFAULT_SPEC, DEFAULT_TARGETS, DesignTargets, EngineSpec
from .design import GLOBAL_BOUNDS, VARIABLE_NAMES, Bounds, Design
from .model import analyse, equality_constraints, inequality_constraints
from .scenarios import (
    DEFAULT_EQUALITY_TOLERANCE,
    DEFAULT_SAMPLES,
    EQUALITY_OUTPUTS,
    INEQUALITY_OUTPUTS,
    build_design_space,
)

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from gemseo.core.discipline import Discipline

FloatArray = NDArray[np.float64]

#: One characteristic magnitude per inequality constraint, in its own units.
#:
#: The five inequalities differ by four orders of magnitude -- a clearance
#: violation is tens of millimetres, a side-load violation is hundredths --
#: so the *smallest* margin of an unscaled vector is always the same row, and
#: maximising it ignores the other four.  Each magnitude below is the constraint's
#: own right-hand side from :class:`~exlink.constants.DesignTargets`, which is
#: the scale on which that constraint was written to matter.
CONSTRAINT_SCALES: dict[str, float] = {
    "rod_angle_margin": 10.0,
    "compatibility_margin": 0.1,
    "tdc_gap_margin": 0.01,
    "clearance_margin": 10.0,
    "side_load_margin": 0.02,
}

#: Floor applied to a scaled margin, in units of that margin's scale.
#:
#: An unanalysable design reports ``g = 1000`` against a bound of 0.01, so its
#: scaled margin is -1e5 and swamps every real row: two plateau points then
#: compare equal on the only number the probe ranks them by, and the ranking
#: is noise.  Clamping makes the plateau one flat value that any analysable
#: point beats, which is the only comparison that has to be right.
MARGIN_FLOOR = -50.0

#: Smallest scaled margin counted as *inside* the feasible set rather than on it.
INTERIOR_MARGIN = 1.0e-3

#: Designs drawn per box before the descent to analysability is attempted.
DEFAULT_PROBES = 512

#: Analysable probes carried forward into a restoration, best margin first.
DEFAULT_RESTARTS = 3


def scaled_margins(
    design: Design,
    samples: int = DEFAULT_SAMPLES,
    spec: EngineSpec = DEFAULT_SPEC,
    targets: DesignTargets = DEFAULT_TARGETS,
    tolerance: Mapping[str, float] = DEFAULT_EQUALITY_TOLERANCE,
) -> FloatArray:
    """Every constraint of the geometric problem as a margin, on one scale.

    Positive means satisfied, and a margin of one means the constraint has a
    whole characteristic magnitude of room.  The two equalities contribute their
    two sides separately, so the vector has nine rows for seven constraints.

    Args:
        design: The design to measure.
        samples: Crank angles per revolution.
        spec: Fixed engine data.
        targets: Constraint right-hand sides.
        tolerance: Half-widths of the two relaxed equality bands.

    Returns:
        Nine scaled margins, floored at :data:`MARGIN_FLOOR`.
    """
    analysis = analyse(design, samples=samples, spec=spec, targets=targets)
    return scaled_margins_of(analysis, targets, tolerance)


def scaled_margins_of(
    analysis: Any,
    targets: DesignTargets = DEFAULT_TARGETS,
    tolerance: Mapping[str, float] = DEFAULT_EQUALITY_TOLERANCE,
) -> FloatArray:
    """The margins of an analysis already computed; see :func:`scaled_margins`."""
    rows = [
        -value / CONSTRAINT_SCALES[name]
        for name, value in zip(
            INEQUALITY_OUTPUTS,
            inequality_constraints(analysis, targets),
            strict=True,
        )
    ]
    for name, residual in zip(
        EQUALITY_OUTPUTS, equality_constraints(analysis, targets), strict=True
    ):
        half_width = tolerance[name]
        rows.append((half_width - residual) / half_width)
        rows.append((half_width + residual) / half_width)
    return np.maximum(np.array(rows, dtype=float), MARGIN_FLOOR)


@dataclass(frozen=True)
class BoxStart:
    """Where one box's sub-problem was told to start, and how that was found."""

    design: Design
    margin: float
    """Smallest scaled margin reached; positive means feasible."""

    stage: str
    """``"probe"``, ``"descent"``, ``"restoration"`` or ``"centre"``."""

    evaluations: int

    @property
    def feasible(self) -> bool:
        """Whether the start is inside the feasible set, not merely analysable."""
        return self.margin >= INTERIOR_MARGIN


class _Geometry:
    """The geometric discipline, its margins and their exact Jacobian.

    One object rather than free functions because the Jacobian is only defined
    after the discipline has run at the same point, and because the evaluation
    count is what the comparisons in §5.11 are read against.
    """

    def __init__(
        self,
        samples: int,
        spec: EngineSpec,
        targets: DesignTargets,
        tolerance: Mapping[str, float],
    ) -> None:
        from .disciplines import ExlinkDiscipline

        self.discipline = ExlinkDiscipline(samples=samples, spec=spec, targets=targets)
        self.samples = samples
        self.spec = spec
        self.targets = targets
        self.tolerance = tolerance
        self.evaluations = 0

    def analyse(self, vector: FloatArray) -> Any:
        """Analyse a design vector, counting the call."""
        self.evaluations += 1
        return analyse(
            Design.from_array(vector),
            samples=self.samples,
            spec=self.spec,
            targets=self.targets,
        )

    def margins(self, vector: FloatArray) -> FloatArray:
        """The scaled margins at a design vector."""
        return scaled_margins_of(self.analyse(vector), self.targets, self.tolerance)

    def margins_jacobian(self, vector: FloatArray) -> FloatArray:
        """The exact Jacobian of :meth:`margins`, nine rows by eleven columns."""
        mapping = Design.from_array(vector).to_mapping()
        self.discipline.execute(mapping)
        jacobian = self.discipline.linearize(mapping, compute_all_jacobians=True)

        def row(name: str) -> FloatArray:
            return np.array(
                [float(np.ravel(jacobian[name][v])[0]) for v in VARIABLE_NAMES]
            )

        rows = [-row(name) / CONSTRAINT_SCALES[name] for name in INEQUALITY_OUTPUTS]
        for name in EQUALITY_OUTPUTS:
            half_width = self.tolerance[name]
            gradient = row(name)
            rows.append(-gradient / half_width)
            rows.append(gradient / half_width)
        return np.array(rows)

    def compatibility(self, vector: FloatArray) -> float:
        """``W - C``, the first gate, which is smooth on the penalised branch."""
        output = self.discipline.execute(Design.from_array(vector).to_mapping())
        self.evaluations += 1
        return float(np.ravel(output["compatibility_margin"])[0])

    def compatibility_jacobian(self, vector: FloatArray) -> FloatArray:
        """The exact gradient of :meth:`compatibility`."""
        mapping = Design.from_array(vector).to_mapping()
        self.discipline.execute(mapping)
        jacobian = self.discipline.linearize(mapping, compute_all_jacobians=True)
        return np.array(
            [
                float(np.ravel(jacobian["compatibility_margin"][v])[0])
                for v in VARIABLE_NAMES
            ]
        )


def find_box_start(
    lower: FloatArray,
    upper: FloatArray,
    incumbent: Design | None = None,
    seed: int = 0,
    probes: int = DEFAULT_PROBES,
    restarts: int = DEFAULT_RESTARTS,
    restoration_iterations: int = 120,
    descent_iterations: int = 40,
    samples: int = DEFAULT_SAMPLES,
    spec: EngineSpec = DEFAULT_SPEC,
    targets: DesignTargets = DEFAULT_TARGETS,
    tolerance: Mapping[str, float] = DEFAULT_EQUALITY_TOLERANCE,
) -> BoxStart:
    """Find a point inside a box that a local solver can start from.

    The three stages of the module docstring, in order, stopping as soon as one
    of them produces a feasible point.  A box in which none of them does still
    returns its best analysable point, and a box with no analysable point at all
    returns its center -- from which the sub-problem will return nothing, which
    is the honest report that this run found the box empty.

    Args:
        lower: Lower corner of the box, over all eleven variables.
        upper: Upper corner of the box.
        incumbent: A design already known to be good, clipped into the box and
            tried before anything is drawn.  Clipping a feasible design into a
            box it is not in generally leaves the feasible set, so this is a
            start for the restoration rather than a start for the sub-problem;
            it is nonetheless the best one available, being the only point in
            the box whose eight unsubdivided components are known to work
            together.
        seed: Random seed, so a box gives the same start every time it is opened.
        probes: Designs drawn uniformly in the box.
        restarts: Analysable probes carried into a restoration.
        restoration_iterations: SLSQP budget per restoration.
        descent_iterations: SLSQP budget for the descent to analysability.
        samples: Crank angles per revolution.
        spec: Fixed engine data.
        targets: Constraint right-hand sides.
        tolerance: Half-widths of the two relaxed equality bands.

    Returns:
        The starting point, its margin and the stage that produced it.
    """
    geometry = _Geometry(samples, spec, targets, tolerance)
    rng = np.random.default_rng(seed)
    centre = 0.5 * (lower + upper)

    # Stage 1: probe.  The incumbent and the centre go in first, so that a box
    # which is fine at either costs one analysis rather than ``probes`` of them.
    seeded = [centre] if incumbent is None else [
        np.clip(incumbent.to_array(), lower, upper),
        centre,
    ]
    drawn = rng.uniform(lower, upper, size=(max(probes - len(seeded), 0), lower.size))
    candidates: list[tuple[float, FloatArray]] = []
    best_compatibility: tuple[float, FloatArray] | None = None
    for vector in (*seeded, *drawn):
        analysis = geometry.analyse(vector)
        if analysis.valid:
            margin = float(
                scaled_margins_of(analysis, targets, tolerance).min()
            )
            candidates.append((margin, np.array(vector)))
            if margin >= INTERIOR_MARGIN:
                return BoxStart(
                    Design.from_array(vector), margin, "probe", geometry.evaluations
                )
            continue
        value = float(analysis.metrics.compatibility)
        if np.isfinite(value) and (
            best_compatibility is None or value < best_compatibility[0]
        ):
            best_compatibility = (value, np.array(vector))

    # Stage 2: no probe was analysable, so descend on the one quantity that is
    # not flat there.
    if not candidates:
        if best_compatibility is None:
            return BoxStart(
                Design.from_array(centre), MARGIN_FLOOR, "centre", geometry.evaluations
            )
        vector = _descend_to_analysability(
            geometry, best_compatibility[1], lower, upper, descent_iterations
        )
        analysis = geometry.analyse(vector)
        if not analysis.valid:
            return BoxStart(
                Design.from_array(centre), MARGIN_FLOOR, "centre", geometry.evaluations
            )
        candidates.append((
            float(scaled_margins_of(analysis, targets, tolerance).min()),
            vector,
        ))
        stage = "descent"
    else:
        stage = "probe"

    # Stage 3: restore, from the most promising analysable points first.
    candidates.sort(key=lambda pair: -pair[0])
    best = BoxStart(
        Design.from_array(candidates[0][1]),
        candidates[0][0],
        stage,
        geometry.evaluations,
    )
    for _, vector in candidates[:restarts]:
        restored, margin = _restore_in_box(
            geometry, vector, lower, upper, restoration_iterations
        )
        if restored is None or margin <= best.margin:
            continue
        best = BoxStart(
            Design.from_array(restored), margin, "restoration", geometry.evaluations
        )
        if best.feasible:
            break
    return BoxStart(best.design, best.margin, best.stage, geometry.evaluations)


def _descend_to_analysability(
    geometry: _Geometry,
    start: FloatArray,
    lower: FloatArray,
    upper: FloatArray,
    max_iterations: int,
) -> FloatArray:
    """Minimise ``W`` inside the box, to get off the penalised branch."""
    from scipy.optimize import minimize

    outcome = minimize(
        geometry.compatibility,
        np.clip(start, lower, upper),
        jac=geometry.compatibility_jacobian,
        method="SLSQP",
        bounds=list(zip(lower, upper, strict=True)),
        options={"maxiter": int(max_iterations), "ftol": 1.0e-9},
    )
    return np.clip(np.asarray(outcome.x, dtype=float), lower, upper)


def _restore_in_box(
    geometry: _Geometry,
    start: FloatArray,
    lower: FloatArray,
    upper: FloatArray,
    max_iterations: int,
) -> tuple[FloatArray | None, float]:
    """Maximise the smallest scaled margin inside the box.

    The epigraph form of :mod:`exlink.restoration`: the scalar ``t`` becomes the
    twelfth unknown and the non-smooth ``min_i c_i`` moves into the constraints.

    What is added here is the bookkeeping.  A step that leaves the analysable
    region lands on a plateau where every margin is :data:`MARGIN_FLOOR` and
    every gradient is zero, and SLSQP does not come back.  So every iterate is
    scored as it is evaluated and the best **analysable** one is what the solve
    returns, which is not in general where it stopped.
    """
    from scipy.optimize import minimize

    best: list[Any] = [None, -np.inf]

    def margins(vector: FloatArray) -> FloatArray:
        clipped = np.clip(vector, lower, upper)
        analysis = geometry.analyse(clipped)
        values = scaled_margins_of(analysis, geometry.targets, geometry.tolerance)
        worst = float(values.min())
        if analysis.valid and worst > best[1]:
            best[0], best[1] = np.array(clipped), worst
        return values

    initial = margins(np.clip(start, lower, upper))
    augmented = np.concatenate([np.clip(start, lower, upper), [float(initial.min())]])
    minimize(
        lambda z: -z[-1],
        augmented,
        jac=lambda z: np.concatenate([np.zeros(z.size - 1), [-1.0]]),
        method="SLSQP",
        bounds=[*zip(lower, upper, strict=True), (MARGIN_FLOOR, 1.0e3)],
        constraints=[
            {
                "type": "ineq",
                "fun": lambda z: margins(z[:-1]) - z[-1],
                "jac": lambda z: np.hstack([
                    geometry.margins_jacobian(np.clip(z[:-1], lower, upper)),
                    -np.ones((initial.size, 1)),
                ]),
            }
        ],
        options={"maxiter": int(max_iterations), "ftol": 1.0e-12},
    )
    return best[0], float(best[1])


def create_restoring_adapter_class(
    subdivision: Any,
    bounds: Bounds = GLOBAL_BOUNDS,
    incumbent: Design | None = None,
    fixed: Mapping[str, float] = {},
    seed: int = 0,
    samples: int = DEFAULT_SAMPLES,
    **start_settings: Any,
) -> type:
    """Return a scenario adapter restoring a startable point inside each box.

    The class is built per run rather than parametrised, because the
    ``Benders`` formulation takes a *class* and instantiates it itself.

    Every box the master opens is restored once and cached under its one-hot
    key: the master re-proposes a box it has already solved, and paying for the
    probe again would charge the method for evaluations it does not need.

    Args:
        subdivision: The :class:`.BoxSubdivision` the scenario was built with.
        bounds: The design box of the original problem, supplying the range of
            the variables that are not subdivided.
        incumbent: A design already known to be good, tried in every box before
            anything is drawn there.
        fixed: Variables an outer loop has already chosen, and their values;
            they are held at that value throughout the restoration.
        seed: Random seed of the probe.
        samples: Crank angles per revolution.
        **start_settings: Forwarded to :func:`find_box_start`.

    Returns:
        The adapter class, for ``BoxSubdivisionScenario(scenario_adapter_cls=...)``.
    """
    from gemseo_bilevel_outer_approximation.disciplines.scenario_adapters.mdo_scenario_adapter_benders import (  # noqa: E501
        MDOScenarioAdapterBenders,
    )

    one_hot_names = subdivision.get_one_hot_names({})
    normalized_names = subdivision.get_normalized_names({})
    index_of = {name: index for index, name in enumerate(VARIABLE_NAMES)}

    class RestoringScenarioAdapter(MDOScenarioAdapterBenders):
        """A Benders adapter starting each box at a restored feasible point."""

        starts: dict[tuple[float, ...], BoxStart] = {}  # noqa: RUF012
        """Every box opened so far, keyed by the one-hot the master chose."""

        solutions: dict[tuple[float, ...], tuple[Design, float] | None] = {}  # noqa: RUF012
        """What each box's sub-problem returned, keyed the same way.

        A box is the unit this method reasons in, so the run has to be readable
        box by box: a master that stops early and a subdivision that had nothing
        to find are the same summary line and different failures.
        """

        current_key: tuple[float, ...] = ()
        """The box being solved, so ``_post_run`` knows what it is recording."""

        current_box: dict[str, tuple[float, float]] = {}  # noqa: RUF012
        """Bounds of the box being solved, per subdivided variable.

        The sub-problem solves in the normalized variables of its box, so
        reading a design back out of its database needs the box it was
        normalized against.
        """

        def _post_run(self) -> None:
            super()._post_run()
            key = type(self).current_key
            if not key or key in type(self).solutions:
                return
            base = (
                GLOBAL_BOUNDS.lower if incumbent is None else incumbent.to_array()
            ).copy()
            for name, value in fixed.items():
                base[index_of[name]] = float(value)
            type(self).solutions[key] = _best_feasible_solution(
                self.scenario.formulation.optimization_problem,
                type(self).current_box,
                index_of,
                base,
            )

        def _pre_run(self) -> None:
            super()._pre_run()
            design_space = self.scenario.formulation.optimization_problem.design_space
            data = self.io.data

            lower = np.array(bounds.lower, dtype=float)
            upper = np.array(bounds.upper, dtype=float)
            for name, value in fixed.items():
                lower[index_of[name]] = upper[index_of[name]] = float(value)

            key: list[float] = []
            boxed: dict[str, tuple[float, float]] = {}
            for variable_name, one_hot_name in one_hot_names.items():
                if one_hot_name not in data:
                    continue
                one_hot = np.asarray(data[one_hot_name], dtype=float)
                key.extend(one_hot.tolist())
                box_lower, box_upper = subdivision.compute_bounds(
                    variable_name, one_hot
                )
                index = index_of[variable_name]
                lower[index] = float(np.ravel(box_lower)[0])
                upper[index] = float(np.ravel(box_upper)[0])
                boxed[variable_name] = (lower[index], upper[index])

            cache_key = tuple(key)
            type(self).current_key = cache_key
            type(self).current_box = dict(boxed)
            start = type(self).starts.get(cache_key)
            if start is None:
                start = find_box_start(
                    lower,
                    upper,
                    incumbent=incumbent,
                    seed=seed,
                    samples=samples,
                    **start_settings,
                )
                type(self).starts[cache_key] = start

            vector = start.design.to_array()
            current: dict[str, Any] = {}
            for variable_name, (box_lower, box_upper) in boxed.items():
                value = float(vector[index_of[variable_name]])
                normalized_name = normalized_names[variable_name]
                if normalized_name in design_space:
                    # The normalized formulation solves for ``u`` in [0, 1],
                    # the mapping being ``x = l + u (u_b - l)``; invert it.
                    span = box_upper - box_lower
                    position = 0.5 if span == 0.0 else (value - box_lower) / span
                    current[normalized_name] = np.array([
                        float(np.clip(position, 0.0, 1.0))
                    ])
                elif variable_name in design_space:
                    # The constraint formulation keeps the variable itself, in
                    # millimetres, with the box imposed as a constraint.
                    current[variable_name] = np.array([value])
            for variable_name in VARIABLE_NAMES:
                # A variable the subdivision left alone stays itself, and its
                # restored value is the start the sub-problem wants.
                if variable_name in boxed or variable_name not in design_space:
                    continue
                current[variable_name] = np.array([
                    float(vector[index_of[variable_name]])
                ])
            if current:
                design_space.set_current_value(current)

    return RestoringScenarioAdapter


def _best_feasible_solution(
    problem: Any,
    box: Mapping[str, tuple[float, float]],
    index_of: Mapping[str, int],
    base: FloatArray | None = None,
) -> tuple[Design, float] | None:
    """The best feasible point of one box's sub-problem, in engineering units.

    Read from the sub-problem's own database rather than re-analysed, for two
    reasons.  Its objective is the one the master is served -- efficiency here,
    range in :func:`build_subdivided_range_scenario`, and the second is not a
    function of the geometry alone.  And a gradient solver on this problem
    routinely stops marginally outside a constraint, so the point it *reports*
    is the least infeasible one it saw rather than a feasible one; the database
    holds the feasible points it passed through on the way.

    Args:
        problem: The sub-problem, after it has been solved.
        box: Bounds of the box, per subdivided variable, for undoing the
            normalization the sub-problem solves in.
        index_of: Position of each variable in the design vector.
        base: Values of the variables the design space does not hold, such as
            the ``I`` the gear pair pinned; the reference design if omitted.

    Returns:
        The best feasible design and the objective value there, or ``None``
        when the box's sub-problem passed through no feasible point.
    """
    from .scenarios import _best_feasible_vector

    objective = problem.objective.name
    vector = _best_feasible_vector(problem, objective)
    if vector is None:
        return None

    from .reference import REFINED_DESIGN

    values = problem.design_space.convert_array_to_dict(vector)
    design = (
        REFINED_DESIGN.to_array() if base is None else np.array(base, dtype=float)
    )
    for name, value in values.items():
        scalar = float(np.ravel(value)[0])
        if name.endswith("_normalized"):
            variable = name[: -len("_normalized")]
            lower, upper = box[variable]
            design[index_of[variable]] = lower + scalar * (upper - lower)
        elif name in index_of:
            design[index_of[name]] = scalar

    entry = problem.database.get(vector)
    value = entry.get(objective) if entry is not None else None
    if value is None:
        return None
    return Design.from_array(design), -float(np.ravel(value)[0])


def _relaxed_equality_disciplines() -> tuple[list[Discipline], tuple[str, ...]]:
    """One discipline per side of each relaxed equality band.

    ``add_constraint(..., constraint_name=...)`` is the natural way to write
    ``|r| <= h`` as two inequalities, and it cannot be used here: the Benders
    adapter differentiates a sub-problem constraint by looking its *name* up
    among the discipline outputs, so a renamed constraint raises ``KeyError``
    on the first linearisation.  Each side therefore becomes an output in its
    own right, through a :class:`.LinearCombination` whose Jacobian is exact
    and constant.
    """
    from gemseo.disciplines.linear_combination import LinearCombination

    disciplines: list[Discipline] = []
    names: list[str] = []
    for name in EQUALITY_OUTPUTS:
        half_width = DEFAULT_EQUALITY_TOLERANCE[name]
        for suffix, sign in (("upper", 1.0), ("lower", -1.0)):
            output_name = f"{name}_{suffix}"
            disciplines.append(
                LinearCombination(
                    [name],
                    output_name,
                    input_coefficients={name: sign},
                    offset=-half_width,
                )
            )
            names.append(output_name)
    return disciplines, tuple(names)


def _negated_disciplines(
    output_names: Sequence[str],
) -> tuple[list[Discipline], tuple[str, ...]]:
    """One discipline per constraint written the wrong way round.

    ``runs_margin`` and ``gear_margin`` are *margins*: positive means the
    engine runs and the mesh fits, so they are attached with ``positive=True``.
    GEMSEO implements that by negating the function and naming the constraint
    ``-runs_margin``, and the Benders adapter then looks that name up among the
    discipline outputs and raises ``KeyError`` on the first linearisation --
    the same limitation the relaxed equalities hit, reached by a different
    route, and it does not appear until the master linearises the adapter, so
    a short run does not show it.

    Negating them in a discipline instead gives a violation, ``<= 0`` like
    every other constraint, under a name the adapter can find.
    """
    from gemseo.disciplines.linear_combination import LinearCombination

    disciplines: list[Discipline] = []
    names: list[str] = []
    for name in output_names:
        violation = name.removesuffix("_margin") + "_violation"
        disciplines.append(
            LinearCombination([name], violation, input_coefficients={name: -1.0})
        )
        names.append(violation)
    return disciplines, tuple(names)


def build_subdivided_scenario(
    n_subdivisions: Mapping[str, int],
    bounds: Bounds = GLOBAL_BOUNDS,
    initial: Design | None = None,
    samples: int = DEFAULT_SAMPLES,
    spec: EngineSpec = DEFAULT_SPEC,
    targets: DesignTargets = DEFAULT_TARGETS,
    settings: Any = None,
    seed: int = 0,
    restoring: bool = True,
    main_level: bool = True,
    **start_settings: Any,
) -> Any:
    """Pose the geometric efficiency problem as a box subdivision.

    The problem is §4.2's: maximise indicated efficiency over the eleven
    linkage variables, under the five inequalities and the two equalities,
    the latter relaxed to bands.  What changes is that the variables named in
    ``n_subdivisions`` are cut into boxes and a mixed-integer master chooses
    between them.

    Every constraint is attached with ``main_level=True``.  A box may hold no
    feasible design at all -- most of them hold none -- and the master has to
    be told so by a feasibility cut rather than left waiting for a sub-problem
    that cannot converge.

    Args:
        n_subdivisions: The variables to subdivide, and how finely.  Naming two
            or three is the scale the measurements were taken at; the master
            grows with the binaries, which is the sum of these numbers.
        bounds: The design box.
        initial: Starting design, used for the variables left unsubdivided.
        samples: Crank angles per revolution.
        spec: Fixed engine data.
        targets: Constraint right-hand sides.
        settings: ``BoxSubdivisionSettings`` or ``SweptBoxSubdivisionSettings``;
            the swept one if omitted, since the convexity margin is in the units
            of the objective and efficiency has no scale worth guessing.
        seed: Random seed of the box probes.
        restoring: Whether to restore a startable point in each box.  ``False``
            leaves the default policy -- the center of the box -- which is what
            §5.11 measures the restoration against.
        main_level: Whether an infeasible box produces a feasibility cut.
            ``False`` leaves the constraints on the sub-problem alone and is
            measured only to show what the master is actually steered by here.
        **start_settings: Forwarded to :func:`find_box_start`.

    Returns:
        A ``BoxSubdivisionScenario`` ready to execute.
    """
    from gemseo_box_subdivision import BoxSubdivisionScenario, SweptBoxSubdivisionSettings

    from .disciplines import ExlinkDiscipline
    from .reference import REFINED_DESIGN

    start = REFINED_DESIGN if initial is None else initial
    sides, side_names = _relaxed_equality_disciplines()
    geometry = ExlinkDiscipline(samples=samples, spec=spec, targets=targets)
    # The default cache keeps the last call only, and what a box subdivision has
    # to be read against is every design its sub-problems evaluated: the run
    # reports one number per box and the cost is the whole tour.
    geometry.set_cache(geometry.CacheType.MEMORY_FULL)
    disciplines: list[Discipline] = [geometry, *sides]
    design_space = build_design_space(bounds, start)
    subdivision_settings = (
        SweptBoxSubdivisionSettings() if settings is None else settings
    )

    adapter = None
    if restoring:
        from gemseo_box_subdivision.subdivisions.box import BoxSubdivision

        adapter = create_restoring_adapter_class(
            BoxSubdivision.from_design_space(
                design_space, n_subdivisions, tuple(n_subdivisions)
            ),
            bounds=bounds,
            incumbent=start,
            seed=seed,
            samples=samples,
            spec=spec,
            targets=targets,
            **start_settings,
        )

    scenario = BoxSubdivisionScenario(
        disciplines,
        "neg_efficiency",
        design_space,
        n_subdivisions=dict(n_subdivisions),
        settings=subdivision_settings,
        scenario_adapter_cls=adapter,
    )
    for name in (*INEQUALITY_OUTPUTS, *side_names):
        scenario.formulation.add_constraint(
            name, constraint_type="ineq", main_level=main_level
        )
    return scenario


@dataclass(frozen=True)
class SubdivisionOutcome:
    """What one subdivided run found, and what it spent finding it."""

    design: Design | None
    """The best feasible design found, or ``None`` if the run found none."""

    efficiency: float
    """Its indicated efficiency; ``nan`` when nothing feasible was found."""

    evaluations: int
    """Analyses charged to the sub-problems."""

    restoration_evaluations: int
    """Analyses charged to the box starts, which the run must also pay for."""

    boxes: int
    """Boxes the master opened."""

    feasible_boxes: int
    """Boxes whose sub-problem returned a design meeting every constraint."""

    seconds: float = 0.0

    per_box: tuple[tuple[tuple[int, ...], float], ...] = ()
    """Each box the master opened, and the efficiency it yielded.

    ``nan`` where the box yielded no feasible design.  A summary line cannot
    tell a master that stopped early from a subdivision that had nothing left
    to find; this can.
    """


def maximise_efficiency_by_subdivision(
    n_subdivisions: Mapping[str, int] | None = None,
    bounds: Bounds = GLOBAL_BOUNDS,
    initial: Design | None = None,
    samples: int = DEFAULT_SAMPLES,
    seed: int = 0,
    restoring: bool = True,
    main_level: bool = True,
    settings: Any = None,
    **start_settings: Any,
) -> SubdivisionOutcome:
    """Run the subdivided geometric problem and report what it found.

    Args:
        n_subdivisions: The variables to subdivide and how finely; ``q_1``,
            ``q_2`` and ``theta_f`` at 2, 2 and 3 if omitted, which is the grid
            §5.11 reports.
        bounds: The design box.
        initial: Starting design for the variables left unsubdivided.
        samples: Crank angles per revolution.
        seed: Random seed of the box probes.
        restoring: Whether to restore a startable point in each box.
        main_level: Whether an infeasible box produces a feasibility cut.
        settings: Settings of the subdivision.
        **start_settings: Forwarded to :func:`find_box_start`.

    Returns:
        The best feasible design found, with its cost.
    """
    import time

    from .scenarios import is_feasible

    grid = {"q_1": 2, "q_2": 2, "theta_f": 3} if n_subdivisions is None else dict(
        n_subdivisions
    )
    began = time.perf_counter()
    scenario = build_subdivided_scenario(
        grid,
        bounds=bounds,
        initial=initial,
        samples=samples,
        settings=settings,
        seed=seed,
        restoring=restoring,
        main_level=main_level,
        **start_settings,
    )
    scenario.execute()

    # The sub-problem's database is in the normalized variables of whichever
    # box was open at the time, so the designs are read back from the
    # discipline's own cache instead, where they are in engineering units.
    adapter = scenario.formulation.sub_problem_scenario_adapter
    cache = _designs_in(adapter.scenario)
    best_design, best_efficiency = None, float("-inf")
    for design in cache:
        analysis = analyse(design, samples=samples)
        if not is_feasible(analysis):
            continue
        efficiency = float(analysis.metrics.efficiency)
        if efficiency > best_efficiency:
            best_design, best_efficiency = design, efficiency

    starts = getattr(type(adapter), "starts", {})
    solutions = getattr(type(adapter), "solutions", {})
    per_box = tuple(
        (
            _decode_one_hot(key, grid),
            float("nan") if solutions.get(key) is None else solutions[key][1],
        )
        for key in starts
    )
    return SubdivisionOutcome(
        design=best_design,
        efficiency=float("nan") if best_design is None else best_efficiency,
        evaluations=len(cache),
        restoration_evaluations=sum(start.evaluations for start in starts.values()),
        boxes=len(starts),
        feasible_boxes=sum(1 for _, value in per_box if not np.isnan(value)),
        seconds=time.perf_counter() - began,
        per_box=per_box,
    )


@dataclass(frozen=True)
class BoxOutcome:
    """One box solved on its own, as the enumeration reference reports it."""

    index: tuple[int, ...]
    start: BoxStart
    design: Design | None
    efficiency: float
    evaluations: int

    @property
    def solved(self) -> bool:
        """Whether the box yielded a design meeting every constraint."""
        return self.design is not None


def enumerate_boxes(
    n_subdivisions: Mapping[str, int],
    bounds: Bounds = GLOBAL_BOUNDS,
    initial: Design | None = None,
    samples: int = DEFAULT_SAMPLES,
    spec: EngineSpec = DEFAULT_SPEC,
    targets: DesignTargets = DEFAULT_TARGETS,
    seed: int = 0,
    sub_problem_max_iter: int = 60,
    algorithm: str = "SLSQP",
    **start_settings: Any,
) -> list[BoxOutcome]:
    """Solve every box of the subdivision, which is what the master must beat.

    The reference the box-subdivision literature compares against, and the only
    way to tell a master that explored badly from a subdivision that had nothing
    to find.  It runs the same box-start policy and the same sub-problem solver,
    so the difference between this and
    :func:`maximise_efficiency_by_subdivision` is the master and nothing else.

    Args:
        n_subdivisions: The variables to subdivide, and how finely.
        bounds: The design box.
        initial: Incumbent, tried in every box before anything is drawn.
        samples: Crank angles per revolution.
        spec: Fixed engine data.
        targets: Constraint right-hand sides.
        seed: Random seed of the box probes.
        sub_problem_max_iter: Budget of the solver inside a box.
        algorithm: The solver inside a box.
        **start_settings: Forwarded to :func:`find_box_start`.

    Returns:
        One outcome per box, in the order of the Cartesian product.
    """
    from itertools import product

    from .reference import REFINED_DESIGN
    from .scenarios import _best_design, build_scenario, is_feasible

    incumbent = REFINED_DESIGN if initial is None else initial
    index_of = {name: index for index, name in enumerate(VARIABLE_NAMES)}
    names = tuple(n_subdivisions)
    outcomes: list[BoxOutcome] = []
    for combination in product(*(range(n_subdivisions[name]) for name in names)):
        lower = np.array(bounds.lower, dtype=float)
        upper = np.array(bounds.upper, dtype=float)
        for name, position in zip(names, combination, strict=True):
            index = index_of[name]
            width = (upper[index] - lower[index]) / n_subdivisions[name]
            edge = float(bounds.lower[index])
            lower[index] = edge + position * width
            upper[index] = edge + (position + 1) * width

        start = find_box_start(
            lower,
            upper,
            incumbent=incumbent,
            seed=seed,
            samples=samples,
            spec=spec,
            targets=targets,
            **start_settings,
        )
        box = Bounds(lower=lower, upper=upper)
        scenario = build_scenario(
            "neg_efficiency",
            bounds=box,
            initial=start.design,
            samples=samples,
            spec=spec,
            targets=targets,
            relax_equalities=True,
        )
        try:
            scenario.execute(algo_name=algorithm, max_iter=sub_problem_max_iter)
        except Exception:  # a box whose solver fails is data, not an error
            outcomes.append(
                BoxOutcome(combination, start, None, float("nan"), start.evaluations)
            )
            continue

        found = _best_design(scenario, objective="neg_efficiency", fallback=start.design)
        analysis = analyse(found, samples=samples, spec=spec, targets=targets)
        solved = is_feasible(analysis, targets)
        outcomes.append(
            BoxOutcome(
                index=combination,
                start=start,
                design=found if solved else None,
                efficiency=(
                    float(analysis.metrics.efficiency) if solved else float("nan")
                ),
                evaluations=start.evaluations
                + len(scenario.formulation.optimization_problem.database),
            )
        )
    return outcomes


def format_enumeration(outcomes: Sequence[BoxOutcome], title: str = "every box") -> str:
    """Render an enumeration as an aligned table."""
    lines = [title, "=" * 64, ""]
    lines.append(
        f"  {'box':<14}{'start':>13}{'margin':>10}{'efficiency':>13}{'evals':>9}"
    )
    for outcome in outcomes:
        label = ",".join(str(i) for i in outcome.index)
        efficiency = "-" if outcome.design is None else f"{outcome.efficiency:.4f}"
        lines.append(
            f"  {label:<14}{outcome.start.stage:>13}{outcome.start.margin:>10.3f}"
            f"{efficiency:>13}{outcome.evaluations:>9}"
        )
    solved = [o for o in outcomes if o.solved]
    lines.append("")
    lines.append(f"  {len(solved)} of {len(outcomes)} boxes yielded a design")
    if solved:
        best = max(solved, key=lambda o: o.efficiency)
        lines.append(f"  best {best.efficiency:.4f} in box {best.index}")
        lines.append(
            f"  total {sum(o.evaluations for o in outcomes)} analyses"
        )
    return "\n".join(lines)


def _decode_one_hot(
    key: tuple[float, ...], n_subdivisions: Mapping[str, int]
) -> tuple[int, ...]:
    """Turn the master's concatenated one-hot back into a box index.

    The master optimises over one binary per subdivision per variable, so the
    key of a box is a flat vector of those binaries; a reader wants ``(1, 0, 2)``.
    """
    index: list[int] = []
    position = 0
    for name in n_subdivisions:
        width = n_subdivisions[name]
        block = key[position : position + width]
        position += width
        chosen = [i for i, value in enumerate(block) if value > 0.5]
        index.append(chosen[0] if chosen else -1)
    return tuple(index)


def _designs_in(scenario: Any) -> list[Design]:
    """Every design a scenario's geometric discipline has evaluated."""
    from .disciplines import ExlinkDiscipline

    designs: list[Design] = []
    seen: set[bytes] = set()
    for discipline in _walk(scenario):
        if not isinstance(discipline, ExlinkDiscipline):
            continue
        for entry in discipline.cache:
            try:
                design = Design.from_mapping(dict(entry.inputs))
            except (KeyError, ValueError):
                continue
            key = design.to_array().tobytes()
            if key not in seen:
                seen.add(key)
                designs.append(design)
    return designs


def _walk(node: Any) -> Sequence[Any]:
    """Every discipline reachable from a scenario, chains included."""
    found: list[Any] = []
    stack = list(getattr(node, "disciplines", ()))
    while stack:
        item = stack.pop()
        found.append(item)
        stack.extend(getattr(item, "disciplines", ()))
    return found


def format_subdivision(outcome: SubdivisionOutcome, title: str = "box subdivision") -> str:
    """Render a run as the aligned block the other scenarios print."""
    lines = [title, "=" * 56, ""]
    total = outcome.evaluations + outcome.restoration_evaluations
    lines.append(f"  {'boxes opened':<24}{outcome.boxes:>12}")
    lines.append(f"  {'boxes with a design':<24}{outcome.feasible_boxes:>12}")
    lines.append(f"  {'analyses, box starts':<24}{outcome.restoration_evaluations:>12}")
    lines.append(f"  {'analyses, sub-problems':<24}{outcome.evaluations:>12}")
    lines.append(f"  {'analyses, total':<24}{total:>12}")
    lines.append(f"  {'seconds':<24}{outcome.seconds:>12.1f}")
    if outcome.design is None:
        lines.append("")
        lines.append("  no feasible design found")
        return "\n".join(lines)
    lines.append(f"  {'efficiency':<24}{outcome.efficiency:>12.4f}")
    if outcome.per_box:
        lines.append("")
        lines.append("  box                 efficiency")
        for index, efficiency in sorted(outcome.per_box):
            label = ",".join(str(i) for i in index)
            value = "-" if np.isnan(efficiency) else f"{efficiency:.4f}"
            lines.append(f"    {label:<18}{value:>12}")
    lines.append("")
    lines.append("  design")
    for name, value in zip(VARIABLE_NAMES, outcome.design.to_array(), strict=True):
        lines.append(f"    {name:<10}{value:>12.3f}")
    return "\n".join(lines)


def _mda_ignoring_its_own_couplings(mda: Any) -> Any:
    """Stop an MDA inside a chain from being differentiated w.r.t. its couplings.

    ``MDOChain`` treats every input no earlier discipline produces as an input
    of the chain, and an MDA both consumes and produces its couplings, so
    ``diameters`` is one.  The chain then asks the MDA to differentiate with
    respect to it, and the Jacobian assembly refuses outright::

        ValueError: Variable diameters is both a coupling and a design variable

    Dropping it is not a workaround, it is the right derivative.  The coupling
    enters the MDA as an initial guess and leaves it converged, and a converged
    fixed point does not depend on where the iteration started, so the
    derivative with respect to it is zero.  Under ``MDF`` the formulation knows
    that and never asks; inside a chain nothing does, so the MDA has to say so
    itself.

    Args:
        mda: The MDA, which is modified in place.

    Returns:
        The same MDA, differentiated with respect to design variables only.
    """
    couplings = frozenset(mda.coupling_structure.all_couplings)

    # The MDA's class is only known at run time -- it is whichever one
    # ``create_mda`` built -- so the base is dynamic by construction.
    class CoupledBlock(type(mda)):  # type: ignore[misc]
        """The MDA as one block of a chain, with its couplings internal."""

        def add_differentiated_inputs(self, input_names: Sequence[str] = ()) -> None:
            super().add_differentiated_inputs([
                name for name in input_names if name not in couplings
            ])

    mda.__class__ = CoupledBlock
    return mda


def build_subdivided_range_scenario(
    n_subdivisions: Mapping[str, int],
    bounds: Bounds = GLOBAL_BOUNDS,
    initial: Design | None = None,
    speed_rpm: float | None = None,
    samples: int | None = None,
    vehicle: Any = None,
    module: float | None = None,
    teeth: int | None = None,
    spec: EngineSpec = DEFAULT_SPEC,
    targets: DesignTargets = DEFAULT_TARGETS,
    material: Any = None,
    safety: Any = None,
    settings: Any = None,
    seed: int = 0,
    restoring: bool = True,
    main_level: bool = True,
    mda_name: str | None = None,
    mda_settings: Mapping[str, Any] | None = None,
    **start_settings: Any,
) -> Any:
    """Pose §4.7's range problem as a box subdivision.

    The headline problem rather than the geometric one: maximise kilometres per
    litre over the ten linkage variables that remain once the gear pair pins
    ``I``, under the geometric constraints, the coupled ones and the two the
    vehicle adds.

    The MDA is kept, and keeping it is the whole difficulty
    ------------------------------------------------------
    ``MDF`` is what makes every point the optimizer sees a self-consistent
    engine: the sizing/inertia fixed point between
    :class:`~exlink.disciplines.DynamicsDiscipline` and
    :class:`~exlink.disciplines.StructureDiscipline` is converged before range
    is computed.  A box subdivision, though, gives its sub-problem the
    ``DisciplinaryOpt`` formulation -- the master owns the formulation, so the
    sub-problem cannot also be an MDF scenario -- and handing it the five
    disciplines in a chain would evaluate them once each, in order, and call
    the result an engine.  It is not one: the diameters that leave
    ``StructureDiscipline`` never come back to ``DynamicsDiscipline``.

    So the MDA is built explicitly and placed in the chain as one discipline.
    The coupling is then converged inside a single chain step, which is what
    MDF does, and the sub-problem sees exactly the function §4.7 optimises.

    What it costs
    -------------
    One point of this problem is 10.4 s against 0.29 ms for the geometric one,
    measured over a forty-iteration solve rather than a handful of steps -- a
    short run reports about half that, because the MDA is warm-started and has
    less to converge near its own last answer.  Between the two problems that
    is four and a half orders of magnitude, and the box subdivision spends its
    budget on boxes rather than on iterations.  The box starts stay geometric --
    :func:`find_box_start` uses the cheap analysis -- because the gate that
    defeats a cold start is geometric: everything downstream of
    :class:`~exlink.disciplines.ExlinkDiscipline` is only reached by a design
    whose linkage closes.

    Args:
        n_subdivisions: The variables to subdivide, and how finely.  ``I`` is
            not among them; the gear pair pins it, as in
            :func:`~exlink.scenarios.build_range_scenario`.
        bounds: The design box.
        initial: Starting design, and the incumbent tried first in every box;
            :data:`~exlink.reference.COUPLED_DESIGN` if omitted.  This is
            **not** the default of
            :func:`~exlink.scenarios.build_range_scenario`, which starts from
            the published design, and the difference is deliberate: here the
            incumbent is also the warm start of every box, and the published
            design does not run at all -- it scores 0 km/L, because friction
            exceeds its indicated work -- so as a warm start it is worth
            nothing and every box pays a full restoration instead.
        speed_rpm: Crankshaft speed [rev/min].
        samples: Crank angles per revolution.
        vehicle: The car.
        module: Gear module [mm], derived from ``initial`` if omitted.
        teeth: Teeth on the small gear, derived likewise.
        spec: Fixed engine data.
        targets: Constraint right-hand sides.
        material: The material.
        safety: The design factors.
        settings: Settings of the subdivision.
        seed: Random seed of the box probes.
        restoring: Whether to restore a startable point in each box.
        main_level: Whether an infeasible box produces a feasibility cut.
        mda_name: Inner MDA.
        mda_settings: Overrides for the MDA settings.
        **start_settings: Forwarded to :func:`find_box_start`.

    Returns:
        A ``BoxSubdivisionScenario`` ready to execute.

    Raises:
        ValueError: When ``I`` is named among the variables to subdivide.
    """
    from gemseo import create_mda
    from gemseo_box_subdivision import BoxSubdivisionScenario, SweptBoxSubdivisionSettings

    from .disciplines import (
        DynamicsDiscipline,
        ExlinkDiscipline,
        RangeDiscipline,
        StructureDiscipline,
    )
    from .gears import lattice_inter_axle, size_pair, tooth_count
    from .materials import DEFAULT_MATERIAL, DEFAULT_SAFETY
    from .reference import COUPLED_DESIGN
    from .scenarios import (
        COUPLED_INEQUALITY_OUTPUTS,
        COUPLED_SAMPLES,
        DEFAULT_MDA,
        DEFAULT_MDA_SETTINGS,
        DEFAULT_SPEED_RPM,
        RANGE_INEQUALITY_OUTPUTS,
        BearingMarginDiscipline,
    )

    if "I" in n_subdivisions:
        msg = (
            "``I`` is pinned by the gear pair and is not a design variable of "
            "the range problem, so it cannot be subdivided."
        )
        raise ValueError(msg)

    speed = DEFAULT_SPEED_RPM if speed_rpm is None else speed_rpm
    crank_samples = COUPLED_SAMPLES if samples is None else samples
    the_material = DEFAULT_MATERIAL if material is None else material
    the_safety = DEFAULT_SAFETY if safety is None else safety

    start = COUPLED_DESIGN if initial is None else initial
    if module is None:
        module = size_pair(start.I, 1000.0).module
    if teeth is None:
        teeth = tooth_count(start.I, module)
    inter_axle = float(lattice_inter_axle(module, teeth))
    start = start.replace(I=inter_axle)

    geometry = ExlinkDiscipline(samples=crank_samples, spec=spec, targets=targets)
    geometry.set_cache(geometry.CacheType.MEMORY_FULL)
    mda = _mda_ignoring_its_own_couplings(create_mda(
        DEFAULT_MDA if mda_name is None else mda_name,
        [
            DynamicsDiscipline(
                speed_rpm=speed,
                samples=crank_samples,
                material=the_material,
                safety=the_safety,
                spec=spec,
            ),
            StructureDiscipline(
                samples=crank_samples,
                material=the_material,
                safety=the_safety,
                spec=spec,
            ),
        ],
        **{**DEFAULT_MDA_SETTINGS, **(mda_settings or {})},
    ))
    sides, side_names = _relaxed_equality_disciplines()
    negated, negated_names = _negated_disciplines(RANGE_INEQUALITY_OUTPUTS)
    disciplines: list[Discipline] = [
        geometry,
        *sides,
        mda,
        BearingMarginDiscipline(limit=targets.max_bearing_load),
        RangeDiscipline(
            speed_rpm=speed,
            samples=crank_samples,
            vehicle=vehicle,
            module=module,
            teeth=teeth,
            material=the_material,
            safety=the_safety,
            spec=spec,
        ),
        *negated,
    ]
    # ``I`` leaves the design space, so every discipline has to be told the
    # value the gear pair pinned it at; they would otherwise fall back to the
    # reference and quietly optimise a different engine.
    for discipline in disciplines:
        if "I" in discipline.input_grammar:
            discipline.default_input_data = {
                **discipline.default_input_data,
                "I": np.array([inter_axle]),
            }

    design_space = build_design_space(bounds, start, fixed=("I",))
    adapter = None
    if restoring:
        from gemseo_box_subdivision.subdivisions.box import BoxSubdivision

        adapter = create_restoring_adapter_class(
            BoxSubdivision.from_design_space(
                design_space, dict(n_subdivisions), tuple(n_subdivisions)
            ),
            bounds=bounds,
            incumbent=start,
            fixed={"I": inter_axle},
            seed=seed,
            samples=crank_samples,
            spec=spec,
            targets=targets,
            **start_settings,
        )

    scenario = BoxSubdivisionScenario(
        disciplines,
        "neg_range",
        design_space,
        n_subdivisions=dict(n_subdivisions),
        settings=SweptBoxSubdivisionSettings() if settings is None else settings,
        scenario_adapter_cls=adapter,
    )
    # Two sign conventions meet here, exactly as in ``build_range_scenario``:
    # the coupled margins are violations and the range margins are margins.
    # Unlike there, the second pair cannot be attached with ``positive=True``
    # -- see :func:`_negated_disciplines` -- so they arrive already negated and
    # every constraint of this problem is attached the same way.
    for name in (
        *INEQUALITY_OUTPUTS,
        *side_names,
        *COUPLED_INEQUALITY_OUTPUTS,
        *negated_names,
    ):
        scenario.formulation.add_constraint(
            name, constraint_type="ineq", main_level=main_level
        )
    return scenario
