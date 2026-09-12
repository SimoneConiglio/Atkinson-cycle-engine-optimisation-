r"""The five-bar's multidisciplinary optimization, posed as the EX-link's is.

§5.10 prices the mechanism §5.9's synthesis produced and finds it indefensible:
its connecting rod runs 60 degrees off the cylinder axis, piston friction eats
over half the indicated work, and its sizing fixed point runs away above about
1000 rpm.  None of that is a property of the *topology*.  The precision-point fit
was asked for a motion and charged for nothing else, so it spent the mechanism's
nine degrees of freedom on the motion alone.

This module asks the question that actually matters: **with the same nine
variables, and now charged for everything, what is this topology worth?**  It is
the same question §4.7 asks of the EX-link, in the same form --

.. math::

    \max_X \; R(X) \quad \text{s.t.} \quad g(X) \le 0, \;
    |STE - 74| \le \delta, \; |\varepsilon - 16| \le \delta, \;
    X \in [X_{lb}, X_{ub}]

-- and it is run through GEMSEO so that the two are solved by the same optimizer
under the same conventions.

Where the multidisciplinary analysis is
---------------------------------------
Inside the discipline, and converged before any output is returned.  The
five-bar's coupling is the same one the EX-link's MDA carries: member masses set
the inertia loads, the inertia loads set the diameters through yield, fatigue and
buckling, and the diameters set the masses.
:func:`exlink.fivebar.solve_sized` iterates that loop to a fixed point, so every
point the optimizer sees is a range computed on a self-consistent structure --
which is what :class:`exlink.disciplines.RangeDiscipline` achieves by sitting
downstream of an explicit MDA.

Carrying it monolithically rather than as a coupled graph costs the MDF/IDF
comparison of §4.5 and nothing else, and it buys the thing that matters here: the
fixed point is allowed to *fail to converge*, which for this mechanism is a real
and frequent outcome rather than a numerical hiccup.  A design whose masses run
away is infeasible, and saying so is more useful than converging a graph onto a
structure that cannot exist.

The objective is a ladder
------------------------
For the same reason :func:`exlink.synthesis.maximise_range_from_target` uses one:
the range is not computable everywhere.  A five-bar can fail to close, close but
not produce a four-stroke motion, produce one but not size.  Each failure is
scored worse than every success and carries whatever gradient is still available,
so the search is pushed back towards designs that run.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

import numpy as np

from .fivebar import (
    SYNTHESISED,
    AssemblyError,
    FiveBar,
    FiveBarPerformance,
    evaluate,
)
from .materials import FloatArray

#: Bounds on the nine variables, wide enough to contain a sensible engine.
BOUNDS: dict[str, tuple[float, float]] = {
    "q_1": (4.0, 60.0),
    "q_2": (4.0, 60.0),
    "I": (20.0, 160.0),
    "theta_r": (-180.0, 180.0),
    "theta_f": (-180.0, 180.0),
    "L_1": (20.0, 260.0),
    "L_2": (20.0, 260.0),
    "e": (40.0, 340.0),
    "x_1": (-200.0, 200.0),
}
"""The design box.

Deliberately generous on the rod lengths and the cylinder offset, because the
point of the exercise is to let the optimizer move the geometry that the
kinematic fit chose badly -- constraining it to a compact envelope in advance
would decide the answer instead of finding it.
"""

#: Outputs the discipline publishes.
OUTPUT_NAMES: tuple[str, ...] = (
    "neg_range",
    "range",
    "engine_mass",
    "brake_efficiency",
    "side_load_ratio",
    "stroke_error",
    "ratio_error",
    "cycle_error",
    "transmission_margin",
    "side_load_margin",
    "assembly_margin",
    "converged",
)

#: What an unanalysable design scores.
RANGE_UNAVAILABLE: float = 1.0e4
CYCLE_PENALTY: float = 1.0e5

#: Largest side-load ratio a design may carry, the EX-link's own constraint.
MIN_TRANSMISSION: float = 0.50
"""Smallest ``|sin mu|`` allowed at the floating pin -- a 30 degree transmission angle.

Alt's classical limit, and here it is load-bearing rather than decorative.  Asked
only to align the connecting rod, the restoration drove this quantity to its floor
and returned a mechanism with a side-load ratio of 0.045 and 2.34 MN on its main
journals: the singularity §5.1 warns about, reached in a different topology by a
different route.  Constraining the transmission angle is what makes the
side-load objective safe to pursue.
"""

#: Millimetres a five-bar must keep between itself and coming apart.
ASSEMBLY_MARGIN: float = 5.0
"""Millimetres a five-bar must keep between itself and coming apart.

Two ways it can: the two rods reaching each other along one line, where the
mechanism locks, and the connecting rod failing to reach the cylinder axis.  Both
are boundaries of the design space rather than penalties on it, and a search that
is not held off them will sit on one, because that is where a rod angle looks
best.
"""

#: Largest side-load ratio a design may carry, the EX-link's own constraint.
SIDE_LOAD_LIMIT: float = 0.30
"""``tan`` of the worst connecting-rod angle.

The EX-link reference sits at 0.018, two orders of magnitude inside this; a
conventional engine stays below about 0.3, and the synthesised five-bar is at
1.77.  Imposing it is the whole difference between the two exercises.
"""


@dataclass(frozen=True)
class FiveBarOutcome:
    """The result of a five-bar range optimization."""

    design: FiveBar
    performance: FiveBarPerformance | None
    evaluations: int
    message: str

    @property
    def feasible(self) -> bool:
        """Whether the final design meets every constraint it was solved under."""
        if self.performance is None or not self.performance.result.converged:
            return False
        result = self.performance.result
        from .topology import requirement_error

        return bool(
            requirement_error(result.kinematics.lam) <= 0.5 + 1.0e-6
            and self.performance.side_load_ratio <= SIDE_LOAD_LIMIT + 1.0e-6
            and result.kinematics.transmission >= MIN_TRANSMISSION - 1.0e-6
        )


def analyse(
    design: FiveBar,
    speed_rpm: float,
    samples: int = 240,
) -> dict[str, float]:
    """Score one five-bar, returning GEMSEO-convention outputs.

    Constraints follow the convention the rest of the package uses: an inequality
    is feasible when ``<= 0`` and an equality when ``== 0``.

    Every failure mode is caught and scored rather than raised, because an
    optimizer differencing a gradient will step into all of them.  The ladder runs
    downwards: a design that cannot even be assembled is worse than one whose
    motion is not a cycle, which is worse than one that will not size, which is
    worse than any design that runs.
    """
    from .cycle import PhaseError
    from .topology import requirement_error

    # Every constraint is reported *violated* in the failure branches, not zero.
    # Zero was tried and it is wrong in a way that quietly ruins the solve: a
    # design whose motion is not a cycle then satisfies every constraint, so
    # GEMSEO declares the point feasible and reports a "solution" with a side
    # load of 1.27 and no computable range at all.  Nothing unanalysable may read
    # as feasible.
    blank = {
        "range": 0.0,
        "engine_mass": 0.0,
        "brake_efficiency": 0.0,
        "side_load_ratio": 0.0,
        "stroke_error": 10.0,
        "ratio_error": 10.0,
        "cycle_error": 100.0,
        "transmission_margin": 10.0,
        "side_load_margin": 10.0,
        "assembly_margin": 1.0,
        "converged": 0.0,
    }
    try:
        performance = evaluate(design, speed_rpm, samples=samples)
    except AssemblyError:
        # The rods do not reach.  Nothing downstream has a value, and the only
        # gradient available is the geometric one, so hand back the worst rung.
        return {**blank, "neg_range": CYCLE_PENALTY}
    except PhaseError:
        return {**blank, "neg_range": 0.5 * CYCLE_PENALTY, "assembly_margin": -1.0}

    result = performance.result
    reach = performance.range_km_per_litre
    # The ladder, and every rung needs a gradient or the search stops on it.
    if not result.converged:
        # A mass spiral that never closed is not a structure.  Lighter designs are
        # the ones that close, so the mass is the tie-break that points the way
        # out -- scaled down so it cannot outweigh the rung itself.
        objective = RANGE_UNAVAILABLE + 0.1 * performance.engine_mass
    elif reach <= 0.0:
        # It sizes, but friction eats more than the gas delivers.  The *signed*
        # brake work says how close it is to delivering anything at all; the
        # clamped efficiency is flat here and would stall the search.
        objective = 0.5 * RANGE_UNAVAILABLE - 1.0e-3 * performance.brake_work
    else:
        objective = -reach

    return {
        "neg_range": float(objective),
        "range": float(reach),
        "engine_mass": float(performance.engine_mass),
        "brake_efficiency": float(performance.brake_efficiency),
        "side_load_ratio": float(performance.side_load_ratio),
        "stroke_error": float(result.expansion_stroke - 74.0),
        "ratio_error": float(result.compression_ratio - 16.0),
        "cycle_error": float(requirement_error(result.kinematics.lam)),
        "transmission_margin": float(MIN_TRANSMISSION - result.kinematics.transmission),
        "side_load_margin": float(performance.side_load_ratio - SIDE_LOAD_LIMIT),
        # A mass spiral that never closed is not a structure, so it is reported as
        # a constraint violation rather than left to the objective's ladder alone.
        "assembly_margin": -1.0 if result.converged else 1.0,
        "converged": 1.0 if result.converged else 0.0,
    }


def maximise_range(
    start: FiveBar | None = None,
    speed_rpm: float = 2000.0,
    samples: int = 240,
    max_iterations: int = 120,
    band: float = 0.5,
    bounds: dict[str, tuple[float, float]] | None = None,
    side_load_limit: float | None = SIDE_LOAD_LIMIT,
) -> FiveBarOutcome:
    """Maximise the five-bar's range under every constraint, through GEMSEO.

    Args:
        start: Initial design; the synthesised one by default, which is a
            feasible *kinematic* point and a hopeless mechanical one.
        speed_rpm: Crankshaft speed [rev/min].  The study's operating point is
            2000, which is where the synthesised geometry cannot be sized at all.
        samples: Crank angles per cycle.
        max_iterations: SLSQP iteration budget.
        band: Largest requirement error allowed [mm].
        side_load_limit: Cap on the connecting-rod side load, or ``None`` to
            leave it to the objective.  ``None`` is the more informative
            setting and not a relaxation: the range already charges the piston
            friction the side load causes, at the study's own coefficient, so
            the cap is a design rule standing in for a cost that is being
            computed anyway.  It is worth switching off here because for this
            topology the cap is not reachable -- holding the four-stroke cycle
            costs a side-load ratio of about 1.2 whatever else is done -- and a
            constraint nothing satisfies tells the optimizer nothing.
        bounds: The design box; :data:`BOUNDS` by default.

    Returns:
        The best design found and what it delivers.
    """
    from gemseo import create_design_space, create_scenario
    from gemseo.core.discipline import Discipline
    from gemseo.typing import StrKeyMapping

    from .fivebar import VARIABLE_NAMES

    start = SYNTHESISED if start is None else start
    box = BOUNDS if bounds is None else bounds
    calls = {"n": 0}

    class FiveBarDiscipline(Discipline):
        """The geared five-bar, with its sizing MDA converged inside."""

        auto_detect_grammar_files: ClassVar[bool] = False

        def __init__(self) -> None:
            super().__init__(name="FiveBar")
            self.input_grammar.update_from_names(VARIABLE_NAMES)
            self.output_grammar.update_from_names(OUTPUT_NAMES)
            self.default_input_data = start.to_mapping()

        def _run(self, input_data: StrKeyMapping) -> StrKeyMapping:
            calls["n"] += 1
            design = FiveBar.from_mapping(dict(input_data))
            scored = analyse(design, speed_rpm, samples=samples)
            return {name: np.array([value]) for name, value in scored.items()}

    space = create_design_space()
    values = start.to_array()
    for index, name in enumerate(VARIABLE_NAMES):
        low, high = box[name]
        space.add_variable(
            name,
            lower_bound=low,
            upper_bound=high,
            value=float(np.clip(values[index], low, high)),
        )

    discipline = FiveBarDiscipline()
    scenario = create_scenario(
        [discipline], "neg_range", space, formulation_name="DisciplinaryOpt"
    )
    # The two stroke requirements as relaxed inequalities, the way §4.7 relaxes
    # them: an equality gives SLSQP nowhere to stand when the range is a ladder.
    scenario.add_constraint("cycle_error", constraint_type="ineq", value=band)
    scenario.add_constraint("transmission_margin", constraint_type="ineq")
    if side_load_limit is not None:
        scenario.add_constraint(
            "side_load_ratio", constraint_type="ineq", value=side_load_limit
        )
    scenario.add_constraint("assembly_margin", constraint_type="ineq")

    # There is no analytic Jacobian to be had here: the outputs come through a
    # fixed point whose iteration count changes with the design, so nothing
    # differentiable is available to chain.  A central difference on the design
    # vector is the honest treatment, and the step is set in *millimetres*
    # because the variables are lengths and angles of very different magnitudes.
    scenario.set_differentiation_method("finite_differences", 1.0e-4)
    scenario.execute(algo_name="SLSQP", max_iter=max_iterations)

    problem = scenario.formulation.optimization_problem
    best = problem.solution.x_opt if problem.solution is not None else values
    design = FiveBar.from_array(
        np.asarray(best, dtype=float), branch=start.branch, piston_branch=start.piston_branch
    )
    performance: FiveBarPerformance | None
    try:
        performance = evaluate(design, speed_rpm, samples=samples)
    except (AssemblyError, Exception):
        performance = None
    return FiveBarOutcome(
        design=design,
        performance=performance,
        evaluations=calls["n"],
        message=str(getattr(problem.solution, "message", "")),
    )


def restore_feasibility(
    start: FiveBar | None = None,
    samples: int = 360,
    band: float = 0.5,
    max_iterations: int = 300,
    bounds: dict[str, tuple[float, float]] | None = None,
) -> FiveBar:
    """Find a five-bar that meets the cycle *and* points its rod down the bore.

    A restoration phase, and the reason for it is that the two things being
    demanded live in different parts of the model.  The strokes, the compression
    ratio and the connecting-rod angle are **pure kinematics** -- no masses, no
    sizing, no fixed point -- while the range needs the whole multidisciplinary
    analysis converged.  Asking for both at once makes the optimizer pay MDA
    prices to discover geometric facts, and from a start with a side-load ratio of
    1.77 against a limit of 0.30 it never gets far enough to earn them: SLSQP
    stalled on the full problem in 210 evaluations without ever reaching a design
    that could be sized.

    So the geometry is fixed first, on kinematics alone, at a few milliseconds an
    evaluation.  What comes out is a mechanism that delivers the specified motion
    with a rod that runs down the cylinder, which is the point
    :func:`maximise_range` can start from.

    Args:
        start: Initial design; the synthesised one by default.
        samples: Crank angles per cycle.  Not a free choice: the assembly and
            transmission margins are read off this grid, so a coarse one can miss
            the angle where the linkage is worst and return a design that comes
            apart when checked more finely.  At 240 angles that happens; 360 is
            the floor used here.  A guarantee would need interval arithmetic over
            the whole revolution rather than a sample of it.
        band: Largest requirement error allowed [mm], as
            :func:`exlink.topology.requirement_error` measures it.
        max_iterations: SLSQP iteration budget.
        bounds: The design box; :data:`BOUNDS` by default.

    Returns:
        The most rod-aligned design found that still meets the cycle.
    """
    from scipy.optimize import minimize

    from .cycle import PhaseError
    from .fivebar import VARIABLE_NAMES, solve
    from .topology import requirement_error

    start = SYNTHESISED if start is None else start
    box = BOUNDS if bounds is None else bounds
    lower = np.array([box[name][0] for name in VARIABLE_NAMES])
    upper = np.array([box[name][1] for name in VARIABLE_NAMES])

    def measure(vector: FloatArray) -> tuple[float, ...] | None:
        design = FiveBar.from_array(
            vector, branch=start.branch, piston_branch=start.piston_branch
        )
        try:
            motion = solve(design, samples=samples)
            error = requirement_error(motion.lam)
        except (AssemblyError, PhaseError):
            return None
        return (
            motion.side_load_ratio,
            error,
            motion.closure,
            motion.reach,
            motion.transmission,
        )

    cache: dict[bytes, tuple[float, ...] | None] = {}

    def cached(vector: FloatArray) -> tuple[float, ...] | None:
        key = np.ascontiguousarray(vector, dtype=float).tobytes()
        if key not in cache:
            if len(cache) > 8 * len(VARIABLE_NAMES):
                cache.pop(next(iter(cache)))
            cache[key] = measure(vector)
        return cache[key]

    def objective(vector: FloatArray) -> float:
        got = cached(vector)
        return 10.0 if got is None else got[0]

    def constraints(vector: FloatArray) -> FloatArray:
        got = cached(vector)
        if got is None:
            return np.full(4, -1.0)
        _side, error, closure, reach, transmission = got
        return np.array(
            [
                # The *whole* requirement, not two of its five parts.  Holding
                # only the expansion stroke and the compression ratio does not
                # specify a four-stroke cycle, and this restoration proved it:
                # asked for those two alone it dropped one top dead centre 47 mm
                # below the other -- a motion find_phases still accepts, whose
                # peak pressure then came out at 0.098 MPa instead of 5.5.
                band - error,
                # Clear of both ways a five-bar comes apart, and of the
                # transmission-angle singularity it is otherwise drawn towards.
                closure - ASSEMBLY_MARGIN,
                reach - ASSEMBLY_MARGIN,
                transmission - MIN_TRANSMISSION,
            ]
        )

    result = minimize(
        objective,
        np.clip(start.to_array(), lower, upper),
        method="SLSQP",
        bounds=list(zip(lower, upper, strict=True)),
        constraints=[{"type": "ineq", "fun": constraints}],
        options={"maxiter": max_iterations, "ftol": 1.0e-8, "eps": 1.0e-5},
    )
    return FiveBar.from_array(
        np.asarray(result.x, dtype=float),
        branch=start.branch,
        piston_branch=start.piston_branch,
    )
