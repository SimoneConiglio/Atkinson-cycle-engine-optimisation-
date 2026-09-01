# EX-link Atkinson-cycle engine optimization

Multidisciplinary design of an **extended-expansion (Atkinson) linkage** for a Shell
Eco-marathon engine, optimised for the only thing the competition scores: **how far the car
gets on a given quantity of fuel**.

Eleven geometric design variables, a discrete gear choice, a strongly coupled
structure/dynamics analysis solved as an MDA, and a chain that carries all of it through to
kilometres per litre - so efficiency, envelope size, torque ripple and structural mass are
priced against each other by physics rather than by weights.

The physics is NumPy, with exact analytic derivatives through the parts where finite
differences are not merely inaccurate but wrong. The optimization is driven by
**[GEMSEO](https://gemseo.readthedocs.io)**. The mechanism is animated with **matplotlib**.

![the mechanism through one revolution](docs/figures/exlink.gif)

## Three findings

**The quasi-static optimum is the worst place to be.** Maximising efficiency without inertia
drives the linkage to its transmission-angle singularity - exactly where the accelerations,
the bearing loads and hence the structure are worst. Backing off costs nothing and halves the
engine.

**A specified constraint cannot be manufactured.** The top-dead-centre gap is bounded at
0.01 mm, and the dimensions that produce it scatter by more than that. Four independent
routes agree, and a reliability formulation puts a number on it: a **64.5 %** chance of
missing at least one requirement.

**The linkage's advantage is firing frequency, not extended expansion.** Against a
slider-crank sized by identical code it wins 24 %; remove the one-revolution cycle and the
advantage falls to **2.8 %**.

Each is pinned by a test, and each is explained in the documentation rather than here.

## Documentation

```bash
make docs          # or: sphinx-build docs docs/_build/html
```

| page | what is in it |
|---|---|
| [The design problem](docs/problem.md) | the mechanism, the specification, why the feasible set is a sliver |
| [Formulation](docs/formulation.md) | objectives, the discipline graph, the mass budget, MDF vs IDF |
| [Sizing and coupling](docs/coupling.md) | inertia in the load path, and the singularity result |
| [Exact derivatives](docs/gradients.md) | forward mode, the envelope theorem, derivatives through the MDA |
| [Discrete variables](docs/discrete.md) | the gear lattice, bi-level outer approximation, manufacturability |
| [Tolerance and reliability](docs/reliability.md) | ISO 286 propagation, and reliability in the formulation |
| [Local optima](docs/search.md) | why uniform multistart cannot work here, and what does |
| [Verification](docs/validation.md) | independent checks, and the slider-crank comparison |
| [Results](docs/results.md) | what the optimizations produced, and what they did not |
| [Theory](docs/theory.md) | the full derivation |


## Install

Any of these works; pick whichever fits your setup.

```bash
# uv (fastest)
uv venv && uv pip install -e ".[dev,moea]"

# pip
python -m venv .venv && . .venv/bin/activate
pip install -e ".[dev,moea]"

# conda
conda env create -f environment.yml && conda activate exlink
pip install -e ".[dev]"

# tox — creates its own envs and runs the suite
tox

# make
make dev && make test
```

Four optional extras, each needed only by the module that uses it:

| extra | pulls in | needed by |
|---|---|---|
| `moea` | `gemseo-pymoo` | NSGA-II fronts (`exlink pareto`) |
| `minlp` | `gemseo-bilevel-outer-approximation` | the mixed-integer gear choice (`exlink.minlp`) |
| `uq` | `gemseo-umdo` | the sampling reliability reference |
| `docs` | Sphinx and friends | building the documentation |

`pip install -e ".[all]"` takes everything.

## Use

```bash
exlink analyse                       # objectives and constraints of the reference design
exlink analyse --design published    # the historical baseline design (see Provenance)
exlink plot -o figures               # motion, p-V cycle, torque, mechanism
exlink animate -o figures/exlink.gif # animated mechanism + live cycle and torque
exlink refine --design published --save refined.json   # augmented Lagrangian
exlink optimize --save best.json     # differential evolution over the full box
exlink pareto --pop-size 200 --max-gen 60              # NSGA-II front
```

The commands chain through JSON design files:

```bash
exlink optimize --save best.json && exlink refine -d best.json --save final.json
exlink animate -d final.json -o final.gif
```

As a library:

```python
from exlink import analyse, PUBLISHED_DESIGN
from exlink.scenarios import refine, pareto_front, format_analysis

print(format_analysis(analyse(PUBLISHED_DESIGN)))

outcome = refine(PUBLISHED_DESIGN)          # augmented Lagrangian polish
print(outcome.feasible, outcome.analysis.metrics.efficiency)

front = pareto_front(pop_size=200, max_gen=60)   # NSGA-II
```

The vehicle-level chain, which is what the study actually optimises:

```python
from exlink import COUPLED_DESIGN, evaluate
from exlink.robustness import failure_probability, format_reliability

outcome = evaluate(COUPLED_DESIGN, speed_rpm=1000.0)
print(outcome.km_per_litre, outcome.engine_mass_kg)
print(outcome.budget.kilograms())          # where the mass actually is

print(format_reliability(failure_probability(COUPLED_DESIGN)))
```

The mixed-integer gear choice needs the `minlp` extra:

```python
from exlink.minlp import candidates_from_design, solve, format_result

points = candidates_from_design(COUPLED_DESIGN, speed_rpm=1000.0)
print(format_result(solve(COUPLED_DESIGN, candidates=points, speed_rpm=1000.0)))
```

## Layout

```
src/exlink/
  constants.py     engine spec, constraint targets, penalty values
  design.py        the 11-variable design vector and its box bounds
  kinematics.py    closed-form loop-closure inversion  →  λ(θ₁)
  cycle.py         four-phase detection + approximate Atkinson cycle
  loads.py         quasi-static force chain  →  crankshaft torque
  metrics.py       η, φ, H, B, mra, W, g, d, γ
  model.py         analyse(): the whole chain, penalising instead of raising
  disciplines.py   GEMSEO Discipline wrappers
  scenarios.py     design space, constraints, and the five workflows
  plots.py         motion / p-V / torque / mechanism / Pareto figures
  animation.py     animated mechanism, alone or with a live dashboard
  cli.py           the `exlink` command

  derivatives.py   exact spectral d/dtheta and d2/dtheta2 of periodic histories
  materials.py     material data, Marin factors, the Goodman line
  dynamics.py      inertia in the load path: the 18x18 equilibrium solve
  sizing.py        internal loads, then static / fatigue / buckling sizing
  coupled.py       the sizing <-> dynamics fixed point, and why it converges
  jacobian.py      exact d/dX of the analysis chain, by forward mode + envelope
  dynamics_jacobian.py   exact derivatives through the coupling itself

  friction.py      journal, ring and mesh losses -> real mechanical efficiency
  gears.py         the 2:1 pair, and the lattice the module imposes on I
  manufacturing.py preferred stock sizes, IT grades, minimum wall thicknesses
  mass_budget.py   the whole engine: crankcase from H x B, flywheel from ripple
  vehicle.py       burn-and-coast road load  ->  kilometres per litre
  performance.py   the full chain, from 11 dimensions to km/L
  robustness.py    ISO 286 tolerance, and reliability in the formulation
  formulations.py  coupling strength from the MDA residuals; MDF vs IDF
  slidercrank.py   the second mechanism, on identical terms
  minlp.py         the gear choice as a MINLP, by bi-level outer approximation
```

Every module opens with why its model is what it is, what was tried and rejected, and where
the approximations are; the [API reference](docs/api.rst) renders those.


## Tests

```bash
pytest -m "not slow"    # fast: physics, grammars, figures
pytest                  # everything, optimizers included
tox                     # across Python 3.10 – 3.12
```

Every claim above is pinned by a test, including the ones that would be embarrassing to get
wrong: that inertia does no net work over a cycle, that burn-and-coast conserves energy, that
the Otto cycle reproduces its closed-form efficiency, that `ρ = 0` at rest, that stock
rounding never shrinks a member, and that the range advantage collapses when the firing
frequency assumption is removed.


## Provenance

The mechanism and the design brief come from an unpublished student study by the author
(Université de Technologie de Compiègne, 2015), which set up the kinematics, the idealised
cycle, the quasi-static load chain and the efficiency measure, and solved the quasi-static
problem in MATLAB. Everything needed to read, run and check this repository is restated here
and in [docs/theory.md](docs/theory.md); the document itself is not a citable reference and
nothing here depends on it.

Two designs carry over from that study as **historical baselines**, and they are labelled as
such wherever they appear:

- `PUBLISHED_DESIGN` — the design vector tabulated there. Re-analysed it violates five
  constraints (notably `g = 8.5 mm` against a 0.01 mm bound), so it is used as a *starting
  point*, not as a result. See `exlink/reference.py` for why it cannot be the design that
  produced the properties reported alongside it.
- `REFINED_DESIGN` — what an augmented Lagrangian makes of it here: feasible, `η = 27.87 %`.
  The quasi-static reference point the rest of the study is measured against.

Everything else — the dynamic load analysis, the sizing disciplines, the coupled MDA, the
analytic derivatives, and the `GRADIENT_DESIGN` and `COUPLED_DESIGN` results — is new work in
this repository.

Mechanism topology after Honda's
[EXlink](https://global.honda/en/power/technology/exlink/), modified as described in
[The design problem](#the-design-problem).


## License

MIT — see [LICENSE](LICENSE).
