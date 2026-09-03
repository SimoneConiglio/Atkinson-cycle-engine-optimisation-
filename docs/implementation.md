# Running the code

Everything reported in this documentation is produced by the package described
here and pinned by its test suite. {doc}`api` documents the modules
themselves; each opens with why its model is what it is, what was tried and
rejected, and where the approximations are.

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

The `moea` extra pulls in [`gemseo-pymoo`](https://gitlab.com/gemseo/dev/gemseo-pymoo) for
NSGA-II. Everything except `exlink pareto` works without it.

---

## Command line and library

```bash
exlink analyse                       # objectives and constraints of the reference design
exlink analyse --design published    # the historical baseline design (see Provenance)
exlink plot -o figures               # motion, p-V cycle, torque, mechanism
exlink animate -o figures/exlink.gif # animated mechanism + live cycle and torque
exlink animate --formulations \
  -o figures/formulations.gif        # each formulation's final design, side by side
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

---

## Module map

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
  animation.py     animated mechanism: alone, with a live dashboard, or one
                   panel per formulation on a common scale
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
  robustness.py    ISO 286 tolerance, propagated first-order and by sampling
  formulations.py  coupling strength from the MDA residuals; MDF vs IDF, and
                   the coupling recounted in a Fourier basis
  slidercrank.py   the second mechanism, on identical terms
  synthesis.py     prescribed-motion targets, the constrained fit, and range
                   under every constraint with the target as a fallback
```

## Reproducing the results

| result in {doc}`results` | produced by |
|---|---|
| the reference designs | `python examples/01_analyse_reference.py` |
| the augmented-Lagrangian polish | `python examples/03_optimize.py` |
| the efficiency/size trade-off | `python examples/04_pareto.py` |
| sizing, dynamics and the singularity result | `python examples/05_sizing_and_dynamics.py` |
| the mass budget, loss breakdown and mechanism comparison | `python examples/06_range.py` |
| the tolerance study, ISO grades and coupling curve | `python examples/07_robustness.py` |

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
and in {doc}`theory`; the document itself is not a citable reference and
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
{doc}`introduction`.

---

Next: [API reference](api.rst)
