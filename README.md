# EX-link Atkinson-cycle engine optimization

**[Read the study →](https://simoneconiglio.github.io/Atkinson-cycle-engine-optimisation-/)**

Multidisciplinary design of an extended-expansion (Atkinson) engine linkage for a Shell
Eco-marathon car, with **range** — distance on a given quantity of fuel — as the objective.

The formulation conventionally used for this mechanism maximises a lever-arm quality measure
subject to envelope bounds. It prices nothing, so it yields a Pareto front and never a design;
its central quantity is not an efficiency, because with no friction nothing is lost in it; and
it cannot see the parts, because nothing in it determines a cross-section. Range makes
efficiency, envelope size, torque ripple and structural mass commensurable at exchange rates
the physics fixes.

Built on **[GEMSEO](https://gemseo.readthedocs.io)**, with exact analytic derivatives through
the parts where finite differences are not merely inaccurate but wrong.

![the mechanism through one revolution](docs/figures/exlink.gif)

## Three results

**The quasi-static optimum is the worst place to be.** Maximising efficiency without inertia
drives the linkage to its transmission-angle singularity — exactly where the accelerations and
bearing loads are worst. That design has no feasible structure above 1000 rpm; one backed off
weighs half as much and goes further.

**A specified constraint cannot be manufactured.** The top-dead-centre gap is bounded at
0.01 mm and the dimensions producing it scatter by more than that. No ISO grade holds it, and
the reference design has a **64.5 %** probability of missing at least one requirement.

**The advantage is firing frequency, not extended expansion.** Against a conventional engine
*optimised* under identical code the linkage wins 15.6 %; remove its one-revolution cycle and
it loses by **4.3 %**.

## The study

Written as a paper, hosted on GitHub Pages, built from `docs/` by
[`.github/workflows/docs.yml`](.github/workflows/docs.yml).

| | |
|---|---|
| [1. Introduction](docs/introduction.md) | the problem and what is at stake |
| [2. State of the art](docs/state_of_the_art.md) | the methods available for each feature, and which apply |
| [3. Methodology](docs/methodology.md) | the formulation, and each method as a consequence of a stated limitation |
| [4. Implementation framework](docs/framework.md) | how the methodology maps onto code, and what is verified |
| [5. The use case](docs/use_case.md) | the mechanism, the specification, the set-up |
| [6. Results and discussion](docs/results.md) | each result stated, supported, discussed |
| [7. Conclusions](docs/conclusions.md) | limitations of the framework and possible improvements |
| [Running the code](docs/implementation.md) | install, CLI, module map, reproducing each result |
| [API reference](docs/api.rst) | every module, and the reasoning inside it |
| [Theory](docs/theory.md) | the derivations |

Build locally with `make docs`.

## Quick start

```bash
pip install -e ".[all]"        # or: make dev
pytest -m "not slow"
```

```python
from exlink import COUPLED_DESIGN, evaluate

outcome = evaluate(COUPLED_DESIGN, speed_rpm=1000.0)
print(outcome.km_per_litre)          # 3338.3
print(outcome.budget.kilograms())    # where the mass actually is
```

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

MIT.
