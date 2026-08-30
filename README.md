# EX-link Atkinson-cycle engine optimization

A Python reconstruction of a 2015 student project at the Université de Technologie de
Compiègne: sizing an **extended-expansion (Atkinson) linkage** for a Shell Eco-marathon
single-cylinder engine, by optimizing eleven geometric design variables against three
competing objectives.

The original study was written in MATLAB. Here the physics is NumPy and the optimization
is driven by **[GEMSEO](https://gemseo.readthedocs.io)**, so the same problem can be handed
to a gradient-based solver, differential evolution, an augmented Lagrangian, or NSGA-II
without rewriting anything. The mechanism is animated with **matplotlib**.

---

## The problem

Honda's EX-link mechanism reaches top dead centre **twice per crankshaft revolution**, with
two *different* bottom dead centres. The short one sets the compression stroke, the long one
the expansion stroke — so the engine expands the burnt gas further than it compressed the
fresh charge, and recovers work an Otto cycle throws away with the exhaust.

Getting a linkage to do that, at a 74 mm expansion stroke and a compression ratio of 16,
while staying small and keeping the piston rod within 10° of vertical, is the design problem.

![the mechanism turning through one cycle](docs/figures/exlink.gif)

*Left: the linkage. Right: piston height, the p–V cycle, and crankshaft torque, with a
marker tracking the crank angle on each. Regenerate with `exlink animate`.*

![piston motion, cycle and torque](docs/figures/overview.png)

The piston reaches the **same** top dead centre twice (`g = 0.006 mm`) but two different
bottom dead centres — `STE = 74.00 mm` against `STC = 55.95 mm`. Torque is strongly positive
through expansion, negative through compression, and flat through intake and exhaust, where
the cylinder is at plenum pressure and the piston carries no gas load.

**Design variables** — `X = (a, c, I, x_b, y_b, x_1, e, q_1, q_2, θ_f, θ_r)ᵀ`

| | |
|---|---|
| `a` | swing rod `QA` |
| `c` | side `AD` of the trigonal link |
| `x_b`, `y_b` | position of `E` in the frame carried by `AD` |
| `e` | piston rod `EP` |
| `q_1`, `q_2` | cranks on the crankshaft and the eccentric shaft |
| `I` | distance between the two shafts |
| `x_1` | lateral offset of the cylinder axis |
| `θ_f`, `θ_r` | crank dephasing, and shaft-axis orientation |

**Formulation**

```
l_b ≤ X ≤ u_b
min  f(X)    = (−η,  H,  B)ᵀ                       efficiency, and the two envelope sizes
s.t. c(X)    = (mra−10, W−0.985, g−0.01, 10−d, γ−0.02)ᵀ ≤ 0
     c_eq(X) = (STE−74, ε−16)ᵀ = 0
```

Most of those constraints exist to make the problem **well posed**, not to express a
specification — and that is the interesting part of the original work:

- **`W ≤ 0.985`** — the two arccosine arguments of the closed-form inversion. If either
  reaches 1, the crankshaft cannot turn through a full revolution; it only rocks. Feeding
  `W` to the optimizer as a *number* rather than letting the analysis throw is what makes the
  global search converge.
- **`g ≤ 0.01 mm`** — the gap between the two top dead centres. Non-zero means burnt gas that
  cannot be scavenged.
- **four monotone phases** — a design whose piston goes up and down only once per revolution
  is a plain Otto engine, not an Atkinson one. Designs failing this are *penalised*
  (`η = 0`, `H = B = 1000`), never rejected.
- **`γ ≤ 0.02`** — piston side load, which drives friction.

---

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

## Use

```bash
exlink analyse                       # objectives and constraints of the reference design
exlink analyse --design published    # the design as tabulated in the 2015 report
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

---

## How the 2015 workflow maps onto GEMSEO

The report worked through four stages, each of which has a direct counterpart here:

| 2015 (MATLAB) | here (GEMSEO) |
|---|---|
| external penalty `F(X) = −η + r⁻²(c_eqᵀc_eq + ⟨c⟩ᵀ⟨c⟩)` | `PenalisedExlinkDiscipline` |
| conjugate gradient / simplex | `NELDER-MEAD`, `SLSQP`, `NLOPT_COBYLA` |
| genetic algorithm, 550 individuals | `DIFFERENTIAL_EVOLUTION` (`maximise_efficiency`) |
| MOEA seeded near a known optimum | `local_pareto` (`PYMOO_NSGA2` in a shrunk box) |
| moving limits on `H` and `B` | `sweep_moving_limits` |
| augmented Lagrangian, final polish | `Augmented_Lagrangian_order_0` (`refine`) |

### Why the multi-objective stage needs a relaxed problem

The report found its MOEA "still disappointing … even with big population I still got
solutions that were even worse than the ones I got with the gradient based methods." Sampling
2000 designs from a box around a good solution shows why:

| constraint | satisfied by |
|---|---|
| `d ≥ 10` | 96 % of analysable designs |
| `W ≤ 0.985` | 76 % |
| `mra ≤ 10` | 28 % |
| `\|ε − 16\|` | 8 % |
| `γ ≤ 0.02` | 6 % |
| `\|STE − 74\|` | 4 % |
| **`g ≤ 0.01`** | **0.1 %** |

`g` is nominally an inequality, but at 0.01 mm it is a **third equality in disguise**. With
three equalities in eleven variables the feasible set is a thin sheet and the joint hit rate
for a random population is of order 1e-7 — so NSGA-II returns a "front" of one point, because
one point is all it ever found.

`exlink` therefore runs the multi-objective stage on a *relaxed* problem (`moea_targets`),
treating it as a scouting device that maps the shape of the trade-off, and finishes with
`refine`, which drives `g`, `STE` and `ε` back onto their true targets:

```python
front = pareto_front(pop_size=80, max_gen=30)   # 16 designs, η 0.247 – 0.282
final = refine(front.design)                    # feasible, η = 0.280
```

`sweep_moving_limits` walks the same trade-off on the *unrelaxed* problem instead — every
step a warm-started local solve, so `g`, `STE` and `ε` are met exactly at each point:

```
   H limit   eta [%]    H [mm]   feasible
       236    28.11      236.0   True
       232    28.23      232.0   True
       228    27.74      228.0   True
```

It needs a real budget per step (the augmented Lagrangian, ~120 outer iterations): each
step must satisfy both equalities and all five inequalities while pushed against a limit it
did not previously meet, and a smaller budget returns "no feasible point" rather than a
trade-off.

---

## Does it reproduce the original?

**The model: yes, and it is independently verified.** The quasi-static force chain is pinned
by the principle of virtual work — in a massless, frictionless mechanism the instantaneous
power in must equal the power out, so `M_r(θ₁) = −P dλ/dθ₁` at every crank angle. The chain
reproduces that to machine precision (`tests/test_loads.py`). Every link also keeps its
length to 1e-9 mm over the revolution.

That check found a **sign slip in the report**: its printed inversion of the trigonal-link
moment equation gives the swing-rod load `A` with the wrong sign, which makes the torque
disagree with virtual work by a factor of about −4. The corrected sign is used here and
documented in `exlink/loads.py`.

**The published design vector: no — and it does not reproduce for the report's own code
either.** Re-analysed as printed, it gives `g = 8.5 mm` against the reported 0.0069. Since
`g` is exactly the quantity the optimizer drove to zero, the printed table cannot be the
design that produced the reported properties. It is printed to four significant figures, and
perturbing each variable by its rounding half-width moves the answers far too little to
explain the discrepancy. The design sits at `W = 0.982`, a hair from the singularity, where
the piston motion is extremely sensitive to the link lengths.

**The published *properties*: yes, closely.** Running the report's own final step — the
augmented Lagrangian, from the published table — lands on a fully feasible design that
matches the reported table almost line for line:

| | this reconstruction | report, 2015 |
|---|---|---|
| `η` efficiency | **27.87 %** | 27.76 % |
| `W` compatibility | 0.9811 | 0.9817 |
| `mra` rod angle | 8.68° | 7.55° |
| `g` TDC gap | 0.0060 mm | 0.0069 mm |
| `STE` expansion stroke | 74.000 mm | 73.98 mm |
| `ε` compression ratio | 16.000 | 15.98 |
| `γ` side load | 0.0181 | 0.02 |
| `B` width | 151.9 mm | 156.2 mm |
| `H` height | 238.6 mm | 256.7 mm |

That design ships as `exlink.reference.REFINED_DESIGN` and is what the CLI uses by default.

`d` (trigonal-link to cylinder clearance) is not compared: the report states the constraint
but not the geometric construction behind it, so it is reconstructed here — monotone in the
right direction and zero on contact, but not expected to match digit for digit.

---

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
```

## Tests

```bash
pytest -m "not slow"    # fast: physics, grammars, figures
pytest                  # everything, optimizers included
tox                     # across Python 3.10 – 3.12
```

## Reference

S. Coniglio, *Exlink Motor Mechanism Optimization*, Université de Technologie de Compiègne,
2015. Mechanism after Honda's
[EXlink](http://world.honda.com/powerproducts-technology/exlink/), with a crank added
between each shaft and its link, and the roles of the crankshaft and eccentric shaft
exchanged.

## License

MIT — see [LICENSE](LICENSE).
