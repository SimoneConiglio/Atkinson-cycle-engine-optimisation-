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
   H limit   eta [%]    H [mm]    B [mm]  feasible      <- examples/04_pareto.py
       236    28.113     236.0     151.7  True
       232    27.982     232.0     151.9  True
       228    27.209     228.0     152.2  False
       224    27.927     224.0     155.1  True
```

![efficiency against the two envelope dimensions](docs/figures/pareto.png)

Every limit is met exactly, and efficiency falls as the envelope shrinks — shortening `H`
from 236 to 224 mm costs about 0.19 points of efficiency and pushes `B` out from 151.7 to
155.1 mm. It is a sequence
of independent local solves rather than one global sweep, so the curve is not perfectly
monotone and a step can land infeasible — the table reports which. It also needs a real
budget per step (~120 augmented Lagrangian outer iterations, ~3 minutes for the four): each
step must satisfy both equalities and all five inequalities while pushed against a limit it
did not previously meet, and a smaller budget returns "no feasible point" instead.

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

## Gradients

The feasible set here is a **sliver**. At the reference design the two equality
constraints leave a band 0.1 mm wide on `STE` and 0.1 wide on `ε`, inside an 11-dimensional
box with sides of tens of millimetres, while `W` and `γ` sit within 0.4 % and 7 % of their
bounds. A derivative-free method cannot work in that — COBYLA returns its starting point
unchanged after 313 s and 120 evaluations, whatever the budget.

Everything the report derives is closed form, so the derivatives are too. `exlink.jacobian`
propagates them forward through the same chain the kinematics evaluates, and gets the
extremum-based metrics from the **envelope theorem**: for `max_θ f(X, θ)` attained at `θ*`,
the derivative is `∂f/∂X` evaluated there, because the term through the moving maximiser
carries `∂f/∂θ = 0`.

That matters more than "faster". These metrics are maxima over the crank angle, so the
sample attaining them *switches* as the design moves, and a difference quotient taken across
the switch is simply wrong — on `γ` at a 1e-4 mm step, 25 % wrong. GEMSEO's own
`check_jacobian` passes on all of them; `tests/test_jacobian.py` pins both the agreement and
that failure mode.

| | COBYLA | SLSQP + finite differences | SLSQP + exact gradients |
|---|---|---|---|
| moved from start | 0 mm | 32 mm | 38 mm |
| time | 313 s | 15 s | **4 s** |

Run from the published design, the gradient-based search reaches **η = 30.77 %**, feasible
at every resolution from 720 samples up, against the report's 27.76 % and the 27.87 % of the
augmented-Lagrangian reproduction:

| | η | `H` | `W` | `g` | feasible |
|---|---|---|---|---|---|
| report, 2015 | 27.76 % | 256.7 mm | 0.9817 | 0.0069 | — |
| augmented Lagrangian (`REFINED_DESIGN`) | 27.87 % | 238.5 mm | 0.9811 | 0.0060 | yes |
| SLSQP + exact gradients (`GRADIENT_DESIGN`) | **30.77 %** | 319.8 mm | 0.9850 | 0.0095 | yes |

Read that comparison carefully: nothing limits the envelope in this single-objective form, so
the extra efficiency is bought partly with size — the mechanism grows to `H = 320 mm`. What
it does show is that the report's final step stopped three points short of the efficiency
available to it, and that the method rather than the mechanism was the limit.

### Through the MDA, too

`exlink.dynamics_jacobian` differentiates the coupling itself, so GEMSEO assembles the
coupled derivative from each discipline's *local* Jacobian rather than differencing the whole
fixed point. Three ideas carry it:

- **The spectral operator is linear.** Accelerations are `Ω²·D²r`, so `da/dp = Ω²·D²(dr/dp)`
  — the same operator applied to the kinematic derivative arrays. No new mathematics.
- **The 18×18 solve** gives `dx/dp = A⁻¹(db/dp − dA/dp·x)`, reusing the factorisation the
  forward solve already made. `A`'s entries are moment arms, so `dA/dp` follows from the
  position derivatives.
- **The sizing bisection is never differentiated.** The diameter is defined implicitly by
  `U(d, N, M) = 1`, so the implicit function theorem gives `dd/dN = −(∂U/∂N)/(∂U/∂d)`.

Verified piece by piece against converged central differences — mass properties and
accelerations to round-off, the equilibrium solve to 5e-7, the internal loads to 4e-7, and
the sizing by directional derivatives — plus GEMSEO's `check_jacobian` on the dynamics
discipline.

The same finite-step trap appears here and is pinned in `tests/test_dynamics_jacobian.py`:
perturb the loads by 1e-4 and the crank angle attaining the fatigue extremum hops across a
near-flat minimum, putting the difference quotient tens of percent out. The implicit function
theorem doesn't care.

**What this buys.** Minimising total moving mass at 1000 rpm, from the augmented-Lagrangian
design, subject to every constraint and a 25 % efficiency floor:

| | COBYLA | SLSQP, FD through the MDA | SLSQP, analytic |
|---|---|---|---|
| result | did not move | did not finish | **1.039 → 0.234 kg** |
| cost | 120 evals, 313 s | timed out at 560 s | **40 evals, 148 s** |

| | refined (quasi-static optimum) | coupled optimum |
|---|---|---|
| total moving mass | 1.039 kg | **0.234 kg** |
| peak bearing load | 12 629 N | 7 504 N |
| `H` | 238.5 mm | 205.7 mm |
| `W` | 0.9811 | 0.9372 |
| `η` | 28.20 % | 25.00 % |

Four times lighter, a third off the bearing load, a smaller envelope, for three points of
efficiency — and it got there by moving off the transmission-angle singularity, which is
where the mass was going. Ships as `reference.COUPLED_DESIGN`, feasible at every crank-angle
resolution from 360 to 2880 samples.

**Still differenced:** `η`, `H`, `B` and the clearance in the analysis discipline. All are
smooth, none is tight, and `η` would additionally need the crank angle of top dead centre
differentiated, because the combustion pressure jump puts moving-boundary terms in its
integral.

---

## Sizing the parts, and the coupling it creates

The report stops before this on purpose — *"to have the masses of the pieces we have to know
their shape, so those passages are for another iteration"*. That iteration is here, and it
closes a loop the quasi-static study does not have:

```
    section diameters ──▶ member masses ──▶ inertia forces ──▶ internal loads
            ▲                                                        │
            └──────────── sizing against yield, ─────────────────────┘
                          fatigue and buckling
```

Neither half can go first. Under GEMSEO that is the **MDF** formulation:
`DynamicsDiscipline` and `StructureDiscipline` are strongly coupled — `diameters` one way,
the internal load histories the other — and an MDA converges them before the optimizer sees
a number.

```bash
exlink size --rpm 1000                    # size every part, inertia included
exlink size --rpm 3000 --plot fig.png     # ...and watch the loop run away
```

```python
from exlink import solve_for_design
sized = solve_for_design(PUBLISHED_DESIGN, speed_rpm=1000)
print(sized.total_mass_kg, sized.diameters, sized.feasible)

from exlink.scenarios import minimise_mass
best = minimise_mass(speed_rpm=1000, min_efficiency=0.25)
```

**Why the report's method cannot simply be extended.** Without inertia every rod is a
two-force member, which is exactly what lets the report eliminate unknowns one body at a
time. Give a rod mass and its end forces are no longer collinear, so nothing can be solved
before its neighbours. Gruebler gives the mechanism one degree of freedom, so the load
problem is statically determinate — **18 unknowns against 6 bodies × 3 equations** — and
`exlink.dynamics` assembles and solves that 18×18 system at every crank angle. Its
determinant is the same quantity condition (4a) protects: at a critical configuration the
matrix goes singular and the forces blow up.

Each member is then sized against **static yield**, **Goodman fatigue** (Marin-corrected
endurance limit, evaluated per extreme fibre) and **Euler buckling**, by bisection — needed
because the fatigue size factor `k_b` itself depends on the diameter being solved for.

### Two results worth the trouble

![bearing load and torque against crank angle, at three speeds](docs/figures/bearing_loads.png)

**Efficiency does not care about speed; the parts care enormously.** At constant speed the
mechanism returns to its starting state each revolution, so inertia does no net work and the
*mean* torque is provably unchanged — verified to 1e-6 across 0–3000 rpm. Only the peaks
move, and with the gas load switched off every reaction scales as exactly `Ω²`. Since mass
goes as the cube of the acceleration level, it goes as the **sixth power of engine speed**:

| speed | moving mass | peak bearing load |
|---|---|---|
| 0 rpm | 0.25 kg | 7.7 kN |
| 1000 rpm | 1.03 kg | 12.5 kN |
| 1500 rpm | 8.43 kg | 245 kN |
| 2000 rpm | *loop runs away — no section is thick enough* | |

![sized sections, mass and binding failure mode](docs/figures/sizing.png)

Fatigue binds on six of the seven members; the long slender piston rod goes to Euler
buckling instead.

**The quasi-static optimum turns out to be the wrong answer.** The report's design sits at
`W = 0.981`, a hair from the singularity, because that is where the quasi-static lever arm is
longest. But that same proximity is what amplifies accelerations — joint `A` sees 75× the
crank pin's. Backing off costs nothing and saves almost everything:

| swing rod | `W` | `η` | `H` [mm] | mass [kg] | peak bearing [N] |
|---|---|---|---|---|---|
| ×1.00 | 0.9811 | 28.20 % | 238.5 | 1.039 | 12 629 |
| ×0.94 | 0.9670 | 27.79 % | 227.8 | 0.610 | 6 541 |
| ×0.88 | 0.9560 | 27.92 % | 218.0 | 0.498 | 6 647 |
| ×0.82 | 0.9488 | 28.56 % | 213.2 | 0.450 | 6 027 |

Half the bearing load, *better* efficiency, a smaller envelope — and less than half the mass,
at 1000 rpm. At 1500 rpm the same comparison spans 8.43 kg against 2.03 kg. The near-singular
design was an artefact of leaving inertia out.

Both tables above are printed by `python examples/05_sizing_and_dynamics.py` (about 30 s).

The coupled problem therefore gains an objective the report could not express — structural
mass, since nothing in its formulation determines a cross-section — and three constraints
that only exist once loads are dynamic: `saturation_margin` (the loop ran away),
`slenderness_margin` (a "rod" thicker than a third of its length is not a beam), and
`bearing_margin`.

Full derivation, including the one idealisation used for the trigonal link and why the loop
converges at all, is in [docs/theory.md](docs/theory.md) §9.

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

  derivatives.py   exact spectral d/dtheta and d2/dtheta2 of periodic histories
  materials.py     material data, Marin factors, the Goodman line
  dynamics.py      inertia in the load path: the 18x18 equilibrium solve
  sizing.py        internal loads, then static / fatigue / buckling sizing
  coupled.py       the sizing <-> dynamics fixed point, and why it converges
  jacobian.py      exact d/dX of the analysis chain, by forward mode + envelope
  dynamics_jacobian.py   exact derivatives through the coupling itself
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
