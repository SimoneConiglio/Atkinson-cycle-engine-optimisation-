# 6. Results and discussion

Each result is stated, then supported, then discussed. Every number is computed
by the code of §4 and pinned by a test; none is transcribed. §7 collects the
limitations.

![each formulation's final design, turning on a common scale](figures/formulations.gif)

*What each objective converged to, from the same starting point, drawn at one
scale and one crank angle. Left to right: the geometric objective under an
augmented Lagrangian and under SLSQP with exact gradients, then minimum coupled
mass, then maximum range. The two geometric optima are long-limbed and stand
their cylinder high — §6.1 is about why — while the two that can see mass are
visibly shorter and squatter. Regenerate with `exlink animate --formulations`.*

## 6.1 The quasi-static optimum is the worst place to be

### Result

Maximising the lever-arm measure without inertia drives the design to
$W = 0.981$, a hair from the transmission-angle singularity, because that is
where the quasi-static lever arm is longest. Restoring inertia makes that the
worst available choice.

| swing rod | $W$ | $\eta$ | $H$ mm | moving mass | peak bearing |
|---|---|---|---|---|---|
| x1.00 | 0.9811 | 28.20 % | 238.5 | 1.039 kg | 12 629 N |
| x0.94 | 0.9670 | 27.79 % | 227.8 | 0.610 kg | 6 541 N |
| x0.88 | 0.9560 | 27.92 % | 218.0 | 0.498 kg | 6 647 N |
| x0.82 | 0.9488 | 28.56 % | 213.2 | 0.450 kg | 6 027 N |

Half the bearing load, a smaller envelope, less than half the mass — at equal or
better efficiency, at 1000 rpm.

### Why

The same proximity that lengthens the lever arm amplifies the accelerations:
joint $A$ sees 75 times the crank pin's. Since $m \sim (Ca)^3$ (§3.2), and every
inertia load scales as $\Omega^2$, structural mass grows as the **sixth power of
speed**:

| speed | moving mass | peak bearing load |
|---|---|---|
| 0 rpm | 0.25 kg | 7.7 kN |
| 1000 rpm | 1.03 kg | 12.5 kN |
| 1500 rpm | 8.43 kg | 245 kN |
| 2000 rpm | *no section is thick enough* | |

### Discussion

This is the clearest result in the study and the least fragile. It rests on the
equilibrium solve, which is verified against virtual work to machine precision
(§4.4), and the mechanism is understood rather than merely observed.

It also generalises. A well-conditioned slider-crank shows the *opposite* sign:
its peak main-bearing load **falls** with speed, 4735 N at rest to 2985 N at
4000 rpm, because the peak gas force lands near top dead centre where the
reciprocating inertia pulls the other way — classic inertia relief. Same physics,
opposite sign, and conditioning decides which.

## 6.2 A specified constraint cannot be manufactured

### Result

At IT8, the top-dead-centre gap bound of 0.01 mm has a process capability of
**0.11** against an industrial target of 1.33, and roughly two thirds of
nominally-conforming builds violate it.

The central finding is about *conditioning* — the mechanism sits near a singularity — and a
design chosen for nominal performance in a badly conditioned region is exactly what a
tolerance study exists to catch. Presenting a deterministic optimum without one would be
negligent.

Tolerances are ISO 286 IT grades, not invented numbers: `i = 0.45·D^(1/3) + 0.001·D` µm, with
IT8 at 25i for a machined member. Errors propagate two ways — **first order from the exact
Jacobians**, so a full assessment costs one extra Jacobian evaluation, and **Monte Carlo** to
check the linearisation, which is precisely what should be distrusted near a singularity.

```
  constraint               nominal   sigma_1st    sigma_MC     Cpk   violated
  expansion_stroke        -0.04992     0.03645     0.02011    0.83      11.0%
  compression_ratio       -0.04998    0.009347    0.005162    3.23       0.0%
  rod_angle                 -1.321    0.005555    0.005448   80.83       0.0%
  compatibility          -0.003854   4.029e-05    4.19e-05   30.66       0.0%
  tdc_gap                -0.004323     0.02173     0.01306    0.11      65.5%
  clearance                 -47.65     0.03584     0.03221  493.12       0.0%
  side_load              -0.001414   4.277e-05   4.098e-05   11.50       0.0%
```

**`g ≤ 0.01 mm` cannot be held.** The dimensions producing the top-dead-centre gap are held to
±0.011–0.031 mm at IT8, and combine to give `g` a standard deviation of 0.013 mm — larger than
the constraint band itself. Process capability is **0.11** against an industrial target of
1.33, and two thirds of nominally conforming builds violate it.

Scanning the IT ladder settles what to do about it. Holding `g` would need a tolerance unit
multiple of **1.25i**, below the tightest grade in the table. *No machining grade fixes it.*
This is a defect in the specification, not in any design that meets it, and the remedy is a
shim at assembly or a relaxed bound — not a better optimizer. Every other constraint is
comfortable.

First order overestimates σ by up to 80 % here, so it is **conservative**, not optimistic. Worth
stating: the opposite would make first-order robust design unusable in this region.

---

### Why

The dimensions producing $g$ are held to +/-0.011 to +/-0.031 mm and combine to
give it a standard deviation of 0.013 mm — larger than the constraint band
itself. Scanning the ISO ladder, holding it would need a tolerance unit multiple
of **1.25i**, below the tightest grade in the table. No machining grade fixes it.

Four independent routes agree:

| perturbation | effect on $g$ (bound 0.01 mm) |
|---|---|
| IT8 machining tolerance | $\sigma = 0.013$ mm |
| snapping $I$ 0.18 mm onto the gear lattice | $0.003 \to 0.058$ mm |
| minimum-norm equality projection | $0.0009 \to 0.0201$ mm |
| crank-angle resolution below 360 samples | 44 % error |

### The reliability statement

§3.8 computes a probability rather than a margin, over the seven constraints
whose uncertainty $\Sigma$ actually carries; it is evaluated *on* the solved
design rather than constrained during the search, so what follows is an audit
of that design and not a target it was held to (§3.10). For the coupled
reference design at IT8:

| | |
|---|---|
| system $P_f$, correlation kept | **0.645** |
| the same assuming independence | 0.563 |
| binding constraint | `tdc_gap` |
| bound needed for a $10^{-3}$ target | **0.054 mm** against 0.01 specified |

Keeping the correlation is not a formality and does not always reassure: the two
largest contributors are *anti*-correlated, so the system probability comes out
**above** the independent estimate.

Relaxing the gap to 0.054 mm moves the problem rather than removing it — the
stroke band becomes binding at $\beta = 0.68$, because that band is itself only
1.7 standard deviations wide and the design sits off-centre in it. Reaching
$10^{-3}$ *on that constraint* needs $\pm 0.09$ mm — but the **system**
probability at those bounds is $2\times10^{-2}$, and $10^{-3}$ for the system
needs $\pm 0.15$ mm. §6.4 tabulates the difference and solves the problem at
the wider bound.

### Most of this probability is avoidable, and free

The 0.645 is not the price of the requirements. It is the price of *ignoring
them while optimising*. Sampling 2500 designs about the reference and checking
the best by reliability against the full constraint set (§3.10) gives:

| | $\beta$ | $P_f$ | range |
|---|---|---|---|
| the reference design | $-0.373$ | 0.645 | 3338.3 km/L |
| best sampled, fully feasible | $+0.502$ | **0.308** | 3341.7 km/L |

The failure probability more than halves, every one of the twenty-five best
candidates is feasible, and the range does not fall — the best is 0.10 % higher
than the design it replaces. **The deterministic optimum is dominated on both
objectives at once.**

That is the standard argument for reliability-based design optimization,
measured on this problem rather than asserted: a deterministic optimizer
converges *onto* its active constraints, because nothing in the formulation
rewards standing off them, and a design sitting exactly on $g = 0$ fails half
the time. Backing off by a few hundredths of a millimetre is nearly free in
range and buys most of the reliability back.

Two things this does *not* say. It does not repair the specification — §6.2's
0.01 mm gap bound is still unattainable at any ISO grade, and $P_f = 0.308$ is
still far from a design target. And it is not a substitute for solving the
reliability-constrained problem: sampling found these points, whereas SLSQP
starting from the deterministic optimum could not move at all (§3.10).

### Discussion

The estimator is first-order and its weakest point is exactly where the finding
is strongest: $g$ is the most nonlinear constraint and FORM under-predicts its
failure probability, 0.42 against 0.54 sampled. That error is in the
conservative direction for the finding — the truth is worse than the estimate —
so the conclusion is not fragile even though the estimator is approximate.

This is a defect in the specification, not in any design meeting it, and it is
the single most useful result here for anyone who would build the engine.

## 6.3 Against an optimised conventional engine the advantage is firing frequency alone

### Result

Both mechanisms are held to the same compression ratio, the same clearance
volume and the same fuel per cycle, and both are sized by identical code. The
conventional engine is then *optimised over the variables it has* -- the rod
length, as the obliquity $r/l$, and the speed -- rather than given textbook
proportions, because comparing an optimised EX-link against a hand-set
slider-crank measures the optimization and not the topology.

| | slider-crank, hand-set | slider-crank, optimised | EX-link (Atkinson) |
|---|---|---|---|
| $r/l$, speed | 0.300, 2000 rpm | **0.195, 2151 rpm** | — , 1000 rpm |
| members / journals | 2 / 3 | 2 / 3 | 7 / 7 |
| indicated efficiency | 0.457 | 0.457 | 0.477 |
| mechanical efficiency | 0.740 | 0.787 | 0.853 |
| brake efficiency | 0.338 | 0.359 | 0.407 |
| engine mass | 19.3 kg | 16.9 kg | 12.2 kg |
| **range** | 2690 km/L | **2888 km/L** | **3338 km/L** |

Optimising the baseline is worth 7.4 % on its own, and all of it is mechanical:
the indicated efficiency is unchanged to four figures, because the compression
ratio and the charge are fixed and the thermodynamics cannot move. What moves
is the rod. A longer rod (obliquity 0.195 against 0.300) reduces the piston
side load, and with it the liner friction, which is the largest single loss in
the budget; the optimum is interior, because past about $r/l = 0.19$ the rod's
own mass and the taller engine start costing more than the friction it saves.

### The comparison, with and without the firing-frequency difference

In this model the EX-link completes all four strokes in one crankshaft
revolution -- that is what {py:func}`exlink.cycle.find_phases` requires of the
piston motion -- while a conventional four-stroke needs two. Per unit of work
the EX-link therefore accumulates half the journal rotation and half the piston
sliding. Removing that advantage means doubling its friction *and* halving its
power, since an engine that fires half as often makes half the power at the
same speed; both are done, and the result is re-scored through the vehicle
rather than compared as an efficiency, because halving the power changes the
operating point the burn-and-coast strategy can use.

| | vs hand-set baseline | vs optimised baseline |
|---|---|---|
| EX-link, as modelled (3338 km/L) | +24.1 % | **+15.6 %** |
| EX-link, as a four-stroke (2765 km/L) | +2.8 % | **−4.3 %** |

### Discussion

Against a hand-set baseline the reading was that extended expansion very nearly
fails to pay for itself. Against an optimised one it does not pay for itself:
**the sign changes**. An EX-link firing at conventional frequency is 4.3 %
*worse* than a conventional engine optimised under the same models, because the
two points of indicated efficiency that extended expansion buys do not cover
four extra journals, a gear train, and the mass of both.

So the +15.6 % that the mechanism does deliver is attributable to the
one-revolution cycle, not to extended expansion. That is a real advantage and
the reason to build the linkage — but it is a *cycle-rate* advantage, and it
would survive any other means of achieving the same firing frequency. It is
also the most fragile part of the result: it rests entirely on the phasing that
`find_phases` extracts from the piston motion, and an engine that did not
achieve it in practice would lose the whole margin and more.

### Both sides, optimised the same way

§6.4 finds that imposing constraints during the search rather than checking
them afterwards is worth 4.9 % to the EX-link. If that were a property of the
method it would apply to the baseline too, and this comparison would be unfair
until the baseline got the same treatment. It is not:

| baseline, optimised by | $r/l$ | speed | range | evaluations |
|---|---|---|---|---|
| Nelder-Mead, infeasible points rejected | 0.195 | 2151 rpm | 2887.7 km/L | 393 |
| SLSQP, both constraints imposed | 0.195 | 2141 rpm | **2887.7 km/L** | **100** |

The same optimum, to the tenth of a km/L, four times cheaper. The comparison
above is therefore unaffected by which method produced the baseline.

The reason is the one §6.4 turns on, seen from the other side: **imposing a
constraint helps only where the constraint binds**. The slider-crank has two
degrees of freedom and two constraints, and at its optimum neither is active —
the peak in $r/l$ is interior, set by rod mass against piston side load, not by
a bound. There is nothing for constraint-handling to buy. The EX-link's
optimum sits on several active constraints at once, which is why the same
change is worth 4.9 % there and nothing here. The 4.9 % is not a general
property of the method.

One conservatism worth naming before the caveats, because it runs the other
way. Neither cycle models gas exchange, so neither engine pays any pumping
work — and that is *not* even-handed. The EX-link opens its exhaust valve on a
charge expanded to $1.23\,p_0$ against the slider-crank's $1.70\,p_0$, so the
loss the model discards is about two and a half times larger for the
conventional engine (§7.2). The advantage reported here is therefore
conservative against the EX-link by an amount this model cannot quantify.

Two things this comparison is not. It is not a claim about extended-expansion
engines in general: {cite:t}`watanabe2006` measure a 4-point indicated-efficiency
gain on hardware, and this model reproduces 2 points, so if anything the
thermodynamic side here is conservative. And it is not a statement that the
optimised baseline is a good engine — it is the best engine *these models* can
build with two degrees of freedom, and its 0.195 obliquity is longer-rodded
than practice, which suggests the model prices rod length only through mass and
not through the packaging that would really limit it.

What makes the comparison worth anything at all is that both sides went through
the same sizing, friction, mass-budget and vehicle code, and both were
optimised; only the kinematics, the equilibrium system and the cycle differ.

:::{note}
The optimised baseline is found by Nelder-Mead restarts on a two-variable box
({py:func}`exlink.slidercrank.optimise_slidercrank`). A derivative-free method
is used here having been rejected in §2.6 for the EX-link, and the difference
is the geometry of the feasible set, not a change of mind: the slider-crank's
feasible set is a full-dimensional box in which every sampled point can be
evaluated, while the EX-link's is a manifold of measure zero on which sampling
finds nothing. The compression ratio is deliberately *not* a design variable —
this model has no knock limit, so left free it would rise without bound and the
comparison would measure a missing sub-model.
:::

## 6.4 The announced problem, solved

### Result

§3.10 states a problem: maximise range, hold every constraint, and constrain a
system probability of failure. Solving *that* problem rather than a
deterministic subset of it needs the bounds relaxed first, because §6.2 shows
the specification as written admits no reliable design at all. With the gap at
the 0.054 mm §6.2 derives and both bands at $\pm 0.15$:

| | start (`COUPLED_DESIGN`) | result |
|---|---|---|
| range | 3337.2 km/L | **3394.9 km/L** |
| system $P_f$ | $1.03\times10^{-3}$ | $1.34\times10^{-3}$ |
| system $\beta$ | 3.082 | **3.001**, against a target of 3.0 |
| engine mass | 12.17 kg | 12.94 kg |
| worst constraint | — | $-2.2\times10^{-7}$ |

**+1.7 %, every constraint satisfied, and the reliability target held on its
boundary.** 1352 evaluations, 61 minutes, the iteration cap reached rather than
a convergence test met.

### Why the gain is small, and why that is the result

Three solves of the same objective, differing only in what is imposed:

| | range | reliability | feasible as specified |
|---|---|---|---|
| constraints bind only at the end (`RANGE_DESIGN`) | 3388 km/L | audited: $P_f = 0.79$ | no, by $1.5\times10^{-4}$ |
| constraints imposed throughout | **3501 km/L** | audited, not constrained | no, by $2\times10^{-4}$ |
| the same, plus $\beta \ge 3$, bounds relaxed | 3394.9 km/L | **constrained: $P_f = 1.3\times10^{-3}$** | no — relaxed bounds |

Imposing the constraints *without* the reliability requirement reaches
3501 km/L — 3 % more than with it. It gets there by converging onto its active
constraints,
which is precisely what §6.2 shows destroys reliability: a design sitting on
$g = 0$ misses its requirements about half the time. The reliability
constraint forbids that, and the 3 % is what standing off the boundary costs.

The two figures are therefore not competing estimates of one quantity.
3501 km/L is the best *nominal* design; 3394.9 km/L is the best design that
also survives its own manufacturing scatter at $P_f \approx 10^{-3}$. Which
one is the answer depends on whether the engine is built once or built.

### The relaxation it needed, and a correction to §6.2

§6.2 reports that a $10^{-3}$ target needs the gap at 0.054 mm and the stroke
band at $\pm 0.09$ mm. That second figure is per-constraint reasoning, and the
*system* probability at those bounds is $2\times10^{-2}$:

| gap | band | system $P_f$ | $\beta$ |
|---|---|---|---|
| 0.010 | 0.05, as specified | 0.645 | $-0.373$ |
| 0.054 | 0.09, §6.2's figures | 0.021 | $+2.037$ |
| **0.054** | **0.15** | **0.00103** | **$+3.082$** |

The system is likelier to fail than any one of its constraints, so a bound
giving one constraint a $10^{-3}$ margin does not give the system one.

This design is consequently **not feasible against the specification as
written** — under the specified 0.01 mm gap and $\pm 0.05$ bands it is not
close. It is the best design under a specification relaxed until a reliable
design exists at all, which is the only form in which §6.2's question has an
answer.

### What it cost to make the problem solvable

Three defects stood between the formulation and a usable answer, and each was
invisible in the aggregates the runs reported:

| | symptom | cause |
|---|---|---|
| centre distance free while the gear pair was pinned | $I = 85.1$ against the 57.6 the pair realises | §3.7's catalogue relation makes $I$ an output; the search treated it as a variable |
| reliability index saturating | $\beta$ pinned at $-8.2095$ to fifteen digits in two runs | the orthant integrates to exactly 1 outside the band, so the index goes flat and its difference quotient is meaningless |
| evaluation cache too small | ~2× the necessary MDA calls | SLSQP differences the objective and the constraints over the same stencil in separate passes |

The second is the one worth carrying away. **A probability makes a poor
constraint wherever it saturates**: outside the band it carries no information
about how far outside, and a finite-difference probe straddling the band sees
an eleven-unit fall over a $10^{-5}$ step. The search is steered on
$\min_i \beta_i = \min_i (-g_i/\sigma_i)$, which is smooth through the band,
and the system probability is *reported* at the solution rather than assumed
from the target. That is weaker than constraining the system index — §3.8 says
why — and it is what makes the problem solvable at all.

This is the third place in the study where a flat region defeats a gradient
method, after §3.10's reliability stall and the objective ladder's middle rung.
Each time the repair is the same in shape: supply something that still varies
where the quantity of interest does not.

### The fallback earns its place only where it is needed

Two measurements say what the prescribed motion contributes, and they point in
opposite directions on purpose:

| start | evaluations that fell back | outcome |
|---|---|---|
| `COUPLED_DESIGN` (runs) | **0 of 836** | pure range maximisation |
| `REFINED_DESIGN` at 1250 rpm (does **not** run) | **26 of 103** | 0 km/L $\to$ 3336 km/L |

From a start that already runs, the ladder never fires and costs nothing. From
one where the engine will not run — friction exceeding indicated work, so km/L
simply does not exist — a quarter of the search is conducted on the target, and
the solve climbs out to a working engine. That is what a fallback should do,
and it is the case a constant penalty cannot handle: there is no gradient in a
constant.

The complementary measurement is the motion residual, which ends at **5.53 mm**
— far from the target. Once the range is computable everywhere the optimizer
abandons the prescribed motion entirely and pursues the objective. If the target
were pulling on the answer this number would be small; it is not, which is how
one tells a fallback from a constraint.

## 6.5 Supporting measurements

### How coupled the problem is

$\rho$ is the Gauss–Seidel contraction factor of §3.5.

| rpm | $\rho$ | sweeps | verdict |
|---|---|---|---|
| 0 | 0.0000 | 2 | weak |
| 500 | 0.1307 | 9 | moderate |
| 1000 | 0.6513 | 28 | strong |
| 1500 | 0.6819 | 42 | strong |

At rest $\rho = 0$ exactly, which is the sharpest available check that the
measure reflects the physics rather than the solver: with no inertia there is no
path from mass to load.

### What the derivatives buy

Minimising total moving mass at 1000 rpm, subject to every constraint and a
25 % efficiency floor:

| | COBYLA | SLSQP + differences | SLSQP + analytic |
|---|---|---|---|
| result | did not move | did not finish | 1.039 → 0.234 kg |
| cost | 120 evals | timed out | 40 evals, 148 s |

### The mixed-integer decomposition

Four gear candidates, 25 SLSQP iterations per sub-problem:

| | chosen pair | range | sub-solves | seconds |
|---|---|---|---|---|
| outer approximation | m=0.8, z=48 | 3366 km/L | **2** | 575 |
| exhaustive | m=1.0, z=39 | 3385 km/L | 4 | 1056 |

Half the sub-solves, 0.6 % short of the best lattice point. The convexification
options of §2.7 were enabled and measured to change nothing here: the master
terminates after two solves, which is less history than the adaptive correction
needs.

### Local optima

Manifold-projected restarts (§3.8) show the single-start efficiency optimum was
local: 30.91 % becomes 36.99 %. That better point is 443 mm tall against 320 and
sits on the $g$ bound — the single-objective efficiency problem is unbounded in
mechanism size, so a stronger search exploits that harder. On the range problem,
which is bounded, 0 of 6 restarts reached feasibility at an affordable budget.

### Reference designs

| design | $\eta$ | $H$ mm | $B$ mm | $W$ | $g$ mm | feasible |
|---|---|---|---|---|---|---|
| `PUBLISHED_DESIGN` | 35.62 % | 283 | 157 | 0.9892 | 8.5236 | no |
| `REFINED_DESIGN` | 27.80 % | 239 | 152 | 0.9811 | 0.0060 | yes |
| `GRADIENT_DESIGN` | 30.91 % | 320 | 159 | 0.9850 | 0.0095 | yes |
| `COUPLED_DESIGN` | 25.00 % | 198 | 131 | 0.9372 | 0.0070 | yes |
| `RANGE_DESIGN` | 25.46 % | 231 | 131 | 0.9319 | 0.0012 | no |

`COUPLED_DESIGN` is the design to compare against: it gives up five points of
$\eta$ to move off the singularity and gets a lighter, faster, longer-ranged
engine for it.

### The range optimization

SLSQP on `neg_range`, gear pair pinned, 1000 rpm, started from the coupled reference:

| | range | engine mass | `g` | strictly feasible |
|---|---|---|---|---|
| start (`COUPLED_DESIGN`) | 3338 km/L | 12.17 kg | 0.0067 mm | **yes** |
| best found (`RANGE_DESIGN`) | 3388 km/L | 12.47 kg | 0.0009 mm | no — see below |

The 1.5 % gain is modest, and the reason it is not simply banked is worth stating rather than
smoothing over.

`RANGE_DESIGN` satisfies every inequality, including the gap, at `g = 0.0009 mm`. It misses
the two *relaxed equalities* by 1.5 × 10⁻⁴ mm and 6.1 × 10⁻⁵ — SLSQP stopping within its own
convergence tolerance of the constraint it was handed. For scale, the tolerance study puts
the machining standard deviation of `STE` at 0.020 mm, **130 times larger**; no real part
would tell the two apart.

The obvious fix is to project it back onto the equality manifold, which
`project_onto_equalities` does exactly, by the minimum-norm Newton step from the analytic
Jacobians. That step is a few hundredths of a millimetre — and it moves `g` from 0.0009 to
0.0201 mm, twice its bound.

So the same wall appears from a fourth direction:

| perturbation | effect on `g` (bound: 0.01 mm) |
|---|---|
| IT8 machining tolerance on the members | `σ = 0.013 mm` |
| snapping `I` 0.18 mm onto the gear lattice | `0.003 → 0.058 mm` |
| minimum-norm equality projection | `0.0009 → 0.0201 mm` |
| tightest ISO grade that would hold it | 1.25i — off the ladder |

The honest reading is that **the specification is over-constrained**. The equality manifold
and the region `g ≤ 0.01` intersect in a sliver too thin for machining, gear selection or a
converged optimizer to land inside reliably. Treat `g` as an assembly adjustment — a shim on
the piston-rod length — or as a quantity to minimise, and `RANGE_DESIGN` is the answer at
3388 km/L. Under the specification as written, the best strictly feasible design is the one
we started from, and the 1.5 % is the price of a constraint that cannot be held.

That is not a result the optimizer could have delivered. It came out of the tolerance study,
and it is the single most useful thing in this repository for anyone who would actually build
the engine.

---

Next: [7. Conclusions](conclusions.md)
