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

## 6.0 Every design in one place

Nine designs appear below. They differ in what was optimised, what was
imposed and what specification was applied, and quoting them apart is how a
document ends up contradicting itself — so they are collected here once and
referred to rather than restated.

| design | range | what it was | feasible as specified |
|---|---|---|---|
| `PUBLISHED_DESIGN` | — | the historical baseline | no, five constraints |
| `REFINED_DESIGN` | — | geometric objective, augmented Lagrangian | yes |
| `GRADIENT_DESIGN` | — | geometric objective, SLSQP | yes |
| `COUPLED_DESIGN` | 3338 km/L | minimum coupled mass; the strictly feasible reference | **yes** |
| `RANGE_DESIGN` | 3388 km/L | range, constraints bound at the end | no, by $1.5\times10^{-4}$ |
| range, constraints imposed | 3501 km/L | §3.10's second form, nominal only | no, by $2\times10^{-4}$ |
| **range + reliability, relaxed bounds** | **3395 km/L** | **§3.10's third form; $P_f = 1.3\times10^{-3}$** | no — relaxed spec |
| slider-crank, optimised | 2888 km/L | baseline, its own limits only | n/a — misses the EX-link's |
| slider-crank, same specification | 2467 km/L | baseline, held to the EX-link's limits | yes, on the cap |

Three advantage figures follow, and all are quoted in this document because they
answer different questions (§6.3 tabulates a fourth):

| question | figure |
|---|---|
| at matched power strokes per minute, 800–1400 | **+13 to +19 %** |
| against a conventional engine optimised as such | **+15.6 %** |
| against one held to the same specification | **+37.6 %** |

Figures are given to four significant digits once, here, and rounded
elsewhere.

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

## 6.2 Tolerance decides which of the stated bounds are real

### Result

The requirements of §5.2 are a mathematical specification: eight numbers written
down before any part existed. A tolerance study at IT8 says which of them this
mechanism can hold, and only two are in question — the top-dead-centre gap $g$
and the band each equality is relaxed into. Widening the gap from 0.01 to
0.1 mm and the bands from $\pm 0.05$ to $\pm 0.15$ takes the reference design
from a 0.645 probability of missing a requirement to $1.9\times10^{-5}$. The
physical price is **0.47 % of range**, all of it from the gap: a wider band
relaxes a constraint and cannot cost anything, whereas a dead-centre mismatch
of 0.1 mm is 2.7 % of the clearance volume and is felt by the cycle.

### Why

Tolerances are ISO 286 IT grades, not invented numbers: $i = 0.45 D^{1/3} +
0.001 D$ µm, with IT8 at $25i$ for a machined member. Errors propagate two ways
— **first order from the exact Jacobians**, so a full assessment costs one extra
Jacobian evaluation, and **Monte Carlo** to check the linearisation, which is
precisely what should be distrusted near a singularity.

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

Two rows are near their bounds and five are not. `tdc_gap` has a standard
deviation of 0.013 mm against a bound of 0.01 mm — the scatter is wider than the
requirement — and `expansion_stroke` has 0.036 mm against a half-band of
0.05 mm. The remaining five run from $C_{pk} = 3.2$ to 493 and do not enter the
discussion again.

First order overestimates $\sigma$ by up to 80 % here, so it is **conservative**,
not optimistic; worth stating, because the opposite would make first-order
robust design unusable in this region.

$g$ is the most sensitive quantity in the problem, and four independent
perturbations agree on its scale:

| perturbation | effect on $g$ |
|---|---|
| IT8 machining tolerance | $\sigma = 0.013$ mm |
| snapping $I$ 0.18 mm onto the gear lattice | $0.003 \to 0.058$ mm |
| minimum-norm equality projection | $0.0009 \to 0.0201$ mm |
| crank-angle resolution below 360 samples | 44 % error |

A bound of 0.01 mm lies below every one of them — below the machining scatter,
below the spacing of the gear catalogue, below the optimizer's own convergence,
and below the discretisation at which $g$ is computed. It is not a requirement
the rest of the model can resolve. A bound of 0.1 mm lies above all four.

### What widening the bounds costs

$g$ is the distance between the two top dead centres, and the cycle feels it
only through the volume trapped above the piston. The clearance volume is
3000 mm³, which over a 32 mm bore is 3.73 mm of head space:

| | mismatch | trapped volume | realised $\varepsilon$ | range |
|---|---|---|---|---|
| as specified | 0.010 mm | $+8.0$ mm³, $+0.27$ % | 15.96 | — |
| §6.4 relaxation | 0.054 mm | $+43.4$ mm³, $+1.45$ % | 15.79 | $-0.25$ % |
| accepted here | 0.100 mm | $+80.4$ mm³, $+2.68$ % | 15.61 | $-0.47$ % |

On one of its two revolutions the engine realises a compression ratio of 15.6
rather than 16.0, and the reference design loses 0.47 % of its range. That is
the entire consequence of the relaxation, computed rather than argued.

The bands on the equalities are the same kind of statement, in different units.
$\varepsilon = 16 \pm 0.15$ is $\pm 0.94$ % of the ratio, which is $\pm 0.035$ mm
of piston height — a shim under the cylinder head. $STE = 74 \pm 0.15$ mm is
$\pm 0.2$ % of the stroke.

### What widening them buys

§3.8 computes a probability rather than a margin, over the seven constraints
whose uncertainty $\Sigma$ actually carries. Evaluated on `COUPLED_DESIGN`:

| gap bound | band | system $P_f$ | $\beta$ | binding |
|---|---|---|---|---|
| 0.010 mm | $\pm 0.05$ | 0.645 | $-0.37$ | `tdc_gap` |
| 0.054 mm | $\pm 0.05$ | 0.251 | 0.67 | `stroke_lower` |
| 0.100 mm | $\pm 0.05$ | 0.250 | 0.67 | `stroke_lower` |
| 0.100 mm | $\pm 0.12$ | $1.0\times10^{-3}$ | 3.09 | `stroke_lower` |
| 0.100 mm | $\pm 0.15$ | $1.9\times10^{-5}$ | 4.12 | `stroke_lower` |

The second row is the useful one: **once the gap is at 0.054 mm it stops
binding, and no further widening of it changes anything.** The two bottom rows
differ from the third only in the band, and they span four orders of magnitude
of failure probability. What sets the reliability of this mechanism is not the
gap but how nearly the expansion stroke is required to equal 74 mm.

Keeping the correlation is not a formality and does not always reassure: at the
first row the two largest contributors are *anti*-correlated, so the system
probability, 0.645, comes out **above** the 0.563 an independence assumption
gives.

### Most of this probability is avoidable, and free

The 0.645 is not the price of the requirements. It is the price of *ignoring
them while optimising*. Sampling 2500 designs about the reference and checking
the best by reliability against the full constraint set (§3.10) gives:

| | $\beta$ | $P_f$ | range |
|---|---|---|---|
| the reference design | $-0.373$ | 0.645 | 3338 km/L |
| best sampled, fully feasible | $+0.502$ | **0.308** | 3342 km/L |

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

Sampling is how those points were found, and it is not a substitute for solving
the reliability-constrained problem: SLSQP started from the deterministic
optimum could not move at all (§3.10). §6.4 solves it.

### Discussion

The estimator is first-order and its weakest point is exactly where the
constraint is tightest: $g$ is the most nonlinear constraint and FORM
under-predicts its failure probability, 0.42 against 0.54 sampled. The error is
in the conservative direction, so the conclusion is not fragile even though the
estimator is approximate.

The transferable statement is about the order of the two studies rather than
about this engine. The bounds were fixed first and the tolerance study run on
the result, and by then the specification contained one requirement finer than
the model's own resolution and one that governed the reliability of everything
else — neither visible in any nominal quantity. Running the tolerance study
against the *specification*, before any design exists, costs one Jacobian and
answers a question the optimizer never asks: which of these numbers the
specification is entitled to contain.

## 6.3 Against a conventional engine

### Result

Both engines at the same compression ratio, clearance volume and fuel per
cycle, sized by identical code, and both optimised. The baseline has two
freedoms, rod obliquity and speed, and is optimised twice: once under the
requirements a conventional engine actually has, and once under the EX-link's
own rod-angle and side-load caps.

| | slider-crank, its own limits | slider-crank, same specification | EX-link |
|---|---|---|---|
| $r/l$, speed | 0.195, 2151 rpm | 0.0959, 1227 rpm | — , 1000 rpm |
| power strokes per minute | 1076 | 614 | 1000 |
| rod angle / $\gamma$ | 11.2° / 0.039 | 5.50° / 0.0200 | within 10° / 0.02 |
| brake efficiency | 0.359 | 0.347 | 0.407 |
| range | 2888 km/L | 2467 km/L | 3395 km/L |
| EX-link advantage | +15.6 % | +37.6 % | — |

### The advantage survives matching the firing rate

The EX-link completes its four strokes in one crankshaft revolution where a
four-stroke needs two, so at equal shaft speed it fires twice as often and the
comparison would be measuring the cycle rate. It is not: the own-limits
baseline's optimum is 2151 rpm, which is **1076 power strokes per minute against
the EX-link's 1000**. The 15.6 % above is therefore already an essentially
rate-matched comparison, and it is the baseline that fires more often.

Matching the rate exactly — the baseline pinned to twice the EX-link's speed,
its obliquity re-optimised there — holds the result across the rates the
EX-link runs at:

| power strokes / min | EX-link, rpm | km/L | slider-crank, rpm | $r/l$ | km/L | advantage |
|---|---|---|---|---|---|---|
| 800 | 800 | 3350 | 1600 | 0.178 | 2822 | **+18.7 %** |
| 1000 | 1000 | 3395 | 2000 | 0.191 | 2884 | **+17.7 %** |
| 1200 | 1200 | 3363 | 2400 | 0.200 | 2878 | **+16.9 %** |
| 1400 | 1400 | 3181 | 2800 | 0.214 | 2823 | **+12.7 %** |

Between 800 and 1400 fires per minute the advantage is +13 % to +19 % and never
changes sign. The EX-link design is the one optimised at 1000 rpm and is merely
re-scored at the other three rates, so those rows understate it; the baseline is
re-optimised at every row and does not.

This is the comparison to quote when the question is whether the mechanism is
worth building, because it is the one that holds the useful output constant.

### Where the advantage comes from is a different question

Firing rate and extended expansion arrive together in this topology, and
{py:func}`~exlink.slidercrank.firing_frequency_sensitivity` separates them by a
counterfactual on the EX-link itself: keep the linkage and its motion, but give
it a conventional two-revolution gas exchange — twice the friction per cycle and
half the power at the same speed — and re-score through the vehicle. The
advantage over the own-limits baseline goes from +15.6 % to **−4.3 %**.

That answers "which feature carries the advantage", and the answer is the firing
rate rather than the extended expansion: the longer expansion alone does not pay
for four extra journals and a gear train. It is not a comparison between two
engines, because the engine it scores does not exist — it is the EX-link with
the property that distinguishes it removed. As an attribution it is
informative; as a verdict it would be wrong, and the matched-rate table above is
the verdict.

### Held to the EX-link's own limits the gap widens

`evaluate_slidercrank` tests convergence, net work and the speed rule, never the
10° rod angle or the 0.02 side-load ratio, and the own-limits optimum sits at
11.2° and 0.039 — roughly twice each cap. Imposing them,
{py:func}`~exlink.slidercrank.optimise_slidercrank_to_specification` reaches
2467 km/L: the side-load cap binds, needs a connecting rod ten times the crank
throw, and costs the baseline 15 % of its range. The advantage becomes
**+37.6 %**.

Two details decide that number and are worth stating, because getting either
wrong moves it by ten points.

*Which $\gamma$.* The EX-link's cap is applied to a quasi-static quantity —
{mod}`exlink.loads` propagates the gas load alone. At the operating speed the
reciprocating inertia opposes the gas side load over part of the cycle, so the
same mechanism reads 0.0248 quasi-statically at $r/l = 0.12$ and 0.0181 at
1500 rpm. Capping the *at-speed* ratio would admit proportions the EX-link is
not allowed and hand the baseline about a tenth of its range.

*Which search.* Range rises with obliquity and the cap binds from above, so the
optimum lies against a boundary in one variable and at an interior maximum in
the other, and the function grids the box before refining. Stopping a fraction
inside the cap, at $r/l = 0.095$ and 1500 rpm, returns 2372 km/L — 4 % low, and
low in the direction that flatters the linkage.

This baseline runs at 614 power strokes per minute against the EX-link's 1000,
so unlike the matched-rate table it is *not* rate-matched. Matching it would
widen the gap further, not narrow it: at 2000 rpm the same caps admit nothing
better than 988 km/L, because the cap is quasi-static while the friction and
inertia that punish the short rod at speed are not.

### The constrained baseline is dominated too

The baseline held to the same specification converges onto its side-load cap, at
$\gamma = 0.02000$ against a bound of 0.02, which is precisely where §6.2 shows
reliability is lost:

| $r/l$ | quasi-static $\gamma$ | range | $P_f$ | $\beta$ |
|---|---|---|---|---|
| **0.09593** (the optimum) | 0.02000 | **2467.5 km/L** | 0.595 | $-0.24$ |
| 0.09590 | 0.01999 | 2467.2 km/L | 0.106 | 1.25 |
| 0.09580 | 0.01997 | 2466.2 km/L | $9\times10^{-9}$ | 5.63 |
| 0.09550 | 0.01991 | 2463.4 km/L | $< 10^{-16}$ | $\ge 8.2$ |

The optimum sits on the constraint, so about three builds in five miss it.
Giving up 0.4 % of the obliquity costs **0.17 % of range** and takes the failure
probability below the estimator's floor. This is §6.2's dominated optimum again,
on a mechanism with two toleranced dimensions instead of eleven and a search
over two variables instead of ten — which is the strongest evidence in the study
that the effect belongs to *deterministic optimization under tolerance* rather
than to this linkage.

Once both are backed off their boundaries the baseline is the more reliable
mechanism, at $\beta \ge 8.2$ against 3.0. Not because it is better designed: it
has two toleranced lengths against eleven, so it has fewer ways to be wrong. Its
binding constraint is a function of $r/l$, and machining error on a 28 mm crank
and a 292 mm rod moves that ratio by about a tenth of a percent.

### Discussion

The measurements point different ways on purpose, and the useful statement is
their conjunction: **the linkage buys range by adding degrees of freedom, and
pays for them in the probability that all of them land in tolerance at once.**
That trade is invisible to any comparison scoring only the nominal design, and
it is the counterweight to §6.1's argument that more freedom buys a better
optimum.

| question | comparison | figure |
|---|---|---|
| is the mechanism worth building? | matched power strokes per minute | **+13 to +19 %** |
| optimum against optimum, each at its own best speed | own-limits baseline | **+15.6 %** |
| which feature carries the advantage? | EX-link without its firing rate | **−4.3 %** |
| what if the baseline had to meet the same brief? | same rod-angle and side-load caps | **+37.6 %** |

Which of the two range comparisons — against the baseline under its own limits,
or against the baseline under the linkage's — is the right question is a
judgement this study does not make. The 10° and 0.02 limits come from the EX-link's brief, and practical
slider-cranks run 14–19° routinely, so holding one to 5.5° may be imposing an
alien specification. What is not defensible is holding one engine to limits the
other is silently exempt from.

Two conservatisms run against the EX-link and are not quantified here: no gas
exchange is modelled, and the loss that omits is about 2.5× larger for the
conventional engine (§7.2); and the reliability figures compare mechanisms of
different dimensionality, which is a real difference rather than an artefact but
is not like-for-like the way the range figures are.

## 6.4 The announced problem, solved

### Result

§3.10 states a problem: maximise range, hold every constraint, constrain a
system probability of failure. Solving *that* needs the bounds §6.2 identifies,
because at the bounds as written no design reaches the target. The run used the
gap at 0.054 mm and both bands at $\pm 0.15$:

| | start (`COUPLED_DESIGN`) | result |
|---|---|---|
| range | 3338 km/L | **3395 km/L** |
| system $P_f$ | $1.0\times10^{-3}$ | $1.3\times10^{-3}$ |
| system $\beta$ | 3.08 | **3.00**, on its target |
| worst constraint | — | $-2.2\times10^{-7}$ |

1352 evaluations, 61 minutes, the iteration cap reached rather than a
convergence test — so this is a lower bound.

### What the reliability requirement costs

Imposing the constraints *without* it reaches 3501 km/L, 3 % more. It gets
there by converging onto its active constraints, which §6.2 shows is what
destroys reliability. The 3 % is the price of standing off the boundary, and
the two figures answer different questions: 3501 is the best nominal design,
3395 the best that also survives its own manufacturing scatter.

### Which relaxation the result depends on

Only the band. Re-scoring the same design against the 0.1 mm gap accepted in
§6.2 gives $P_f = 1.339\times10^{-3}$ against $1.344\times10^{-3}$ at 0.054 mm —
a difference of four parts in a thousand, because the gap sits at $\beta = 8.1$
either way and contributes nothing. The band is a different matter: at
$\pm 0.05$ this design is not admissible at all, and at $\pm 0.12$ the target
would have to be $\beta = 3.09$ rather than 3.00.

| | value at the solution | $\sigma$ | $\beta$ |
|---|---|---|---|
| `stroke_upper` | $-0.054$ mm | 0.018 | **3.00** |
| `ratio_upper` | $-0.032$ | 0.005 | 6.94 |
| `tdc_gap` (at 0.1 mm) | $-0.100$ mm | 0.012 | 8.12 |
| every other constraint | — | — | $> 13$ |

One constraint is active in the reliability sense and the rest are spectators,
which is the shape a reliability-constrained optimum should have and a useful
check that the index is being steered rather than merely reported.

### Three defects stood in the way, and one generalises

Each was invisible in the aggregates the runs reported:

| symptom | cause |
|---|---|
| $I = 85.1$ against the 57.6 its gear pair realises | §3.7 makes $I$ an output of the catalogue choice; the search treated it as a variable |
| $\beta$ pinned at $-8.2095$ to fifteen digits in two runs | the orthant integrates to exactly 1 outside the band, so the index goes flat |
| ~2× the necessary MDA calls | SLSQP differences objective and constraints over the same stencil in separate passes |

The second is the one to carry away. **A probability makes a poor constraint
wherever it saturates**: outside the band it says nothing about how far
outside, and a difference quotient straddling the band sees an eleven-unit fall
over a $10^{-5}$ step. The search is steered on
$\min_i \beta_i = \min_i(-g_i/\sigma_i)$, smooth through the band, and the
system probability is *reported* at the solution rather than assumed from the
target. That is weaker than constraining the system index — §3.8 says why — and
it is what makes the problem solvable.

### The fallback earns its place only where it is needed

| start | fell back to the target | outcome |
|---|---|---|
| `COUPLED_DESIGN` (runs) | 0 of 836 | pure range maximisation |
| `REFINED_DESIGN` at 1250 rpm (does **not** run) | 26 of 103 | 0 km/L $\to$ 3336 km/L |

From a start that runs, the ladder never fires and costs nothing. From one
where the engine will not run and km/L does not exist, a quarter of the search
is conducted on the target. The motion residual ends at 5.53 mm — far from the
target — because once the range is computable the optimizer abandons the
prescribed motion entirely. That is how a fallback differs from a constraint.

## 6.5 Supporting measurements

The results above rest on properties of the problem and of the implementation
that are asserted where they are used and measured here: how strongly the
disciplines couple, what the analytic derivatives are worth, what the
decomposition costs against enumeration, whether the optima are global, and what
each reference design is. They are collected rather than interleaved because
none of them is a finding about the engine.

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

The 1.5 % gain is modest, and what it takes to bank it is the point of the exercise.

`RANGE_DESIGN` satisfies every inequality, including the gap, at `g = 0.0009 mm`. It misses
the two *relaxed equalities* by 1.5 × 10⁻⁴ mm and 6.1 × 10⁻⁵ — SLSQP stopping within its own
convergence tolerance of the constraint it was handed. For scale, the tolerance study puts
the machining standard deviation of `STE` at 0.020 mm, **130 times larger**; no real part
would tell the two apart.

Projecting it back onto the equality manifold is exact and cheap:
`project_onto_equalities` takes the minimum-norm Newton step from the analytic Jacobians. The
step is a few hundredths of a millimetre, it lands the equalities at 1.0 × 10⁻⁴ and
3.1 × 10⁻⁵, and it moves `g` from 0.0009 to **0.0201 mm**:

| | equality residuals | `g` | worst inequality, gap at 0.01 mm | at 0.1 mm | range |
|---|---|---|---|---|---|
| `RANGE_DESIGN` | $5.0\times10^{-2}$, outside the band | 0.0009 mm | $-0.0006$ | $-0.0006$ | 3388 km/L |
| projected onto the equalities | $1.0\times10^{-4}$, $3.1\times10^{-5}$ | 0.0201 mm | $+0.0101$ | $-0.0075$ | 3388 km/L |

Which of the two bounds is written down decides whether the projected design exists. At
0.01 mm the equality manifold and the region `g ≤ 0.01` intersect in a sliver thinner than the
Newton step that reaches the manifold, so the best strictly feasible design remains
`COUPLED_DESIGN` and the 1.5 % is unreachable. At the 0.1 mm bound of §6.2 the projected design is
strictly feasible and the 1.5 % is simply banked, at a cost of 0.47 % of range in the
relaxation — a net gain of about one per cent.

That trade is not one the optimizer could have found. It came out of the tolerance study, and
it is the argument for running one against the specification before the design rather than
against the design afterwards.

---

Next: [7. Conclusions](conclusions.md)
