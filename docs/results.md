# 6. Results and discussion

Each result is stated, then supported, then discussed. Every number is computed
by the code of §4 and pinned by a test; none is transcribed. §7 collects the
limitations.

![each formulation's final design, turning on a common scale](figures/formulations.gif)

*What each objective converged to, from the same starting point, drawn at one
scale and one crank angle. Left to right: the geometric objective under an
augmented Lagrangian and under SLSQP with exact gradients, then minimum coupled
mass, then the study's result — maximum range under a reliability constraint,
3395 km/L. The two geometric optima are long-limbed and stand their cylinder
high — §6.1 is about why — while the two that can see mass are visibly shorter
and squatter. Regenerate with `exlink animate --formulations`.*

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
| **range + reliability, relaxed bounds** | **3395 km/L** | **§3.10's third form; $P_f = 9.4\times10^{-3}$ (§6.8)** | no — relaxed spec |
| slider-crank, optimised | 2888 km/L | the baseline of §6.3, optimised over its own two freedoms | its own limits |
| slider-crank, capped | 2467 km/L | the same, held to the linkage's limits; §6.2 only | yes, on the cap |

The comparison of §6.3 uses the first of those two and reports one figure,
**+17.6 %** — the study's result against the baseline's. The second slider-crank
row is not a competing comparison: it exists because a baseline forced onto an
active constraint is what §6.2 needs to show that the dominated-optimum effect
is not peculiar to the linkage.

Every range in the table is at one operating point, with solid round members.
Two later sections re-score the first comparison under a changed assumption and
report their own pairs, which are not interchangeable with the figures above:
§6.5 over a schedule of four speeds (3260 against 2733 km/L, **+19.3 %**) and
§6.6 with tubular members on both engines, each re-optimised (3425 against
2929 km/L, **+17.0 %**).

Speeds are quoted at the crankshaft, which turns twice per cycle on both
mechanisms; the ``speed_rpm`` the code takes for the linkage is the half-speed
shaft's, at half the quoted figure (§5.3). Figures are given to four significant
digits once, here, and rounded elsewhere.

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
better efficiency, at 2000 rpm.

### Why

The same proximity that lengthens the lever arm amplifies the accelerations:
joint $A$ sees 75 times the crank pin's. Since $m \sim (Ca)^3$ (§3.2), and every
inertia load scales as $\Omega^2$, structural mass grows as the **sixth power of
speed**:

| speed | moving mass | peak bearing load |
|---|---|---|
| 0 rpm | 0.25 kg | 7.7 kN |
| 2000 rpm | 1.03 kg | 12.5 kN |
| 3000 rpm | 8.43 kg | 245 kN |
| 4000 rpm | *no section is thick enough* | |

### Discussion

This is the clearest result in the study and the least fragile. It rests on the
equilibrium solve, which is verified against virtual work to machine precision
(§4.5), and the mechanism is understood rather than merely observed.

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
from a 0.664 probability of missing a requirement to $1.9\times10^{-5}$. The
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

§3.8 computes a probability rather than a margin, over the constraints whose
uncertainty $\Sigma$ actually carries. Evaluated on `COUPLED_DESIGN`:

| gap bound | band | system $P_f$ | $\beta$ | binding |
|---|---|---|---|---|
| 0.010 mm | $\pm 0.05$ | 0.664 | $-0.42$ | `tdc_gap` |
| 0.054 mm | $\pm 0.05$ | 0.329 | 0.44 | `stroke_lower` |
| 0.100 mm | $\pm 0.05$ | 0.328 | 0.44 | `stroke_lower` |
| 0.100 mm | $\pm 0.12$ | $1.0\times10^{-3}$ | 3.09 | `stroke_lower` |
| 0.100 mm | $\pm 0.15$ | $1.9\times10^{-5}$ | 4.12 | `stroke_lower` |

The second row is the useful one: **once the gap is at 0.054 mm it stops
binding, and no further widening of it changes anything.** The two bottom rows
differ from the third only in the band, and they span four orders of magnitude
of failure probability. What sets the reliability of this mechanism is not the
gap but how nearly the expansion stroke is required to equal 74 mm.

Keeping the correlation is not a formality, though on this design it works in
the reassuring direction: at the first row the system probability is 0.664
against the 0.680 an independence assumption gives, and by the last row the two
agree to three figures. What the correlation costs is one matrix product; what
it buys is that the number is a probability rather than a bound.

### Most of this probability is avoidable, and free

The 0.664 is not the price of the requirements. It is the price of *ignoring
them while optimising*. Sampling 1200 designs about the reference and checking
the best by reliability against the full constraint set (§3.10) gives:

| | worst $\beta_i$ | system $P_f$ | range |
|---|---|---|---|
| the reference design | $+0.213$ | 0.664 | 3338 km/L |
| best sampled, fully feasible | $+0.610$ | **0.393** | 3340 km/L |

The failure probability falls by **41 %**, all twenty-five of the best
candidates are feasible against the full model, and the range does not fall —
the best is 0.05 % higher than the design it replaces. **The deterministic
optimum is dominated on both objectives at once**, and the margin on range is
small enough to say the reliability is had for nothing rather than bought.

That is the standard argument for reliability-based design optimization,
measured on this problem rather than asserted: a deterministic optimizer
converges *onto* its active constraints, because nothing in the formulation
rewards standing off them, and a design sitting exactly on $g = 0$ fails half
the time. Backing off by a few hundredths of a millimetre is nearly free in
range and buys most of the reliability back.

Sampling is how those points were found, and it is not a substitute for solving
the reliability-constrained problem: SLSQP started from the deterministic
optimum could not move at all (§3.10). §6.4 solves it.

### The same happens to the slider-crank, which settles what the effect is

A result measured on one mechanism could be a property of that mechanism. It is
not. {py:func}`~exlink.slidercrank.optimise_slidercrank_to_specification`
optimises the conventional baseline under the linkage's own rod-angle and
side-load caps — a specification it has no reason to meet, imposed here only so
that it has an active constraint to converge onto. It converges onto it, at
$\gamma = 0.02000$ against a bound of 0.02:

| $r/l$ | quasi-static $\gamma$ | range | $P_f$ | $\beta$ |
|---|---|---|---|---|
| **0.09593** (the optimum) | 0.02000 | **2467.5 km/L** | 0.595 | $-0.24$ |
| 0.09590 | 0.01999 | 2467.2 km/L | 0.106 | 1.25 |
| 0.09580 | 0.01997 | 2466.2 km/L | $9\times10^{-9}$ | 5.63 |
| 0.09550 | 0.01991 | 2463.4 km/L | $< 10^{-16}$ | $\ge 8.2$ |

Three builds in five miss the requirement at the optimum; giving up 0.4 % of the
obliquity costs **0.17 % of range** and takes the failure probability below the
estimator's floor. Two toleranced dimensions instead of eleven, two design
variables instead of ten, a different topology and a different optimizer — and
the same shape of result. **The effect belongs to deterministic optimization
under tolerance, not to this linkage.**

That is also why this baseline is not the one §6.3 compares against: 2467 km/L
is what a conventional engine reaches when held to a specification written for
something else, which answers a question about constraints rather than about
mechanisms.

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

Both engines complete their four strokes in 720° of the shaft that power is
taken from, at the same compression ratio, the same clearance volume and the
same fuel per cycle, and are sized by identical structural and tribological
code. Both are
optimised: the baseline over the two freedoms it has, rod obliquity and speed.

| | slider-crank, optimised | EX-link (`RELIABLE_DESIGN`) |
|---|---|---|
| $r/l$ | 0.195 | — |
| crankshaft speed | 2151 rpm | 2000 rpm |
| power strokes / min | 1076 | 1000 |
| indicated efficiency | 0.457 | 0.480 |
| mechanical efficiency | 0.787 | 0.865 |
| brake efficiency | 0.359 | 0.416 |
| engine mass | 16.9 kg | 12.9 kg |
| range | 2888 km/L | **3395 km/L** |

**+17.6 %**, and that is the whole comparison. There is no firing-rate
correction to make, because taking the power off the shaft that turns twice per
cycle leaves both engines with the same relation between speed and cycles; and
there is no speed mismatch to argue about, because re-scoring the baseline at
the linkage's own 2000 rpm gives 2883 km/L rather than 2888, which moves the
figure to +17.8 %.

Against `COUPLED_DESIGN`, the strictly feasible minimum-mass design rather than
the study's result, the same comparison gives +15.6 %. Both are quoted in §7.5;
the 17.6 % is the one that compares each study's best.

### Why

The linkage is ahead on all three terms of the objective at once, which is worth
separating because the three come from different physics.

| | slider-crank | EX-link | why |
|---|---|---|---|
| indicated efficiency | 0.457 | 0.480 | expansion through 20.8 volumes against 16.0 |
| mechanical efficiency | 0.787 | 0.865 | 9.7° of rod angle against 11.2°, and a side-load ratio of 0.018 against 0.039 |
| engine mass | 16.9 kg | 12.9 kg | a flatter torque curve needs less flywheel |

The thermodynamic term is the one the topology exists for and it is the
smallest of the three: five per cent of indicated efficiency. The other two are
consequences of where the optimizer put the linkage rather than of extended
expansion, and §6.1 is the reason it could go there at all — a mechanism with
eleven dimensions can be placed off its singularity, and a slider-crank with
two cannot be placed anywhere its obliquity does not already put it.

### Discussion

Two conservatisms run against the linkage and are not quantified here. No gas
exchange is modelled, and §7.2 measures the loss that omits as about 2.5× larger
for the conventional engine, because an over-expanded charge reaches the exhaust
valve nearer ambient. And the baseline meets neither the 10° rod-angle cap nor
the 0.02 side-load cap the linkage is held to; those come from the linkage's own
brief and practical slider-cranks run 14–19° routinely, so imposing them on it
would be imposing an alien specification. §6.2 reports what happens when they
are imposed anyway — the point there is about reliability, not about the range
comparison.

One conservatism runs the other way, and it is the reason the mass column is not
quite like-for-like: the baseline carries a flywheel sized for a single-cylinder
four-stroke's turning-moment diagram, which is the worst case in this class, and
the linkage's flatter curve is a genuine property of the topology rather than a
modelling artefact — but both flywheels are sized by the same rule at their own
shaft speeds, so the comparison is at least consistent.

## 6.4 The announced problem, solved

### Result

§3.10 states a problem: maximise range, hold every constraint, constrain a
system probability of failure. Solving *that* needs the bounds §6.2 identifies,
because at the bounds as written no design reaches the target. The run used the
gap at 0.054 mm and both bands at $\pm 0.15$:

| | start (`COUPLED_DESIGN`) | result (`RELIABLE_DESIGN`) |
|---|---|---|
| range | 3338 km/L | **3395 km/L** |
| worst constraint | — | $-2.2\times10^{-7}$ |
| system $\beta$ **as the solve measured it** | 3.08 | 3.00, on its target |
| system $\beta$ **as §6.8 corrects it** | 4.12 | **2.35** |
| system $P_f$, corrected | $1.9\times10^{-5}$ | $9.4\times10^{-3}$ |

1352 evaluations, 61 minutes, the iteration cap reached rather than a
convergence test.

The two $\beta$ rows are the important part of this table, and the second is
the one to believe: it is confirmed against 150 000 exact builds at
$9.28\times10^{-3}$ (§6.8). The solve did hold its target — on the constraint
set it was given, which linearised the expansion stroke on the branch that
attains the maximum. It is the *other* branch the parts breach, and on the
corrected set the result stands at $\beta = 2.35$ and **misses the
$10^{-3}$ target by an order of magnitude**.

Two things follow, and they pull apart. The range is unaffected: 3395 km/L is
a measurement of the design, not of the estimator, and every figure in §6.3,
§6.5 and §6.6 stands. What does not stand is the claim that this design meets a
reliability requirement. The re-solve under the corrected constraint is below.

### What the reliability requirement costs

Imposing the constraints *without* it reaches 3501 km/L, 3 % more. It gets
there by converging onto its active constraints, which §6.2 shows is what
destroys reliability. The 3 % is the price of standing off the boundary, and
the two figures answer different questions: 3501 is the best nominal design,
3395 the best that also survives its own manufacturing scatter.

### Which relaxation the result depends on

Only the band. Re-scoring the same design against the 0.1 mm gap accepted in
§6.2 gives $P_f = 9.382\times10^{-3}$ against $9.393\times10^{-3}$ at 0.054 mm —
a difference of one part in a thousand, because the gap sits at $\beta = 8.1$
either way and contributes nothing. The band is a different matter: at
$\pm 0.05$ this design is not admissible at all, at $\pm 0.12$ it reaches only
$\beta = 0.19$, and at $\pm 0.20$ it would reach 4.53.

| | value at the solution | $\sigma$ | $\beta$ |
|---|---|---|---|
| `stroke_upper_2` | $-0.054$ mm | 0.023 | **2.36** |
| `stroke_upper_1` | $-0.054$ mm | 0.018 | 3.00 |
| `ratio_upper_2` | $-0.032$ | 0.006 | 5.14 |
| `tdc_gap` (at 0.1 mm) | $-0.100$ mm | 0.012 | 8.12 |
| every other constraint | — | — | $> 10$ |

Two constraints are active in the reliability sense and they are the same
requirement measured from the two top dead centres — which is the whole of
§6.8. The solve steered on the second row and the first is what governs. Everything
below $\beta = 5$ is a spectator, which is otherwise the shape a
reliability-constrained optimum should have.

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
| `REFINED_DESIGN` at 2500 rpm (does **not** run) | 26 of 103 | 0 km/L $\to$ 3336 km/L |

From a start that runs, the ladder never fires and costs nothing. From one
where the engine will not run and km/L does not exist, a quarter of the search
is conducted on the target. The motion residual ends at 5.53 mm — far from the
target — because once the range is computable the optimizer abandons the
prescribed motion entirely. That is how a fallback differs from a constraint.

## 6.5 The result over a schedule, not a point

### Result

Everything above is reported at one crankshaft speed. §7.2 lists that as a
limitation, and the objection is specific: an optimum found at one point may be
an artefact of the point. The test is to score the same design over a schedule
of points, and to score it as one engine — which is the part that is easy to get
wrong.

A four-point schedule spanning 0.8 to 1.4 times the design speed, weighted
towards it, gives:

| crankshaft rpm | share | EX-link km/L | slider-crank km/L |
|---|---|---|---|
| 1600 | 0.20 | 3347 | 2811 |
| 2000 | 0.40 | 3315 | 2772 |
| 2400 | 0.25 | 3243 | 2705 |
| 2800 | 0.15 | 3046 | 2589 |
| **cycle** | | **3260** | **2733** |
| worst / best | | 0.910 | 0.921 |
| one engine | | 18.85 kg | 29.53 kg |
| of which flywheel | | 15.84 kg | 28.29 kg |

Neither design was tuned for the schedule; both are the point optima of §6.3
scored over it. The EX-link gives up 4.0 % of its single-point range and the
slider-crank 5.4 % of its own, so the advantage does not close but **widens,
from +17.6 % to +19.3 %**.

### Why the aggregation is a harmonic mean

The figure of merit is distance per unit fuel, so it is *fuel per unit
distance* that adds. Over points carrying distance shares $w_i$,

$$\frac{1}{R} = \sum_i \frac{w_i}{R_i}, \qquad \sum_i w_i = 1$$

which is the distance-weighted harmonic mean, not the arithmetic mean of the
$R_i$. The distinction is not cosmetic: the arithmetic mean flatters a design
whose range collapses at one point, and catching exactly that is what a
schedule is for. Here it costs 1.5 km/L on the EX-link — small, because no
point of this schedule collapses, which is itself the finding.

### Why one engine, and which point sizes what

Scoring each point with `evaluate` and averaging is wrong, and wrong in the
optimistic direction, because it lets the engine change between points: the
sections shrink where the inertia is small, the flywheel shrinks where the
speed is high. The sweep in §6.1 reads that way — 12.9 kg at 2000 rpm, 8.1 kg
at 2800 — and it is a sweep of *different engines*, which is legitimate for
locating an optimum and illegitimate for scoring a schedule.

One engine has to satisfy every point, and the two ends of the schedule size
different parts of it, in opposite directions:

- **the fastest point sizes the structure**, because every inertia load grows
  as $\omega^2$ and the sections follow it;
- **the slowest point sizes the flywheel**, because the inertia needed to hold
  a given speed fluctuation goes as $\omega^{-2}$.

`score_cycle` therefore builds the structure at 2800 rpm, takes the flywheel
requirement from whichever point demands most — 1600 rpm, for both mechanisms —
and scores every point against that one mass. Composing the two is exact rather
than conservative in this model: the flywheel is concentric with its shaft, so
a larger one adds no inertia force and cannot feed back into the member loads.

The composition is what the schedule costs. The EX-link's structure at 2800 rpm
is 3.0 kg against 2.8 kg at the design point — the sections barely notice — but
its flywheel goes from 10.1 kg to 15.8 kg, and that is the whole of the 4.0 %.
An engine specified over a range of speeds is a flywheel problem.

### Why the gap widens

Not because the EX-link's flywheel advantage grows: the two wheels stand in the
same ratio at every speed, 18.1 kg against 10.1 kg at 2000 rpm and 28.3 against
15.8 over the schedule, because both scale as $\omega^{-2}$ from a torque curve
whose shape does not move. The ratio is fixed at 1.79 and the schedule does not
touch it.

What the schedule changes is the *mix*. The two engines split their mass in
opposite ways:

| | EX-link | slider-crank |
|---|---|---|
| structure, at 2800 rpm | 3.01 kg | 1.25 kg |
| flywheel, over the schedule | 15.84 kg | 28.29 kg |
| total | 18.85 kg | 29.53 kg |

The slider-crank is the *lighter* mechanism — two members and no gear pair
against a trigonal link, a rocker and a pair of gears — and it is the heavier
engine, because the item it loses on is the flywheel and the flywheel is most
of both engines. Taking the schedule down to 1600 rpm makes the flywheel a
larger share of each total: 84 % of the EX-link's mass and 96 % of the
slider-crank's, against 78 % and 92 % at each engine's own single point. So the comparison is
weighted further towards the one item the EX-link wins by a factor of 1.79, and
the total-mass ratio moves from 1.31 at a point to 1.57 over the schedule.

That is a general property, not an accident of this schedule: **the wider the
speed range an engine must cover, the more of its mass is flywheel, and the
more a flat torque curve is worth.** The mechanism whose advantage is
smoothness gains from being asked to work over a range rather than at a point,
and the slider-crank's structural simplicity buys less the wider that range is.

### What the schedule does not settle

The baseline was given its one remaining freedom back: re-optimising $r/l$
against the schedule rather than the point moves it from 0.195 to 0.198 and the
cycle range from 2733.2 to 2733.7 km/L — four parts in ten thousand. Its point
optimum was already its cycle optimum. The EX-link was **not** re-optimised
over the schedule, which would be an eleven-variable solve rather than a scalar
one, so +19.3 % is a lower bound on what a cycle-aware design would reach.

The weights themselves are an assumption, not a measurement. Without a surveyed
track there is no defensible way to derive them, and a different spread would
give a different number. What the schedule establishes is not a more accurate
range but the absence of an artefact: the ratio of worst point to best is 0.91,
the design that wins at 2000 rpm wins at every point of the schedule, and no
part of §6.3's conclusion depended on the speed it was measured at.

## 6.6 What tubes are worth, and why the answer is a speed

### Result

§7.2 named solid round bars the largest single modelling conservatism and put
the cost at "perhaps 30 % of mass". Boring the members — the rods and the
trigonal link, not the crank throws, which are forgings — measures it, and the
guess is wrong in both directions at once.

It is wrong about the mass, because the members are not where the mass is:

| item | mass | share |
|---|---|---|
| flywheel | 10.13 kg | 78.3 % |
| crankcase | 1.01 kg | 7.8 % |
| shafts | 0.62 kg | 4.8 % |
| bearings | 0.46 kg | 3.6 % |
| cylinder head | 0.26 kg | 2.0 % |
| gears | 0.25 kg | 1.9 % |
| **the seven members** | **0.17 kg** | **1.3 %** |
| piston | 0.04 kg | 0.3 % |

Thirty per cent of 1.3 % is four grams of a 12.9 kg engine. If tubes mattered
here it would have to be for some reason other than their own weight.

They do matter, and it is for another reason. Against crankshaft speed, with
`RELIABLE_DESIGN` unchanged and the bore at $k = 0.6$:

| crankshaft rpm | solid km/L | tubular km/L | gain | members solid | members tubular |
|---|---|---|---|---|---|
| 1200 | 3206 | 3205 | −0.0 % | 123 g | 134 g |
| 1600 | 3350 | 3350 | 0.0 % | 130 g | 134 g |
| 2000 | 3395 | 3411 | +0.5 % | 167 g | 134 g |
| 2400 | 3363 | 3423 | +1.8 % | 250 g | 149 g |
| 2800 | 3181 | 3385 | +6.4 % | 405 g | 193 g |
| 3200 | 2482 | 3249 | **+30.9 %** | 702 g | 276 g |

Nothing below the design speed, everything above it.

### Why

The member sections are set by two loads that scale differently. The gas load
does not change with speed and a bore does not reduce it, so at the bottom of
the range the sections are what they were and the bore only forces the wall
floor to grow them slightly. The inertia load is the member's own mass times an
acceleration proportional to $\Omega^2$, and there a bore is compound interest:
lighter members need less section, less section is lighter still, and the
sizing/inertia fixed point of §3.5 converges somewhere else entirely.

Which is to say the tube is not a mass improvement that happens to help at
speed. It is a *conditioning* improvement, and it acts on exactly the mechanism
§6.1 is about.

### It does not repeal §6.1

Sizing the quasi-static optimum `REFINED_DESIGN`, which sits at $W = 0.981$:

| crankshaft rpm | solid, moving mass | peak bearing | tubular, moving mass | peak bearing |
|---|---|---|---|---|
| 0 | 0.16 kg | 7 726 N | 0.16 kg | 7 726 N |
| 2000 | 0.96 kg | 12 629 N | 0.37 kg | 6 356 N |
| 3000 | 8.34 kg | 244 805 N | 2.66 kg | 74 225 N |
| 4000 | *no section is thick enough* | | *no section is thick enough* | |

At 3000 rpm the tube cuts the mass by a factor of three and the bearing load by
3.3, and at 4000 rpm it changes nothing at all: the design is still
unbuildable. That is the correct shape of the result. The sixth-power
divergence of §6.1 is a property of the *conditioning* — the acceleration ratio
of 75 between joint $A$ and the crank pin — and a constant factor on the member
mass moves where the divergence bites without removing it. **The clearest
finding in the study survives its largest modelling conservatism being
removed**, which is more than could be said for it before this was measured.

### Why the bore ratio has an interior optimum

A bore ratio is not free to grow, because a bored member has a wall, and a wall
has a minimum. At $k$ and a floor $t$, no member may be drawn below
$2t/(1-k)$; at $t = 1.5$ mm this is 6 mm at $k = 0.5$ and 15 mm at $k = 0.8$,
and most of this linkage's members are smaller than that when solid. So past
some ratio the floor rather than the load sizes the light members and the
section grows back:

| $k$ | wall floor | km/L at 2400 rpm | members |
|---|---|---|---|
| 0.0 | — | 3363 | 250 g |
| 0.3 | 4.3 mm | 3387 | 210 g |
| 0.5 | 6.0 mm | 3415 | 161 g |
| **0.6** | 7.5 mm | **3423** | 149 g |
| 0.7 | 10.0 mm | 3406 | 174 g |
| 0.8 | 15.0 mm | 3322 | 256 g |

Optimising range over the bore ratio and the speed together lands at
$k = 0.579$ at 2276 rpm, giving 3425 km/L. Against the solid optimum of
3395 km/L at 2000 rpm that is **+0.9 %** — and the engine is 10.4 kg rather
than 12.9 kg, which is where the improvement mostly goes: not into range, into
being able to run 300 rpm faster for the same range.

### The comparison, with both engines tubed

A modelling change that flatters one mechanism and not the other would move
§6.3 without being about the engines at all, so the baseline gets the same
bore, under the same rule — its rod may be tube, its crank throw may not — and
the same freedom to re-optimise:

| | slider-crank | EX-link |
|---|---|---|
| solid | 2888 km/L at $r/l$ = 0.195, 2151 rpm | 3395 km/L at 2000 rpm |
| tubular, re-optimised | 2929 km/L at $r/l$ = 0.176, 2232 rpm | 3425 km/L at 2276 rpm |
| the bore is worth | +1.4 % | +0.9 % |

so the advantage narrows slightly, **+17.6 % to +17.0 %**. The bore helps the
conventional engine marginally more, and for the reason §6.3 already gives: its
two members are a larger share of a smaller mechanism, so removing metal from
them is a larger relative change. It is a sixth of a percentage point of the
17, which is the right size for a conservatism that acts on 1.3 % of one engine
and 0.4 % of the other.

### Discussion

Three things are worth separating.

The first is that a modelling conservatism was priced rather than argued about,
and priced in the currency of the objective. "Members 30 % lighter" was true and
useless; "0.9 % of range, and only above the design speed" is the statement a
design decision can be made against.

The second is that the *shape* of the answer was not guessable from the
conservatism itself. A bore acts on 1.3 % of the engine mass and moves the
optimum operating speed by 300 rpm, because in a coupled problem the size of an
effect is not the size of the thing it acts on. Anything that reads the mass
budget alone would have ruled tubes out.

The third is a limit. The model prices a tube through mass, stiffness and a wall
floor, and not through anything else: it has no local buckling of the wall, no
joint or end-fitting mass where a tube must become solid to take a pin, and no
account of what a bored trigonal link costs to make. All three of those work
against the tube, so +0.9 % is an upper bound on a real one.

## 6.7 Reliability against more than the dimensions

### Result

§6.2 and §6.4 price manufacturing scatter on the eleven dimensions, and §3.10
records the limitation that follows: five of the thirteen constraints get no
probability at all, because they are not functions of those eleven. Widening
the uncertain vector to seventeen — adding the material's yield strength,
ultimate strength, stiffness and density, the friction coefficient, and the
explosion ratio that sets the gas load — answers two questions, and the answers
point in opposite directions.

**Nothing that was already priced moves.** Every one of the eight geometric
constraints turns out to be 100 % dimensional to the last figure reported:

| constraint | $\sigma$, dimensions only | $\sigma$, widened |
|---|---|---|
| rod_angle | 0.00637 | 0.00637 |
| compatibility | 0.00005 | 0.00005 |
| tdc_gap | 0.01231 | 0.01231 |
| side_load | 0.00012 | 0.00012 |
| stroke_upper | 0.01814 | 0.01814 |
| ratio_upper | 0.00463 | 0.00463 |

That is exact and not approximate, and deliberately so: the widened model takes
these rows from the same analytic Jacobian rather than differencing them, so a
disagreement here would have been physics rather than a step size. Widening the
uncertain vector therefore changes nothing about the geometric constraints —
which is not to say those constraints were right, and §6.8 shows one of them
was not. The two findings are independent: this section says the *inputs* were
complete enough, and §6.8 says one of the *constraints* was linearised at a
kink.

**What was not priced decides the design.** The system probability goes from
$9.4\times10^{-3}$ — the geometric figure, §6.8's corrected one — to
$5.0\times10^{-1}$, and all of it is one constraint:

| constraint | $\beta$ | dominant source of its variance |
|---|---|---|
| saturation | 460 | explosion ratio, 87 % |
| slenderness | 63 | ultimate strength, 94 % |
| bearing | 36 | explosion ratio, 100 % |
| runs | 26 | friction 51 %, explosion ratio 49 % |
| **gear** | **0.00** | explosion ratio 71 %, yield strength 29 % |

The gear pair the mixed-integer master picked — $m = 0.8$, $z = 48$ — sits
*exactly* on its face-width limit, at 11.99997 against a bound of 12. Whether
it fits is a coin flip.

### Why the narrow model could not have found this

Not because it was careless about the gears, but because the face width is not
a function of any of the eleven dimensions. It is set by the torque the pair
transmits and the strength of what it is cut from, and the pair itself is
chosen by the master problem from a catalogue. A tolerance study on the
linkage's lengths and angles has nothing to say about it, and would not have
said anything however finely those lengths were held.

The general shape is worth naming, because it is not specific to gears. **A
reliability model prices what its uncertain vector contains, and reports
silence as safety on everything else.** §3.10's honest statement that only
seven constraints could carry a probability reads, in hindsight, as a warning
that was not loud enough: the constraint that governs was in the other group.

### What it costs to fix: nothing

The mixed-integer master of §6.8 chose $m = 0.8$, $z = 48$; the exhaustive
search over the same four candidates chose $m = 1.0$, $z = 39$, and was already
known to be 0.6 % better at its own optimum. Re-scoring `RELIABLE_DESIGN`
itself on the two:

| gear pair | range | engine mass | $P_f$ | binding constraint |
|---|---|---|---|---|
| $m = 0.8$, $z = 48$ (pinned, the master's) | 3394.9 km/L | 12.94 kg | $5.0\times10^{-1}$ | gear, $\beta = 0.00$ |
| $m = 1.0$, $z = 39$ (the exhaustive answer) | 3395.3 km/L | 12.91 kg | $9.4\times10^{-3}$ | stroke, $\beta = 2.36$ |
| $m = 2.0$, $z = 19$ (chosen by the sizing rule) | 3396.0 km/L | 12.86 kg | $9.4\times10^{-3}$ | stroke, $\beta = 2.36$ |

The ordering is the same on every column. The pinned pair is the shortest
ranged, the heaviest, and the only one whose gear constraint binds; either of
the other two recovers the geometric probability exactly, because with the gear
constraint at $\beta \ge 4$ the widened model has nothing left to add. That
geometric figure is $9.4\times10^{-3}$ rather than the $1.3\times10^{-3}$ this
section originally reported, for a reason that has nothing to do with widening
the uncertain vector — see §6.8.

There is no trade to negotiate here; the pinned choice was simply worse on all
three counts, and the range gap of 0.03 % is far too small to have revealed it.
Only the widened model could.

The recommendation is one line: **the study's result should not be built with
the pair the master returned.** Its range is 3395 km/L whichever pair is
fitted, to four figures, and its reliability is either $9.4\times10^{-3}$ or a
coin flip depending on which.

### Where each uncertainty actually goes

The variance split is exact rather than attributed — at first order and with
the entries independent, the variance is a plain sum of $(\partial g/\partial
u_j \, \sigma_j)^2$, with no interaction term to allocate:

| constraint | dimensions | yield | ultimate | friction | explosion ratio |
|---|---|---|---|---|---|
| the eight geometric | 100 % | — | — | — | — |
| saturation | 0 % | — | 13 % | — | 87 % |
| slenderness | 0 % | — | 94 % | — | — |
| bearing | 0 % | — | — | — | 100 % |
| runs | 0 % | — | — | 51 % | 49 % |
| gear | 0 % | 29 % | — | — | 71 % |

Two entries are worth reading. `slenderness` rides almost entirely on the
*ultimate* strength, not the yield, which says the members are sized by fatigue
rather than by static yield — the endurance limit is $0.5 S_u$ before Marin
corrections, so the scatter that matters is the one in $S_u$. And `runs`
splits evenly between friction and the gas load, which is the mechanical
efficiency of §6.3 restated as a variance: the margin between indicated and
brake work is a difference of two uncertain quantities of comparable size.

Stiffness and density earn their place by not mattering: 0 % and 6 % of one
constraint between them, which is what a CoV of 0.03 and 0.01 buys. Carrying
them costs two of the eighteen analyses and settles the question.

### Discussion

The finding here is not really about gears. It is that **a reliability model
prices what its uncertain vector contains and reports silence as safety on
everything outside it**, and that the constraint which governs a design has no
obligation to be inside. §3.10 said plainly that only some constraints could
carry a probability; that statement was true, correctly worded, and read as a
limitation rather than as the warning it was.

Two limits of what is done here. The uncertain parameters are taken
independent, which is conservative for the two strengths — a low-yield heat
usually has a low ultimate — and neutral for the rest. And the model remains
first order: the widened constraints include ``saturation``, which is a
threshold on a fixed-point iteration and is about as far from linear as a
constraint gets, though at $\beta = 460$ its linearisation error is not what
decides anything.


## 6.8 What sampling found in the reliability estimate

### Result

Every probability in §6.2, §6.4 and §6.7 is first order: the constraints are
linearised at the nominal design and the failure probability read off a
multivariate-normal orthant. §7.2 lists that as a limitation and §7.3 asks for
a sampling check. The check was run — 150 000 exact builds drawn from the same
covariance, each analysed in full — and it does not agree:

| design | specification | FORM | sampled, 150 000 builds | ratio |
|---|---|---|---|---|
| `COUPLED_DESIGN` | as written | 0.6454 | $0.6650 \pm 0.0012$ | 0.97 |
| `RELIABLE_DESIGN` | §6.2's bounds | $1.35\times10^{-3}$ | $(9.28 \pm 0.25)\times10^{-3}$ | **0.145** |

(Standard errors, so the second row's 95 % interval is
$[8.79, 9.77]\times10^{-3}$ — nowhere near FORM's figure.)

FORM is accurate on one design and optimistic by a factor of seven on the
other, and the difference between them turns out to be the whole finding.

### Why

The expansion stroke is
$STE = \max(\lambda_{tdc,1}, \lambda_{tdc,2}) - \min \lambda$: the piston has
two top dead centres per cycle and the stroke is measured from the higher,
because that is the one that sets the clearance volume. **A maximum of two
smooth functions is not smooth at the tie**, and linearising it uses whichever
branch attains it at the nominal design.

`COUPLED_DESIGN` has its two top dead centres 6.7 μm apart, which is comfortably
more than the scatter can bridge, so one branch attains the maximum in
essentially every build and linearising it is right. `RELIABLE_DESIGN` has them
**0.107 μm** apart — and it has them there because driving the top-dead-centre
gap to zero is exactly what its own $g$ constraint rewards. The parts, whose
dimensions scatter by some 8 μm, straddle that tie in every single build.

Measuring the two branches separately, analytically and by sampling, settles it:

| | $\sigma$, analytic | $\sigma$, sampled | skew | $\beta$ | $P$ |
|---|---|---|---|---|---|
| stroke from TDC 1 | 0.01814 | 0.01817 | $-0.01$ | 3.00 | $1.4\times10^{-3}$ |
| stroke from TDC 2 | 0.02304 | 0.02306 | $-0.02$ | 2.36 | $9.1\times10^{-3}$ |
| the maximum of the two | — | 0.02016 | $+0.11$ | — | $9.3\times10^{-3}$ |

Each *branch* is Gaussian to two decimal places in its skew, which is to say
FORM is exactly right on either one. What FORM did was linearise **branch 1**,
which attains the maximum by a tenth of a micron and is 27 % less sensitive than
branch 2. The parts breach branch 2.

So the estimator was never wrong. The *formulation* was: it presented a maximum
as if it were a single differentiable function, and handed the linearisation the
branch that happened to be on top.

### The fix costs one Jacobian row

Carry both branches as separate constraints instead of the one that attains the
maximum. For the upper bounds this is not merely safe but exactly right: the
realised stroke exceeds its bound when *either* branch does, which is a union
over correlated normals, and a union over correlated normals is precisely what
the orthant integral already computes. (For the two lower bounds the honest
statement is an intersection, so carrying both is conservative there; both sit
above $\beta = 10$ at every design in this study, so it costs nothing.)

With that one change:

| design | FORM, one branch | FORM, both branches | sampled | error |
|---|---|---|---|---|
| `COUPLED_DESIGN` | 0.6454 | 0.6636 | 0.6650 | **0.2 %** |
| `RELIABLE_DESIGN` | $1.35\times10^{-3}$ | $9.38\times10^{-3}$ | $9.28\times10^{-3}$ | **1.1 %** |

No sampling, no second derivatives, no extra analyses — one more row of a
Jacobian that was already being computed.

### Second derivatives are the wrong instrument, and the wrongness is measurable

The natural first guess is curvature, and §7.3 asks for $\nabla^2 g$ to correct
for it. That guess is wrong here and it can be shown to be wrong.

Differencing the *analytic* gradient — which is the well-conditioned way to a
Hessian, since the first derivatives are exact — separates the maximum from its
branches in one table. Norms are in the standard normal space, so they are
comparable across rows:

| step | $\lVert\nabla^2 g\rVert$ | asym. | norm × step | | $\lVert\nabla^2 g\rVert$ | asym. | norm × step |
|---|---|---|---|---|---|---|---|
| | **the maximum** | | | | **branch 2 alone** | | |
| 0.50 | 0.006568 | 1.00 | 0.0032840 | | 0.000010 | 0.00 | 0.0000050 |
| 0.20 | 0.016416 | 1.00 | 0.0032832 | | 0.000010 | 0.00 | 0.0000020 |
| 0.10 | 0.032830 | 1.00 | 0.0032830 | | 0.000010 | 0.00 | 0.0000010 |
| 0.05 | 0.065657 | 1.00 | 0.0032828 | | 0.000010 | 0.00 | 0.0000005 |
| 0.02 | 0.164135 | 1.00 | 0.0032827 | | 0.000010 | 0.00 | 0.0000002 |
| 0.01 | 0.328270 | 1.00 | 0.0032827 | | 0.000010 | 0.00 | 0.0000001 |

On the maximum, *norm × step* is constant to five significant figures across a
fiftyfold range of steps and the matrix is 100 % asymmetric at every one of
them. That is a fixed jump in the gradient divided by a shrinking step — there
is no second derivative there to find. A finer crank-angle grid does not change
it, which rules out discretisation.

On the branch, the norm is constant, the matrix is exactly symmetric, and the
value is five orders of magnitude smaller: a genuine second derivative, and a
negligible one. The branches are as close to linear over the tolerance as
anything in this model.

The jump is the branch switch: perturbing $a$ by 0.3 μm moves
$\partial STE/\partial a$ from 0.0093 to 0.52, a factor of 56. **The constraint
surface is not curved, it is kinked**, and no order of Taylor expansion repairs
a kink. What was needed was not a second derivative but a first one, of the
other branch.

### Discussion

Three things generalise past this problem.

**A sampling check is not a refinement, it is a control.** It was run to put an
error bar on a first-order estimate and it found a factor of seven. Nothing in
the FORM output — not the index, not the correlation, not the per-constraint
breakdown — carried any sign that the estimate was wrong, because from inside
the linearisation it was self-consistent.

**Optimisation drives designs onto the non-smooth features of their own
constraints.** The tie was not bad luck. §6.4's objective rewards a small
top-dead-centre gap, so the optimizer closed it to a tenth of a micron, which is
the same as saying it put the design exactly on the ridge where the stroke
constraint stops being differentiable. Any deterministic optimum should be
suspected of sitting on whatever non-smoothness its formulation contains, and
that is the first place a reliability estimate will fail.

**Diagnosis beats correction.** The correction the limitation list called for —
second derivatives — would have cost eleven extra Jacobian evaluations per
design point and produced noise. Finding out *why* the estimate was wrong cost
one sampling run and produced a fix that is exact, free, and applies at every
design rather than only at the one checked.


## 6.9 A third point on the dimensionality axis

### Result

§7.1 explains the linkage's 17.6 % as a consequence of "having eleven
dimensions to place rather than two", and §7.3 asks for a third topology to
turn that contrast into a trend. The explanation is testable without inventing
one: give the conventional engine a third freedom and see whether the range
moves in proportion.

The wrist-pin offset is the cheapest freedom a real engine has. Optimised over
its own variables, at the same compression ratio and through the same code:

| | freedoms | $r/l$ | speed | $d/r$ | mech. eff. | mass | range |
|---|---|---|---|---|---|---|---|
| slider-crank, centred | 2 | 0.1954 | 2151 rpm | — | 0.7867 | 16.93 kg | 2887.7 km/L |
| slider-crank, offset | 3 | 0.1943 | 2157 rpm | $+0.156$ | 0.7921 | 16.84 kg | **2902.5 km/L** |
| EX-link | 11 | — | 2000 rpm | — | 0.865 | 12.94 kg | 3394.9 km/L |

The third freedom is worth **+0.51 %**.

### The dimension count is not the mechanism

Take the explanation literally and it predicts badly. One extra freedom buys
0.51 %, so nine of them, extrapolated, would buy about 4.6 % — against the
17.6 % the linkage actually reaches. Dimension-counting under-predicts the
topology by a factor of nearly four, which is enough to say that *counting
dimensions is not the explanation*, only a proxy for it.

What the offset does is instructive precisely because it is so narrow. It moves
one term of the three §6.3 decomposes the advantage into, and moves it a
little:

| | centred → offset | centred → EX-link |
|---|---|---|
| indicated efficiency | 0.457 → 0.457 | 0.457 → 0.480 |
| mechanical efficiency | 0.787 → **0.792** | 0.787 → **0.865** |
| engine mass | 16.93 → 16.84 kg | 16.9 → 12.9 kg |

The offset rearranges the side load, so the mechanical efficiency improves and
the mass follows it slightly. It cannot touch the indicated efficiency, because
the piston still makes two identical strokes per revolution however far the pin
is moved: there is no extended expansion to be had from a slider-crank at any
number of dimensions. And it cannot touch the flywheel, because there is no
half-speed shaft to hang one on.

So the honest statement is narrower than §7.1's and stronger for it. **The
eleven dimensions are what let the linkage exploit its topology, not what give
it the advantage.** A conventional engine handed more freedoms improves along
the one axis its topology leaves open, and stops.

### What this is not

It is a third point on the dimensionality axis, and not a third
extended-expansion mechanism. Two topologies still establish the contrast of
§6.3, and a genuine third — one carrying its own internal 2:1 ratio, which is
what unequal strokes require — would still be needed to make it a trend. What
has been settled here is only that the trend cannot be *assumed* to run with
the number of design variables, which is what the limitation list implied when
it asked for a third topology and what §7.1 implied when it explained the
result by counting them.


## 6.10 Getting the search where it needs to go

Two instruments were built for two obstructions §7.3 records. One removes
its obstruction; the other turns out not to have one to remove, and both
answers are here because a method that was not needed is as much a result as
a method that was.

### Restoration: 0 of 6 becomes 6 of 6

§6.11 reports that manifold-projected restarts reached feasibility in **0 of 6**
attempts on the range problem, which is why §3.9's multistart is inconclusive
there: not that the restarts found worse optima, but that they never found the
feasible set, so there was nothing to compare. §7.3 asks for a restoration
phase. Six starts drawn the same way — 5 % scatter about `COUPLED_DESIGN` —
put through the epigraph problem {py:func}`~exlink.restoration.restore` solves:

| start | worst margin as drawn | after restoration | evaluations |
|---|---|---|---|
| 0 | $-1.54$ | $+7.84\times10^{-3}$ | 302 |
| 1 | $-14.9$ | $+5.81\times10^{-3}$ | 296 |
| 2 | $-6.64$ | $+5.64\times10^{-3}$ | 282 |
| 3 | $-5.01$ | $+5.14\times10^{-3}$ | 285 |
| 4 | $-9.39$ | $+6.46\times10^{-3}$ | 289 |
| 5 | $-5.28$ | $+2.43\times10^{-3}$ | 312 |

**0 of 6 feasible as drawn, 6 of 6 after**, at about 290 evaluations each —
1766 in total, 42 minutes.

### Why it works when the range solve does not

The two solves are asked different questions. A range-maximising SLSQP run
started outside the feasible set has to improve an objective *and* find the
set, and its QP subproblem is built around a point where no step satisfies the
linearised constraints; §3.10 records that it responds by reporting a positive
directional derivative and stopping. The restoration problem has no objective
to trade against — only

$$\max_{X,\,t} \; t \quad \text{s.t.} \quad c_i(X) \ge t$$

— so every constraint pulls the same way and the QP always has somewhere to go.
The starts here begin as far as 14.9 outside, and it makes no difference: the
worst of the six needed 312 evaluations rather than 282.

### What the restored points are, and are not

They are inside, and barely. The margins land between $2.4\times10^{-3}$ and
$7.8\times10^{-3}$, which is what §3.4's thinness means in practice — the
interior of this feasible set is a few thousandths of a millimetre wide, and a
restoration that maximises the smallest margin still stops there. That is
enough for a gradient method to start from, which is the whole purpose, and it
is not enough to call the design robust: §6.2 is about exactly that distinction.

What has **not** been measured here is the second half of the claim. Six
restored starts make the multistart *answerable*; whether the answers agree —
whether the range optimum of §6.4 is global — needs six range solves from these
points, which is six times the cost of §6.4's own run and is not attempted.
The entry in §7.3 is therefore struck for the phase it asked for and the
question behind it stays open.

### Continuation: an instrument for an obstruction that is no longer there

§7.3's other request was a way to *reach* the reliable region, on the evidence
that SLSQP with a reliability target attached returns its starting point
unchanged — not for a demanding target but for a step of 0.17. A continuation
in $\beta$, warm-starting each rung from the last, is the natural instrument:
it replaces one impossible step with several possible ones.

Measured against a single solve at the same iteration budget, from
`RANGE_DESIGN` at the bounds §6.2 settles on:

| | evaluations | steered $\min_i\beta_i$ | range | verdict |
|---|---|---|---|---|
| one solve, 150 iterations | 1435 | **3.000** | 3388.4 km/L | on target |
| five rungs, 30 iterations each | 2060 | 2.651 | 3402.1 km/L | short |

The ladder is beaten on both counts: 44 % more evaluations for the same 150
iterations, and it ends 0.35 short of a target the single solve reaches
exactly. The rungs themselves behave as designed — each moves, and the climb is
monotone in reliability from $\beta = 0.001$ to 2.651 while the range *rises*
from 3388 to 3402 km/L — but thirty iterations is not enough for any rung to
converge, so the budget is spent restarting rather than arriving.

The reason is not that continuation is a bad instrument. It is that **the
obstruction it was built for is absent**. §3.10's non-movement was measured
under the first of that section's three forms, with the coupled and vehicle
constraints bound only at the end; under the second form, with everything
imposed and the branch-aware constraint set of §6.8, the search moves freely
from the same start and reaches its target unaided. Given a search that can
move, subdividing its target only fragments its budget.

That is worth recording precisely because §7.3 asked for it. The entry
diagnosed a search failure and proposed a search fix; the fault was in the
constraint set, and §6.8's one extra Jacobian row did more for reachability
than a homotopy does.


## 6.11 Supporting measurements

The results above rest on properties of the problem and of the implementation
that are asserted where they are used and measured here: how strongly the
disciplines couple, what the analytic derivatives are worth, what the
decomposition costs against enumeration, whether the optima are global, and what
each reference design is. They are collected rather than interleaved because
none of them is a finding about the engine.

### How coupled the problem is

$\rho$ is the Gauss–Seidel contraction factor of §3.5.

| crankshaft rpm | $\rho$ | sweeps | verdict |
|---|---|---|---|
| 0 | 0.0000 | 2 | weak |
| 1000 | 0.1307 | 9 | moderate |
| 2000 | 0.6513 | 28 | strong |
| 3000 | 0.6819 | 42 | strong |

At rest $\rho = 0$ exactly, which is the sharpest available check that the
measure reflects the physics rather than the solver: with no inertia there is no
path from mass to load.

### What the derivatives buy

Minimising total moving mass at 2000 rpm, subject to every constraint and a
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

The 0.6 % is not the whole cost of stopping early. §6.7 shows the pair the
master chose sits exactly on its face-width limit, so the design built with it
has a 50 % chance of not fitting; the pair the exhaustive search chose is
0.4 km/L better *and* has $\beta = 3.97$ on that constraint. The master's
answer was worse on both counts, and the range gap alone did not say so.

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
| `RELIABLE_DESIGN` | 25.14 % | 238 | 128 | 0.9364 | 0.00011 | yes, at the §6.2 bounds |

`COUPLED_DESIGN` is the design to compare against: it gives up five points of
$\eta$ to move off the singularity and gets a lighter, faster, longer-ranged
engine for it. `RELIABLE_DESIGN` is the study's result and the design every
figure here is drawn from; it reaches its top dead centres within
$10^{-4}$ mm of each other, so the bound §6.2 spends its length on is not what
constrains it — the band on the stroke is.

### The range optimization

SLSQP on `neg_range`, gear pair pinned, 2000 rpm, started from the coupled reference:

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
