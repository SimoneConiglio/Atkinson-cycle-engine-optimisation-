# Appendix C. Supporting measurements

The results of §5 rest on properties of the problem and of the implementation
that are asserted where they are used and measured here. None is a finding about
the engine.

## C.1 How strongly the disciplines couple

$\rho$ is the Gauss–Seidel contraction factor of §4.3.

| crankshaft rev/min | $\rho$ | sweeps | verdict |
|---|---|---|---|
| 0 | 0.0000 | 2 | weak |
| 1000 | 0.1307 | 9 | moderate |
| 2000 | 0.6513 | 28 | strong |
| 3000 | 0.6819 | 42 | strong |

At rest $\rho = 0$ exactly, which is the sharpest available check that the
measure reflects the physics rather than the solver: with no inertia there is no
path from mass to load.

## C.2 What the analytic derivatives buy

Minimising total moving mass at 2000 rev/min, subject to every constraint and a
25 % efficiency floor:

| | COBYLA | SLSQP + differences | SLSQP + analytic |
|---|---|---|---|
| result | did not move | did not finish | 1.039 → 0.234 kg |
| cost | 120 evaluations | timed out | 40 evaluations, 148 s |

Analytic Jacobians agree with converged central differences to about $10^{-6}$
relative, and GEMSEO's `check_jacobian` passes on every discipline (§4.9).

## C.3 The mixed-integer decomposition

Four gear candidates, 25 SLSQP iterations per sub-problem:

| | chosen pair | range | sub-solves | seconds |
|---|---|---|---|---|
| outer approximation | $m = 0.8$, $z = 48$ | 3366 km/L | **2** | 575 |
| exhaustive | $m = 1.0$, $z = 39$ | 3385 km/L | 4 | 1056 |

Half the sub-solves, 0.6 % short of the best lattice point. The convexification
options of §2.4 were enabled and measured to change nothing here: the master
terminates after two solves, which is less history than the adaptive correction
needs. The 0.6 % is not the whole cost of stopping early — §5.3 shows the pair
the master chose sits exactly on its face-width limit, so the design built with
it has a 50 % chance of not fitting, while the pair the exhaustive search chose is
0.4 km/L better *and* has $\beta = 3.97$ on that constraint. The master's answer
was worse on both counts and the range gap alone did not say so.

## C.4 Local optima

Manifold-projected restarts (§4.6) show the single-start *efficiency* optimum was
local: 30.91 % becomes 36.99 %. That better point is 443 mm tall against 320 and
sits on the $g$ bound — the single-objective efficiency problem is unbounded in
mechanism size, so a stronger search exploits that harder. On the *range* problem,
which is bounded, 0 of 6 restarts reached feasibility at an affordable budget;
§5.5 removes that obstacle with a restoration phase.

## C.5 Reference designs

| design | $\eta$ | $H$ mm | $B$ mm | $W$ | $g$ mm | feasible |
|---|---|---|---|---|---|---|
| `PUBLISHED_DESIGN` | 35.62 % | 283 | 157 | 0.9892 | 8.5236 | no |
| `REFINED_DESIGN` | 27.80 % | 239 | 152 | 0.9811 | 0.0060 | yes |
| `GRADIENT_DESIGN` | 30.91 % | 320 | 159 | 0.9850 | 0.0095 | yes |
| `COUPLED_DESIGN` | 25.00 % | 198 | 131 | 0.9372 | 0.0070 | yes |
| `RANGE_DESIGN` | 25.46 % | 231 | 131 | 0.9319 | 0.0012 | no |
| `RELIABLE_DESIGN` | 25.14 % | 238 | 128 | 0.9364 | 0.000107 | yes, at §5.2's bounds |

`COUPLED_DESIGN` is the design to compare against: it gives up five points of
$\eta$ to move off the singularity and gets a lighter, faster, longer-ranged
engine for it. `RELIABLE_DESIGN` is the study's result and the design every
figure is drawn from.

## C.6 The range optimization, and why the projection matters

SLSQP on `neg_range`, gear pair pinned, 2000 rev/min, started from the coupled
reference, reaches 3388 km/L at 12.47 kg against the start's 3338 km/L at
12.17 kg. The 1.5 % gain is modest and what it takes to bank it is the point of
the exercise.

`RANGE_DESIGN` satisfies every inequality, including the gap at $g = 0.0009$ mm.
It misses the two *relaxed equalities* by $1.5\times10^{-4}$ mm and
$6.1\times10^{-5}$ — SLSQP stopping within its own convergence tolerance of the
constraint it was handed. For scale, the machining standard deviation of
$\mathrm{STE}$ is 0.020 mm, 130 times larger; no real part would tell the two
apart. Projecting back onto the equality manifold with
`project_onto_equalities` is exact and cheap, a minimum-norm Newton step from the
analytic Jacobians of a few hundredths of a millimetre. It lands the equalities
at $1.0\times10^{-4}$ and $3.1\times10^{-5}$ — and moves $g$ from 0.0009 to
0.0201 mm.

Which of the two gap bounds is written down therefore decides whether the
projected design exists at all. At 0.01 mm the equality manifold and the region
$g \le 0.01$ intersect in a sliver thinner than the Newton step that reaches the
manifold, so the best strictly feasible design remains `COUPLED_DESIGN` and the
1.5 % is unreachable. At the 0.1 mm bound of §5.2 the projected design is
strictly feasible and the 1.5 % is banked, at a cost of 0.47 % of range in the
relaxation — a net gain of about one per cent. That trade is not one the
optimizer could have found; it came out of the tolerance study, and it is the
argument for running one against the specification before the design rather than
against the design afterwards.

## C.7 What the prescribed motion established

The fallback objective of §4.7 was developed as a route in its own right before
being demoted to a rung. Four things it established generalise past this problem.

**A target on the manifold need not be reachable.** Both equalities are
functionals of $\lambda$ alone, so a motion with the right two strokes satisfies
them exactly before any linkage exists — measured at 74.0000000000 mm and
16.0000000000. But being on the manifold is a property of the *target* and being
fittable is a property of the *mechanism*. A two-harmonic target is exact and
unreachable: the closest a seven-bar gets is 1.16 mm RMS, which carries the
fitted design outside both tolerance bands. Seeding the target from a motion a
real design produces recovers three orders of magnitude of residual.

**Fitting reaches the manifold; it does not give multistart.** Uniform sampling
finds a design on the equality manifold in 0 of 12 000 draws. Sampling and *then
fitting* reaches it in 22 of 30, and those fits are feasible against the whole
constraint set — but every fit converges to the same linkage whatever start it is
given, the motion very nearly determining the mechanism, so diversity has to come
from varying the target rather than the start.

**What a solve violates is what the formulation left out.** Four apparent
limitations of the prescribed-motion route — that it needed a feasibility
restoration phase, that it was bounded by reachability, that its designs were
infeasible, that its fits could not supply diverse starts — each turned out to be
a constraint absent from the sub-problem rather than a property of the method.
The clearest case is measurable: against a target the mechanism cannot reach, the
unconstrained fit leaves the geometric set by ten to fifteen units, while the same
fit with the inequalities imposed holds every one of them at its boundary.

**A flat region defeats a gradient method, three times over.** The reliability
constraint stalls where $P_f$ saturates; the range objective has nothing to
descend where $R$ does not exist; and the system index is pinned at $-8.2095$ for
every violating design. Each repair has the same shape — supply something that
still varies where the quantity of interest does not: the motion residual as the
objective's middle rung, and $\min_i \beta_i$ in place of the system index.

**The decomposition this suggests.** $\lambda(\theta)$ is the *only* quantity the
linkage sends downstream, so it is the coupling variable, and IDF on it would put
a master choosing $\lambda^\star$ against a sub-problem fitting the linkage to it.
§4.3 rejected IDF because the coupling has 45 367 components — but that count is
*pointwise*, and $\lambda$ is smooth:

| design | harmonics for RMS $< 0.1$ mm | $< 0.01$ mm |
|---|---|---|
| coupled (minimum mass) | 10 | 14 |
| gradient (geometric, SLSQP) | 15 | 23 |

Fourteen to twenty-three coefficients reproduce the motion to 0.01 mm RMS,
tighter than the 0.020 mm machining scatter. A functional IDF would carry of order
30 to 50 consistency variables rather than 45 367, so §4.3's "IDF is unavailable"
is true of the coupling *as parameterised there* and false in a basis matched to
its smoothness: the architecture was selected by a representation choice rather
than by the physics. The measurements both support and constrain that. Supporting:
the fit is a contraction onto a single design, so the sub-problem has an
essentially unique solution, which is what a master/sub split needs. Constraining:
past about 2 mm of added harmonic content there is no linkage to be found, so a
master would need the reachable set as an explicit trust region. That is a study
of its own.
