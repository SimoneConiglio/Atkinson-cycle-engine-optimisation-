# 7. Conclusions

## 7.1 What the study establishes

**Imposing a constraint and checking it are different searches.** The same
SLSQP, on the same problem, reaches 4.9 % further when the coupled and vehicle
constraints are held during the search rather than verified at the end, and
1.7 % further again when a reliability target is held too, which is a different
and better design rather than a smaller number (§6.4).
Nothing about the algorithm changed; the problem was posed better.

**The objective matters more than the algorithm.** The conventional formulation
prices nothing, its central quantity is not an efficiency, and it cannot see the
parts. Replacing it with range changes the answer qualitatively — not by a few
per cent, but from *the singularity is optimal* to *the singularity is the worst
place to be*.

**Conditioning decides the sign of the inertia effect.** A quasi-statically
optimised linkage drifts to its transmission-angle singularity, where the lever
arm is longest and the accelerations are worst; it has no feasible structure
above 2000 rpm. A well-conditioned slider-crank does the opposite, its peak
bearing load falling with speed by inertia relief. Same physics, opposite sign.

**A relaxation made for numerical reasons is a promise about tolerance.** The
equalities make the feasible set measure zero, forcing a relaxation into bands;
those bands are then only 1.7 standard deviations wide against the scatter of
the parts, and the reference design has a 66.4 % chance of missing at least one
requirement. Which bounds are responsible is not visible in any nominal
quantity: the top-dead-centre bound is set finer than the model's own
resolution, and once it is widened the band on the expansion stroke governs
everything, over four orders of magnitude of failure probability. Widening both
costs 0.47 % of range. Most of the rest is self-inflicted — a deterministic
optimizer converges onto its active constraints, and designs beside this one
halve the probability at no cost in range (§6.2).

**The topology is worth 17.6 %, and all three terms of the objective
contribute.** Against a conventional engine sized by identical models and
optimised over its own degrees of freedom rather than proportioned by hand, the
linkage reaches 3395 km/L against 2888. Both take 720° of their crankshaft per
cycle, so the comparison is at equal speed and equal firing rate with nothing to
correct. Extended expansion — the feature the topology exists for — is the
smallest of the three contributions, five per cent of indicated efficiency; the
larger two are a lower side load and a lighter flywheel (§6.3).

An earlier version of this paragraph attributed those two to "having eleven
dimensions to place rather than two", and §6.9 shows that reading is wrong by a
factor of four. Handing the conventional engine a third freedom buys 0.51 %,
which extrapolates to 4.6 % over nine, not 17.6 %. The dimensions are what let
the linkage exploit its topology; they are not what give it the advantage.

**The largest stated modelling conservatism is worth 0.9 %, and only above the
design speed.** Solid round bars looked like the study's weakest assumption, and
the members turn out to be 1.3 % of the engine. Boring them buys nothing at all
below 2000 rpm and 31 % of range at 3200, because what a bore takes out is not
weight but inertia in the load path, which is the quantity §6.1's sizing fixed
point is made of. Given to both engines the comparison of §6.3 moves from
+17.6 % to +17.0 %, and §6.1's divergence still bites — three times later, and
just as hard (§6.6). In a coupled problem the size of an effect is not the size
of the thing it acts on.

**A reliability model reports silence as safety.** Widening the uncertain
vector from the eleven dimensions to seventeen — material strength, stiffness,
density, friction, gas load — leaves every geometric constraint identical to
the last figure and reveals that the binding constraint of the whole design was
never in the model: the gear pair the mixed-integer master chose sits exactly
on its face-width limit, at $\beta = 0.00$. The honest note in §3.10 that only
some constraints could carry a probability was a warning about the answer, not
a footnote to it, and the pair that fixes it costs nothing at all (§6.7).

**A schedule of speeds rewards smoothness, and the advantage widens to 19.3 %.**
Scored over four speeds as one engine rather than at one point, both designs
lose range and the conventional engine loses more. The reason is structural: a
wider speed range makes the flywheel a larger share of the mass — 84 % of the
linkage's and 96 % of the baseline's — and the flywheel is the item the flat
torque curve wins by a factor of 1.79. The slider-crank is the lighter
*mechanism* and the heavier *engine*, and asking it to cover a range rather
than a point makes that worse (§6.5).

**Decomposition buys structure, not speed.** Bi-level outer approximation halves
the sub-solves against enumeration and lands 0.6 % short, on a bound that is not
valid because the sub-problem is nonconvex. What it buys is a mixed-integer
statement, principled handling of infeasible lattice points, and a stopping
criterion in place of a guessed budget.

## 7.2 Limitations of the framework

Grouped by what would have to change to remove them.

### Modelling

| limitation | effect if relaxed |
|---|---|
| ~~solid round bars for every member~~ | **measured** (§6.6): the members are 1.3 % of the engine, so a bore is worth 0.9 % of range and only above the design speed, where it acts on the inertia in the load path rather than on weight. It moves where §6.1's divergence bites, not whether it happens. |
| Coulomb friction with constant coefficients | absolute FMEP uncertain by ~30 %; rankings robust, since comparisons are at equal coefficients |
| instantaneous combustion, no heat transfer | indicated efficiency optimistic by several points, equally for both mechanisms |
| **no gas exchange** | optimistic for both, but **not equally** — see below; it flatters the conventional engine and understates §6.3 |
| reliability compared across mechanisms of different dimensionality | the slider-crank's two toleranced lengths against the EX-link's eleven is a real difference, not an artefact, but it means §6.2's reliability figures are not like-for-like in the way §6.3's range figures are |
| constant crankshaft speed | the flywheel sizing already prices the fluctuation this assumes away |
| pin-jointed trigonal link | small; it is a stiff triangle either way |

#### Neglecting gas exchange is not a neutral simplification

Both cycles hold the intake and exhaust strokes at plenum pressure, so the
gauge pressure is zero across them and **neither engine pays any pumping
work**. Stated that way it sounds even-handed. It is not, and the direction
matters for §6.3.

An over-expanded engine opens its exhaust valve later, on a charge that has
expanded further and therefore sits closer to ambient. Measured on the two
mechanisms this study compares, at the same compression ratio:

| | expansion ratio | $p$ at exhaust-valve-open | $p/p_0$ |
|---|---|---|---|
| EX-link | 20.8 | 0.1475 | **1.23** |
| slider-crank | 16.0 | 0.2040 | **1.70** |

The work still recoverable from that charge — expanding it isentropically to
ambient, which is the theoretical maximum and not an achievable figure — is
9.5 % of indicated work for the EX-link and 26.4 % for the slider-crank. The
model discards both.

So the simplification removes a loss that is roughly two and a half times
larger for the conventional engine, and **§6.3's comparison is conservative
against the EX-link by some margin** — as, separately, is the fact that it lets
the baseline violate two limits the linkage is held to. How large a margin is
not established here: a real engine recovers a fraction of the theoretical maximum, and that
fraction depends on valve timing and port design the model does not represent.
What can be said is the sign: it runs in the linkage's favour, so §6.3's
17.6 % is a lower bound on that account.

Modelling it properly needs a valve-timing model and a pumping loop, which is a
larger change than any other item in this table.

### Method

**The reliability estimator is first order, and first order is enough here —
once it is applied to the right function.** §6.8 found FORM optimistic by a
factor of seven at the study's result and traced it to a maximum being
linearised at a tie, not to curvature; with both branches carried, FORM agrees
with 150 000 sampled builds to 1 %. What is *not* settled is that the same trap
is absent elsewhere: the check has been run at two designs, and there is no
automatic warning when a constraint is a maximum evaluated near its tie.

**The uncertain vector was too narrow, and now is not.** §6.7 widens it from
eleven dimensions to seventeen and prices every constraint. What is left is
that the six added parameters are taken independent, and that the widened
model costs eighteen analyses rather than one Jacobian, so it remains a
reporting tool rather than something the optimizer can call.

**The mixed-integer bound is not a bound.** Outer approximation's guarantee
requires a convex sub-problem, which this problem violates comprehensively.

**The global optimum is not established.** Uniform multistart is inapplicable
(§3.4); manifold-projected restarts showed the *efficiency* optimum was local
but reached feasibility in 0 of 6 attempts on the *range* problem at an
affordable budget. §6.10 removes the obstacle — restoration takes the same six
starts to 6 of 6 — without answering the question, which needs a range solve
from each. It is open, and now cheaply askable.

**The schedule is assumed, not measured.** §6.5 answers the single-point
objection — the design is scored over a four-point schedule as one engine and
keeps its advantage — but the distance weights are stated rather than derived,
because no survey of the track exists. A different spread gives a different
cycle range; what it does not plausibly change is the ranking, since the
EX-link leads at every point of the one tested.

**Single mechanism family.** Two topologies establish a contrast; three would
establish a trend. §6.9 adds a third point on the dimensionality axis rather
than a third topology, which settles that the trend does not simply run with
the variable count but leaves the topological question open.

### Scope

The results at $\beta = 3$ are stated at a widened specification: the
top-dead-centre gap at 0.1 mm and both equality bands at $\pm 0.15$. §6.2 prices
that widening at 0.47 % of range and shows what it buys, but whether those
bounds are acceptable is a question for the customer, not the optimizer. At the
bounds as written the mechanism reaches no reliable design at all.

## 7.3 Possible improvements

In rough order of value per unit of effort:

1. ~~**A drive cycle** in place of the single operating point.~~ **Done**
   (§6.5). The design holds across a four-point schedule scored as one engine —
   worst point at 0.91 of the best — and the advantage over the conventional
   engine widens from +17.6 % to +19.3 %, because a schedule makes the flywheel
   a larger share of both engines and the flywheel is where the flat torque
   curve pays. What remains open is re-optimising the linkage *for* a schedule
   rather than scoring it over one, which would make +19.3 % a lower bound
   rather than the answer.
2. ~~**Tubular sections.**~~ **Done** (§6.6), and the answer was not the one
   the entry assumed. The members are 1.3 % of the engine, so boring them is
   worth 0.9 % of range — but all of it above the design speed, because what a
   bore removes is inertia in the load path, not weight. The optimum operating
   speed moves up 300 rpm and the engine falls from 12.9 kg to 10.4 kg. §6.1
   survives it: the quasi-static optimum is still unbuildable at 4000 rpm with
   tubes. What is still missing is a wall-buckling check, end-fitting mass, and
   a bore ratio treated as a design variable rather than scanned — all three of
   which work against the tube, so 0.9 % is an upper bound.
3. ~~**A widened uncertainty model**, carrying material, load and friction
   scatter alongside the dimensional tolerances.~~ **Done** (§6.7), and it
   found the design's binding constraint. Widening $\Sigma$ to seventeen
   entries leaves all eight geometric constraints identical to the last figure
   — §6.2 and §6.4 stand — and prices the five that had no probability at all.
   One of them, the gear pair's face width, sits exactly on its limit at
   $\beta = 0.00$, so the system probability is $5\times10^{-1}$ rather than
   $1.3\times10^{-3}$. The fix costs nothing: the pair the exhaustive search
   preferred is 0.4 km/L better and puts that constraint at $\beta = 3.97$.
   What remains is a correlated model of the two strengths (independence is
   conservative here) and second-order treatment of ``saturation``, which is a
   threshold on a fixed point and the least linear constraint in the set.
4. **A way to reach the reliable region**, which is now known to be worth
   reaching. The loop itself is closed --
   ``build_range_scenario(beta_target=...)`` constrains the reliability index
   -- but SLSQP from the deterministic optimum does not move, while sampling
   beside it finds fully feasible designs that halve $P_f$ at no cost in range
   (§3.10, §6.2). The deterministic optimum is dominated, so this is the
   cheapest improvement on the list in engineering terms and the most clearly
   algorithmic in nature: what is missing is not the constraint, nor a
   trade-off to negotiate, but a search able to cross a thin feasible region --
   a restoration phase, a continuation in $\beta$ from a sampled start, or the
   prescribed-motion generator of §7.4 supplying starts already on the
   manifold.
5. ~~**Sampling-based reliability as an outer check.**~~ **Done** (§6.8), and
   it did not confirm the first-order estimate — it overturned it. 150 000
   exact builds put the study's result at $9.3\times10^{-3}$ against FORM's
   $1.3\times10^{-3}$. The cause is that the expansion stroke is a *maximum*
   over two top dead centres, and the optimizer had driven those two to within
   0.107 μm of each other — a fifth of a wavelength of light, and
   seventy-five times finer than the scatter of the parts — so the
   linearisation used the branch that
   attained the maximum while the parts breached the other. Carrying both
   branches — one extra row of a Jacobian already being computed — brings FORM
   to within 1 % of sampling at every design tested. What remains is to run the
   check at more designs, and to automate it: nothing in a FORM output signals
   that it is being evaluated at a kink.
6. ~~**A feasibility-restoration phase before each restart**~~ **Done**
   (§6.10), for the phase itself: six starts drawn as §6.11's were, 0 of 6
   feasible as drawn and **6 of 6** after, at about 290 evaluations each. The
   epigraph form — maximise $t$ subject to $c_i \ge t$ — works where the range
   solve does not because it has no objective to trade against the constraints,
   so its QP always has an admissible step; starts as far as 14.9 outside
   restore no more slowly than starts 1.5 outside. What this does *not* settle
   is the question behind the entry: six restored starts make the multistart
   answerable, and answering it needs six range solves from them, at six times
   the cost of §6.4's own run. Whether that optimum is global is still open.
7. ~~**Second-order derivatives of the constraints.**~~ **Attempted, and the
   wrong instrument** (§6.8). Differencing the analytic gradient gives a matrix
   whose norm scales exactly as $1/h$ over two decades of step and is 100 %
   asymmetric at every step — the signature of a discontinuous gradient rather
   than of a second derivative. The constraint surface at the study's result is
   not curved but *kinked*, and no order of Taylor expansion repairs a kink.
   The fix was a first derivative of the other branch. Second derivatives may
   still be worth having for the exact $\partial\beta/\partial x$ of §3.10's
   steered quantity, but they must be taken away from the tie, and that is a
   different piece of work from the one this entry asked for.
8. ~~**A third mechanism topology**, to turn the contrast of §6.3 into a
   trend.~~ **Partly done** (§6.9), and it undercut the reason the entry was
   written. A third *topology* is still missing — a genuine one needs its own
   internal 2:1 ratio, which is what unequal strokes require — but a third
   point on the *dimensionality* axis was cheap, and it shows the trend cannot
   be assumed to run with the number of design variables. A wrist-pin offset
   gives the conventional engine a third freedom and buys 0.51 %; extrapolated,
   nine such freedoms would buy 4.6 % against the linkage's 17.6 %. The offset
   improves the one term its topology leaves open — mechanical efficiency,
   0.787 to 0.792 — and cannot touch the other two at any number of
   dimensions.
9. **Converging §6.4.** Both solves there stopped at their iteration cap, not
   at a convergence test, so 3395 km/L and 3501 km/L are lower bounds on what
   the formulation reaches. Running them to convergence, and from several
   starts, is the cheapest remaining gain in the study. The functional IDF §7.4
   sets out is a larger question again and would need its own study.

## 7.4 What the prescribed motion taught

The formulation is §3.10's and its result is §6.4. What is left to record is
what the detour established, because most of it generalises past this problem.

**A target on the manifold need not be reachable.** Both equalities are
functionals of $\lambda$ alone, so a motion with the right two strokes satisfies
them exactly before any linkage exists — measured at 74.0000000000 mm and
16.0000000000. But being on the manifold is a property of the *target* and
being fittable is a property of the *mechanism*. A two-harmonic target is
exact and unreachable: the closest a seven-bar gets is 1.16 mm RMS, which
carries the fitted design outside both tolerance bands. Seeding the target from
a motion a real design produces
({py:func}`~exlink.synthesis.target_from_design`) recovers three orders of
magnitude of residual. The Fourier table below predicted this and was not read
that way.

**Fitting reaches the manifold; it does not give multistart.** Uniform sampling
finds a design on the equality manifold in 0 of 12 000 draws. Sampling and
*then fitting* reaches it in 22 of 30, and those fits are feasible against the
whole constraint set. But every fit converges to the same linkage whatever
start it is given — the motion very nearly determines the mechanism — so
diversity has to come from varying the target, not the start.

**What a solve violates is what the formulation left out.** Four apparent
limitations of the prescribed-motion route — that it needed a feasibility
restoration phase, that it was bounded by reachability, that its designs were
infeasible, that its fits could not supply diverse starts — each turned out to
be a constraint absent from the sub-problem rather than a property of the
method, and each was removed by restoring that constraint. The clearest case is
measurable: against a target the mechanism cannot reach, the unconstrained fit
leaves the geometric set by ten to fifteen units, while the same fit with the
inequalities imposed holds every one of them at its boundary.

**A flat region defeats a gradient method, three times over.** The reliability
constraint stalls where $P_f$ saturates (§6.4); the range objective has nothing
to descend where $R$ does not exist; and the system index is pinned at
$-8.2095$ for every violating design. Each repair has the same shape: supply
something that still varies where the quantity of interest does not — the
motion residual as the objective's middle rung, and $\min_i \beta_i$ in place
of the system index.

### The decomposition this suggests

$\lambda(\theta)$ is the *only* quantity the linkage sends downstream, so it is
the coupling variable, and IDF on it would put a master choosing
$\lambda^{\star}$ against a sub-problem fitting the linkage to it. §3.6 rejected
IDF because the coupling has 45 367 components — but that count is *pointwise*,
and $\lambda$ is smooth ({py:func}`exlink.formulations.motion_harmonics`):

| design | harmonics for RMS $< 0.1$ mm | $< 0.01$ mm |
|---|---|---|
| coupled (minimum mass) | 10 | 14 |
| gradient (geometric, SLSQP) | 15 | 23 |

Fourteen to twenty-three coefficients reproduce the motion to 0.01 mm RMS,
tighter than the 0.020 mm machining scatter. A functional IDF would carry of
order 30 to 50 consistency variables, not 45 367 — so §3.6's "IDF is
unavailable" is true of the coupling *as parameterised there* and false in a
basis matched to its smoothness. The architecture was selected by a
representation choice, not by the physics.

The measurements both support and constrain that. Supporting: the fit is a
contraction onto a single design, so the sub-problem has an essentially unique
solution, which is what a master/sub split needs. Constraining: past about 2 mm
of added harmonic content there is no linkage to be found, so a master would
need the reachable set as an explicit trust region. That is a study of its own.

## 7.5 Headline numbers

Full provenance for every design is in §6.0.

**What the linkage achieves**

| | |
|---|---|
| best strictly feasible design, specification as written | 3338 km/L, 12.2 kg |
| best nominal design, all constraints imposed | 3501 km/L |
| best design that also holds $P_f \le 10^{-3}$, bounds relaxed | **3395 km/L** |
| what the reliability requirement costs | −3 % |

**Against a conventional engine**

| | |
|---|---|
| optimised as a conventional engine, both at 720° per cycle | 2888 km/L |
| against the study's result, 3395 km/L | **+17.6 %** |
| against `COUPLED_DESIGN`, 3338 km/L | +15.6 % |
| over a four-point schedule, one engine each | 2733 km/L vs 3260 km/L, **+19.3 %** |
| with tubular members on both, each re-optimised | 2929 km/L vs 3425 km/L, **+17.0 %** |
| indicated efficiency | 0.457 → 0.480 |
| mechanical efficiency | 0.787 → 0.865 |
| engine mass | 16.9 kg → 12.9 kg |
| reliability index at IT8, linkage vs baseline off its cap | 2.35 vs **8.2** |
| the same baseline given a third freedom (wrist offset) | 2903 km/L, +0.51 % |

**What the bounds cost**

| | |
|---|---|
| probability the reference design misses a requirement | 66.4 % |
| the same for the best design sampled beside it | 30.8 %, at +0.10 % range |
| gap bound above which the gap stops binding | 0.054 mm; 0.1 mm adopted |
| stroke band the *system* then needs | ±0.15 mm against ±0.05 |
| range given up by widening both | −0.47 % |
| the same design against all thirteen constraints | $5\times10^{-1}$, binding on the gear pair |
| with the gear pair the exhaustive search preferred | $1.4\times10^{-3}$, at +0.4 km/L |
| failure probability bought | 0.664 → $1.9\times10^{-5}$ |

---

Next: [Running the code](implementation.md)
