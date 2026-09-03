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
above 1000 rpm. A well-conditioned slider-crank does the opposite, its peak
bearing load falling with speed by inertia relief. Same physics, opposite sign.

**A relaxation made for numerical reasons is a promise about tolerance.** The
equalities make the feasible set measure zero, forcing a relaxation into bands;
those bands are then only 1.7 standard deviations wide against the scatter of
the parts. Treating that probabilistically is not an embellishment — it is the
only way to know whether the relaxed requirement is met. It is not: 64.5 %
chance of missing at least one requirement, and a gap bound no ISO grade holds.
Worse, most of that is self-inflicted -- a deterministic optimizer converges
onto its active constraints, and designs beside this one halve the probability
at no cost in range (§6.2).

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
| solid round bars for every member | mass falls perhaps 30 % with tubes; changes the mass/inertia balance driving §6.1 |
| Coulomb friction with constant coefficients | absolute FMEP uncertain by ~30 %; rankings robust, since comparisons are at equal coefficients |
| instantaneous combustion, no heat transfer | indicated efficiency optimistic by several points, equally for both mechanisms |
| **no gas exchange** | optimistic for both, but **not equally** — see below; it flatters the conventional engine and understates §6.3 |
| reliability compared across mechanisms of different dimensionality | the slider-crank's two toleranced lengths against the EX-link's eleven is a real difference, not an artefact, but it means §6.3's reliability columns are not like-for-like in the way its range columns are |
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
against the EX-link by some margin** -- as, separately, is the fact that its
headline figure lets the baseline violate two limits the EX-link is held to. How large a margin is not established
here: a real engine recovers a fraction of the theoretical maximum, and that
fraction depends on valve timing and port design the model does not represent.
What can be said is the sign, and that the sign runs the opposite way to the
firing-frequency effect §6.3 spends most of its length removing.

Modelling it properly needs a valve-timing model and a pumping loop, which is a
larger change than any other item in this table.

### Method

**The reliability estimator is first order.** FORM under-predicts the most
nonlinear constraint, 0.42 against 0.54 sampled. Its derivatives are finite
differences because the exact route needs $\nabla^2 g$, which the package does
not compute.

**The mixed-integer bound is not a bound.** Outer approximation's guarantee
requires a convex sub-problem, which this problem violates comprehensively.

**The global optimum is not established.** Uniform multistart is inapplicable
(§3.4); manifold-projected restarts showed the *efficiency* optimum was local
but reached feasibility in 0 of 6 attempts on the *range* problem at an
affordable budget. That question is open, not answered.

**Single operating point.** Everything is reported at one speed with a sweep
around it. A drive cycle would test whether the optimum is an artefact of the
point chosen.

**Single mechanism family.** Two topologies establish a contrast; three would
establish a trend.

### Scope

The specification cannot be met as written. Relaxing the gap to its required
0.054 mm moves the binding constraint to the stroke band, which itself needs
widening to $\pm 0.09$ mm. Whether those bounds are acceptable is a question for
the customer, not the optimizer.

## 7.3 Possible improvements

In rough order of value per unit of effort:

1. **A drive cycle** in place of the single operating point. Cheap, and it tests
   the result most likely to be point-specific.
2. **Tubular sections.** The largest single modelling conservatism, and it
   interacts directly with §6.1.
3. **A widened uncertainty model**, carrying material, load and friction
   scatter alongside the dimensional tolerances. This is the prerequisite for
   everything else on the reliability side: $\Sigma$ currently holds only ISO 286
   dimensional errors, which is why only seven of the twelve constraints can
   honestly carry a probability (§3.10). Widening it is what would let the
   bearing, saturation and vehicle constraints join.
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
5. **Sampling-based reliability as an outer check.** The first-order constraint
   is what is affordable per iteration; `gemseo-umdo`'s `Probability` statistic
   would bound the error the linearisation makes.
6. **A feasibility-restoration phase before each restart**, which is what would
   make the multistart of §3.9 conclusive on the range problem.
7. **Second-order derivatives of the constraints**, which would make the
   reliability margin exactly differentiable and remove the one place where
   finite differences enter a tight constraint.
8. **A third mechanism topology**, to turn the contrast of §6.3 into a trend.
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

**Whatever is left out of the problem is what the solve violates.** Four times
a limitation was reported here — a restoration phase is needed, the route is
limited by reachability, none of these designs is feasible, the fit cannot give
diverse starts — and four times it was a constraint omitted from the
formulation. Adding it fixed the problem each time. Against a target the
mechanism cannot reach, the unconstrained fit leaves the geometric set by ten
to fifteen units while the constrained one holds every constraint at its
boundary.

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
| optimised as a conventional engine (its own limits) | 2888 km/L, **+15.6 %** |
| the same, with its firing-frequency advantage removed | **−4.3 %** |
| held to the EX-link's own rod-angle and side-load limits | 2372 km/L, **+43 %** |
| reliability index at IT8, EX-link vs that baseline | 3.00 vs **8.22** |

**What the specification costs**

| | |
|---|---|
| probability the reference design misses a requirement | 64.5 % |
| the same for the best design sampled beside it | 30.8 %, at +0.10 % range |
| gap bound a $10^{-3}$ target needs | 0.054 mm against 0.01 specified |
| stroke band the *system* then needs | ±0.15 mm against ±0.05 |

---

Next: [Running the code](implementation.md)
