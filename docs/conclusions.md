# 7. Conclusions

## 7.1 What the study establishes

**Imposing a constraint and checking it are different searches.** The same
SLSQP, on the same problem, reaches 4.9 % further when the coupled and vehicle
constraints are held during the search rather than verified at the end (§6.4).
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
against the EX-link by some margin**. How large a margin is not established
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
9. **Converging the imposed form of §3.10.** The formulation behind §6.4's 3501 km/L
   stopped at its iteration cap, not at a convergence test, so that figure is a
   lower bound on what it reaches. Running it to convergence, and from several
   starts, is the cheapest remaining gain in the study; the
   duplicate-evaluation fix has since halved what each iteration costs. The
   functional IDF §7.4 sets out is a larger question again and would need its
   own study.

## 7.4 A prescribed-motion formulation, measured

The formulation of §3 asks for the best mechanism and lets the piston motion
$\lambda(\theta)$ fall out of it. Kinematic synthesis conventionally does the
opposite: it *prescribes* a target motion $\lambda^{\star}(\theta)$ and fits the
linkage to it,

$$\min_X\; J(X) = \sum_k \bigl(\lambda_k(X) - \lambda^{\star}_k\bigr)^2 ,$$

which is the same move as minimising $\lVert u - u^{\star}\rVert^2$ instead of
compliance in compliant-mechanism topology optimization. This section was
written as a proposal and is now a measurement:
{py:mod}`exlink.synthesis` implements it, and the result is more interesting
than the proposal was, because the part that looked incidental turned out to
be the whole difficulty.

### What is simpler, and it is genuinely simpler

Both equalities are functionals of $\lambda$ alone: $STE$ is the span from top
dead centre to the deeper bottom dead centre, and $\varepsilon$ follows from the
shallower one through the swept volume. So a target built to have the right two
strokes satisfies both requirements before any linkage exists.
{py:func}`~exlink.synthesis.target_motion` solves two harmonics against the two
strokes, and {py:func}`~exlink.synthesis.describe_target` measures the result
with `find_phases` — the same code the constraints use, not the formula the
target was built from:

| | target | required |
|---|---|---|
| $STE$ | 74.0000000000 mm | 74 mm |
| $\varepsilon$ | 16.0000000000 | 16 |
| $g$ | 0.0 mm | — |

The objective also becomes a sum of squares, so Gauss–Newton applies and the
residual Jacobian supplies a positive-semidefinite Hessian approximation for
free; and the residual lives on a fixed $\theta$-grid, so it contains no maximum
over the revolution and needs no envelope-theorem argument (§3.5).

### The difficulty is attainability, and it is decisive

The proposal claimed the equalities would then be satisfied "by construction".
**That claim was wrong**, and the measurement says so. The fit is
over-determined — eleven variables against several hundred residuals — so
$\lambda(X) = \lambda^{\star}$ is unattainable, and what matters is whether the
residual is small enough to land inside the tolerance bands of §3.4. Fitting
the reference design to a two-harmonic target:

| target | correction | fit RMS | $\lvert\Delta STE\rvert$ | $\lvert\Delta\varepsilon\rvert$ | inside the bands |
|---|---|---|---|---|---|
| two harmonics, built from nothing | $A = 32.33$, $B = 9.02$ | 1.1621 mm | 0.0735 mm | 0.5594 | **no** |
| the same, seeded from a real motion | $A = 0.038$, $B = 0.042$ | **0.0011 mm** | 0.0013 mm | 0.0001 | **yes** |

A target can be exactly on the equality manifold and still be useless, because
being on the manifold is a property of the target while being *fittable* is a
property of the mechanism. The two-harmonic target misses the compression-ratio
band by a factor of eleven. In hindsight the Fourier table below predicted it:
a real motion needs ten to twenty-three harmonics, so a two-harmonic target sits
far outside the reachable set.

{py:func}`~exlink.synthesis.target_from_design` is the repair. It seeds the
target from a motion an actual design produces — reachable by definition — and
adds only the two harmonics that move the strokes onto the requirement. The
correction is three orders of magnitude smaller, and so is the residual.

### What it is good for: feasible points, not multistart

With a reachable target the generator does the one thing the main formulation
cannot. Uniform sampling of the design box finds a design on the equality
manifold in **0 of 12 000** draws (§3.4). Sampling and *then fitting* is a
different operation, and most fits land inside both bands:

| scatter about the reference | starts | fits obtained | inside the bands | distinct designs |
|---|---|---|---|---|
| 10 % | 30 | 22 | 22 | **1** |
| 25 % | 30 | 13 | 13 | **1** |
| 50 % | 30 | 9 | 8 | **1** |
| 80 % | 30 | 3 | 2 | **1** |

The last column is the finding. Every fit converges to the *same* linkage,
whatever start it is given. That is not a defect of the optimizer: the motion
very nearly determines the mechanism, so fitting to one target yields one
design. **Fitting to a single target solves the feasibility problem and not the
multistart problem**; diversity has to be asked for separately, and the rest of
this section is about how.

Diversity therefore has to come from varying the *target*. Perturbing it with
third-to-fifth-harmonic content and re-solving the two stroke harmonics — so the
target stays exactly on the manifold and only its shape changes — gives:

| perturbation | targets built | fits | inside the bands | distinct designs |
|---|---|---|---|---|
| 0.2 mm | 12 | 12 | 4 | **4** |
| 0.5 mm | 12 | 12 | 1 | **1** |
| 1.0 mm | 12 | 12 | 3 | **3** |
| 2.0 mm | 10 | 10 | 0 | 0 |

Distinct now equals in-band at every level: varying the target gives genuinely
different linkages where varying the start gave one. Past about 2 mm of added
harmonic content the target leaves the reachable set and no fit lands in band,
which is the attainability limit measured from the other side.

**Are these points usable, or only near the manifold?** Usable. Of 30 random
starts fitted to the fixed reachable target, 26 produce a fit, all 26 land
inside both bands, and **all 26 are feasible against the complete constraint
set** — every inequality, coupled and vehicle constraints included. The
generator does not merely approach the equality manifold; it lands on usable
designs.

That is a property of *this* target rather than of the method, and the
distinction matters for the jittered targets below, whose unconstrained fits
do sometimes leave the geometric set. The repair is not a separate restoration
phase: it is to stop discarding the inequalities. Absorbing the equalities into
the target removed the measure-zero part of the feasible set, so what remains is
an ordinary box-and-inequality problem that SQP handles directly,

$$\min_X \; \lVert \lambda(X) - \lambda^{\star} \rVert^2
  \quad \text{s.t.} \quad g(X) \le 0, \; X \in [X_{lb}, X_{ub}] ,$$

which is what {py:func}`~exlink.synthesis.fit_within_constraints` solves.
Dropping $g$ was never justified by the reformulation — only the *equalities*
were absorbed — and keeping it costs one SQP in place of one least-squares
solve. Measured, the difference is not marginal:

| target | fit | worst $g$ | inside the bands |
|---|---|---|---|
| reachable | unconstrained | $-0.0026$ | yes |
| reachable | constrained | $-0.0026$ | yes |
| jittered, 0.5 mm | unconstrained | $+10.0$ to $+15.7$ | mostly no |
| jittered, 0.5 mm | constrained | $\mathbf{+0.0000}$ | no |

On a reachable target the two agree exactly, so the constraints cost nothing.
On an unreachable one the unconstrained fit violates by *ten to fifteen units*
while chasing a motion no buildable linkage produces, and the constrained fit
holds every constraint at its boundary instead. Discarding $g$ was not a
simplification, it was a way of not noticing.

The last row also shows the constrained fits still outside the stroke bands,
because holding $g \le 0$ while chasing an unreachable motion trades one
against the other. But the bands are constraints too, and the same argument
puts them in the problem rather than in a check afterwards
({py:func}`~exlink.synthesis.fit_within_constraints` with ``hold_bands``).
Doing so on six jittered targets:

| | |
|---|---|
| constraints held ($\max g$) | $+0.0000$ on all six |
| inside both bands | six of six, at the boundary |
| feasible against the *whole* set | **two** |
| distinct designs among those | **two** |

So the diversity route does work, and the earlier reading that it is limited by
reachability was premature. Reachability limits how closely the motion can be
tracked; it does not prevent the solve from returning usable designs, because
the constrained problem simply trades motion accuracy for feasibility and still
lands somewhere buildable.

The four that remain infeasible are instructive in the same way as everything
above: they satisfy every constraint that was *given* to them — $\max g = 0$,
both bands held — and fail the coupled and vehicle constraints, which are not
in the fit. The pattern has now repeated three times. Each apparent limitation
of the reformulation turned out to be a constraint that had been left out of
it, and adding the constraint fixed it. Those two are the next to add; they are
excluded here only because each evaluation would then carry an MDA.

One caution that does not dissolve this way. Every one of these fits sits at
$g = 0$ exactly — SQP converging onto its active constraints — and §6.2 shows
that sitting on active constraints is precisely what makes a design unreliable.
A generator built this way should hand its output to a reliability-aware stage
rather than be trusted as a final design.

### Where it ended up

The argument that took the inequalities and then the bands into the fit does not
stop there, and its endpoint is not a synthesis problem at all: if every
constraint belongs in the problem, and the fit is only a means of reaching a
good design, then the objective should be the one the application scores. That
formulation is the second form in §3.10, its result is §6.4, and it is the best design this study
produces. What remains here is what the prescribed motion taught on the way,
which is not the same thing as what it achieved.

**The lesson that generalises.** Four times in this exercise a limitation was
reported — a restoration phase is needed, the route is limited by reachability,
none of these designs is feasible, the fit cannot give diverse starts — and four
times it turned out to be a constraint that had been left out of the problem.
Adding the constraint fixed it each time. The habit worth keeping is to ask,
of any apparent limitation of a reformulation, whether it is a property of the
problem or of what was omitted from it; here it was never the former.

**What the target is, precisely.** It is a fallback objective, not a
constraint, and §6.4 measures both halves of that claim: it fires for a quarter
of the search from a start where the range does not exist, and not at all from
one where it does, ending 5.53 mm from a motion it was free to abandon. A
constant penalty cannot do the first — there is no gradient in a constant — and
a constraint would have prevented the second.

### The functional decomposition this suggests

$\lambda(\theta)$ is the *only* quantity the linkage sends downstream: the cycle
takes $\lambda \to V \to p$, the dynamics take $\ddot\lambda$, the vehicle takes
the resulting work. So $\lambda$ **is** the coupling variable, and IDF on it
reads: a master choosing $\lambda^{\star}$ in a finite basis to maximise range,
a sub-problem fitting the linkage to it, and
$\lVert \lambda(X) - \lambda^{\star} \rVert$ as the consistency constraint.

§3.6 rejected IDF because the coupling has 45 367 components — but that count is
*pointwise*, and $\lambda$ is smooth ({py:func}`exlink.formulations.motion_harmonics`):

| design | harmonics for RMS $< 0.1$ mm | $< 0.01$ mm | $< 0.001$ mm |
|---|---|---|---|
| refined (geometric, AL) | 12 | 18 | 25 |
| gradient (geometric, SLSQP) | 15 | 23 | 30 |
| coupled (minimum mass) | 10 | 14 | 19 |
| range (vehicle-level) | 10 | 16 | 20 |

Fourteen to twenty-three coefficients reproduce the motion to 0.01 mm RMS,
tighter than the 0.020 mm machining scatter §6.2 puts on $STE$. A functional IDF
would carry of order **30 to 50** consistency variables, not 45 367 — three
orders of magnitude fewer. §3.6's "IDF is unavailable" is true of the coupling
*as parameterised there* and false of the same coupling in a basis matched to
its smoothness: the architecture was selected by a representation choice, not by
the physics.

The measurements above both support and constrain that. Supporting it: the fit
is a contraction onto a single design, so the sub-problem has an essentially
unique solution — precisely what a well-posed master/sub split needs.
Constraining it: the master cannot be allowed to ask for any $\lambda^{\star}$ in
the basis, because past ~2 mm of harmonic content there is no linkage to be
found and the consistency residual will not close. A master would need the
reachable set as an explicit trust region, calibrated on the sub-problem's
achieved residual. That is the shape of the study this suggests, and it is
not a small one.

## 7.5 Headline numbers

| | |
|---|---|
| best strictly feasible design | 3338 km/L, 12.2 kg |
| best design overall (§6.4), $2\times10^{-4}$ outside one band | **3501 km/L** |
| what imposing every constraint during the search is worth | **+4.9 %** |
| range of an optimised conventional engine, identical code | 2888 km/L |
| advantage over it | +15.6 % |
| the same with the firing-frequency advantage removed | **-4.3 %** |
| probability the reference design misses a requirement | 64.5 % |
| the same, for the best sampled design beside it | **30.8 %**, at +0.10 % range |
| gap bound needed for a $10^{-3}$ target | 0.054 mm against 0.01 specified |

---

Next: [Running the code](implementation.md)
