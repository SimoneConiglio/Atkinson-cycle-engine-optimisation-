# 7. Conclusions

## 7.1 What the study establishes

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
| instantaneous combustion, no heat transfer or gas exchange | indicated efficiency optimistic by several points, equally for both mechanisms |
| constant crankshaft speed | the flywheel sizing already prices the fluctuation this assumes away |
| pin-jointed trigonal link | small; it is a stiff triangle either way |

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
4. **Closing the reliability loop**, by attaching
   {py:class}`exlink.robustness.FailureProbabilityDiscipline` as a constraint
   rather than reporting $P_f$ after the fact. The discipline exists and is
   tested; no scenario uses it, so every design in §6 was obtained
   deterministically and audited afterwards.
5. **Sampling-based reliability as an outer check.** The first-order constraint
   is what is affordable per iteration; `gemseo-umdo`'s `Probability` statistic
   would bound the error the linearisation makes.
6. **A feasibility-restoration phase before each restart**, which is what would
   make the multistart of §3.9 conclusive on the range problem.
7. **Second-order derivatives of the constraints**, which would make the
   reliability margin exactly differentiable and remove the one place where
   finite differences enter a tight constraint.
8. **A third mechanism topology**, to turn the contrast of §6.3 into a trend.
9. **A prescribed-motion sub-problem** as a feasible-point generator for the
   multistart of item 6, and — more speculatively — as the basis of a
   functional IDF. §7.4 sets both out, with the measurement that makes the
   second plausible.

## 7.4 A prescribed-motion formulation, and what it would change

The formulation of §3 asks for the best mechanism and lets the piston motion
$\lambda(\theta)$ fall out of it. Kinematic synthesis conventionally does the
opposite: it *prescribes* a target motion $\lambda^{\star}(\theta)$ and fits the
linkage to it,

$$\min_X\; J(X) = \int_0^{2\pi} \bigl(\lambda(\theta; X) - \lambda^{\star}(\theta)\bigr)^2
   \, w(\theta)\, d\theta ,$$

which is the same move as minimising $\lVert u - u^{\star}\rVert^2$ instead of
compliance in compliant-mechanism topology optimization. It is worth setting
out what that would and would not buy here, because the answer is not uniform.

### Where it is genuinely simpler

**It dissolves the measure-zero feasible set.** Both equalities are functionals
of $\lambda$ alone: $STE$ is the span of the expansion stroke and $\varepsilon$
is fixed by the extreme volumes. A target built to have the right stroke and
compression ratio satisfies them *by construction*, so they leave the
constraint set and reappear as part of the residual. The feasible set becomes a
full-dimensional box — which restores, at a stroke, every method §2.6 and §2.9
had to rule out: sampling, population methods and uniform multistart all become
admissible, because there is no longer a manifold to miss.

**The objective becomes a sum of squares**, so Gauss–Newton applies and
$J^{\mathsf{T}}J$ supplies a positive-semidefinite Hessian approximation for
free — materially better conditioned than the current problem near the
singularity.

**The derivatives get easier.** The residual is evaluated on a *fixed*
$\theta$-grid, so it contains no maximum over the revolution and needs no
envelope-theorem argument (§3.5). The existing kinematic Jacobian is the whole
derivative.

### Where it is not

**It presupposes its own answer.** In four-bar path synthesis $\lambda^{\star}$
comes from the application — trace this curve, dwell here. Here the best
piston motion is precisely what the study set out to find, so prescribing one
replaces an optimization over performance with an interpolation onto a guess.

**It reintroduces the proxy problem.** A least-squares residual has no exchange
rate with mass, friction or range — the same objection §2.3 raises against
$\eta$. Two mechanisms with equal residual can differ substantially in bearing
load.

**The residual floor is uninterpretable.** Eleven variables against several
hundred residuals is over-determined, so the fit is generically imperfect and
$J^{\star} > 0$ measures the distance from $\lambda^{\star}$ to the reachable
set, which is a property of the target rather than of the design.

### How it adapts to the multidisciplinary problem

Two uses, of quite different depth.

**(a) As a feasible-point generator.** This is the missing piece of §7.2's
"uniform multistart is inapplicable". Fitting to a $\lambda^{\star}$ constructed
with $STE = 74$ and $\varepsilon = 16$ lands *on* the equality manifold by
construction, and is an unconstrained box problem that can be run from
thousands of random starts cheaply. Each converged fit is then a feasible start
for the real range optimization. The target is used as a sampling device, not
as an objective, so none of the objections above applies — and it would turn
the multistart question from unanswerable into merely expensive.

**(b) As a functional decomposition — which would overturn §3.6.** Note what
$\lambda(\theta)$ is in the dependency graph: the *only* quantity the linkage
sends downstream. The cycle takes $\lambda \to V \to p$; the dynamics take
$\ddot\lambda$; the vehicle takes the resulting work. The eleven dimensions
affect nothing except through $\lambda$. So $\lambda(\theta)$ **is** the coupling
variable of the whole problem, and IDF on it reads:

- *master*: choose $\lambda^{\star}$ in a finite basis to maximise range, with no
  linkage in the loop at all;
- *sub*: given $\lambda^{\star}$, fit the linkage — classical dimensional
  synthesis, box-constrained least squares;
- *consistency*: drive $\lVert \lambda(X) - \lambda^{\star} \rVert$ to zero.

§3.6 rejected IDF because the coupling has 45 367 components. But that count is
*pointwise*, and $\lambda$ is smooth and periodic. Expanded in a Fourier basis
it is not high-dimensional at all:

| design | harmonics for RMS $< 0.1$ mm | $< 0.01$ mm | $< 0.001$ mm |
|---|---|---|---|
| refined (geometric, AL) | 12 | 18 | 25 |
| gradient (geometric, SLSQP) | 15 | 23 | 30 |
| coupled (minimum mass) | 10 | 14 | 19 |
| range (vehicle-level) | 10 | 16 | 20 |

Fourteen to twenty-three complex coefficients reproduce the piston motion to
0.01 mm RMS — below the 0.020 mm machining scatter §6.2 puts on $STE$, so
tighter than the part can be made. A functional IDF would therefore carry of
order **30 to 50** consistency variables, not 45 367: three orders of magnitude
less, and entirely tractable.

That is a concrete correction to a stated conclusion of this study. §3.6's
"IDF is unavailable" is true of the coupling *as parameterised there*, and
false of the same coupling in a basis matched to its smoothness. The general
lesson is the one §2.10 claims for the feasible set, in a second instance: the
architecture was selected by a representation choice, not by the physics.

The obvious hazard is attainability. Not every periodic function is the piston
motion of a seven-bar linkage, so a master free to ask for any $\lambda^{\star}$
will ask for unfittable ones and the consistency residual will not close — the
same difficulty as prescribing an unreachable displacement field in topology
optimization. The remedies are the usual ones: keep the basis small, impose the
known necessary conditions ($\lambda > 0$, two top and two bottom dead centres
per revolution, the $STE$ and $\varepsilon$ identities), and trust-region the
master on the residual the sub-problem actually achieved.

Of the two, (a) is cheap and should be done; (b) is the interesting one and
would need its own study.

## 7.5 Headline numbers

| | |
|---|---|
| range of the coupled reference design | 3338 km/L |
| engine mass at that design | 12.2 kg |
| range of an optimised conventional engine, identical code | 2888 km/L |
| advantage over it | +15.6 % |
| the same with the firing-frequency advantage removed | **-4.3 %** |
| probability the reference design misses a requirement | 64.5 % |
| gap bound needed for a $10^{-3}$ target | 0.054 mm against 0.01 specified |

---

Next: [Running the code](implementation.md)
