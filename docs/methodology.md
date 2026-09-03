# 3. Methodology

Each subsection states a limitation of the formulation as it stands, then
derives what resolves it. Nothing here depends on the implementation; §4
describes that separately.

## 3.1 Design variables and objective

Eleven continuous dimensions describe the linkage:

$$X = (a,\; c,\; I,\; x_b,\; y_b,\; x_1,\; e,\; q_1,\; q_2,\; \theta_f,\; \theta_r)^{\mathsf T}$$

§2.3 rejected the conventional objective. The objective used here is the
quantity the application scores,

$$\max_X \; R(X) \quad [\mathrm{km/L}]$$

which requires the chain

$$X \;\longrightarrow\; \lambda(\theta_1) \;\longrightarrow\; p(V)
  \;\longrightarrow\; \text{loads} \;\longleftrightarrow\; \text{sections}
  \;\longrightarrow\;
  \begin{cases} W_{\text{brake}} \\ m_{\text{engine}} \end{cases}
  \;\longrightarrow\; R$$

Range prices the competing quantities at rates the physics fixes rather than the
designer: a point of brake efficiency is worth a fixed distance through the fuel
burnt per unit work; a millimetre of envelope is worth a fixed mass of
crankcase, and a kilogram a fixed distance through rolling resistance; a
newton-millimetre of torque ripple is worth a fixed mass of flywheel.

## 3.2 Evaluating the objective

The objective of §3.1 is not a formula but a chain. Each link is a first-order
model; each is checked in §4.4 against a result computed independently of it.

### 3.2.1 Motion and cycle

The loop-closure equations are inverted analytically rather than solved with
Newton–Raphson, which matters twice: the evaluation is fast, and the two
arccosine arguments it exposes become the well-posedness conditions

$$\delta_{c1} \le 1, \qquad \delta_{c2} \le 1, \qquad W = \max(\delta_{c1}, \delta_{c2})$$

that keep the linkage assemblable. $W \to 1$ is the transmission-angle
singularity that §6.1 turns on.

Before pressure is assigned, $\lambda(\theta_1)$ must be shown to have exactly
four monotone phases — two maxima and two distinct minima — which is what makes
the motion Atkinson rather than a stair-stepped artefact. The cycle is then
adiabatic compression and expansion about a constant-volume heat release,

$$Q = \frac{V_0 (P_3 - P_2)}{\gamma - 1}$$

with $Q$ fixing the fuel consumed. Deriving the fuel from the heat release keeps
the range model consistent with the thermodynamics rather than bolting a second
combustion model alongside it.

### 3.2.2 Mechanical losses

§2.3 noted that $\eta$ loses nothing. The real losses come from quantities the
equilibrium solve already produces — joint reactions $R_j$ turning through
relative angles $\Delta\theta_j$, and the liner reaction $D$ that the side-load
constraint already bounds:

$$W_{\text{bearings}} = \sum_j \mu_j r_j \oint |R_j| \left|\frac{d(\Delta\theta_j)}{d\theta_1}\right| d\theta_1$$

$$W_{\text{piston}} = \mu_p \oint \big(|D| + F_{\text{ring}}\big)
   \left|\frac{d\lambda}{d\theta_1}\right| d\theta_1$$

This is what makes the constraint set mean something: without it, the side-load
and bearing-load bounds are assertions about wear that never enter any
objective. With it, a design that leans on the liner burns its fuel on the
liner.

### 3.2.3 Engine mass

The sized members are 0.15 kg; no engine weighs 0.15 kg. Optimising that number
optimises a tail. The budget carries eight items, of which two change the shape
of the problem rather than its scale.

**The crankcase converts the envelope into kilograms.** A box encloses the
mechanism and its walls scale with $H \times B$, so the two envelope objectives
of §2.3 become mass at a rate the physics fixes rather than the designer.

**The flywheel converts torque ripple into kilograms.** A single cylinder needs
rotating inertia to carry it through compression,

$$J = \frac{\Delta E}{\delta\, \omega^2}, \qquad
\Delta E = \max_\theta E(\theta) - \min_\theta E(\theta), \qquad
E(\theta) = \int_0^\theta (M_r - \bar M_r)\, d\theta_1$$

so a linkage whose torque curve is flatter is lighter — a driver no geometric
constraint expresses, pushing against a long lever arm. The turning-moment
diagram here must be the **gas** torque: the inertia part is energy traded with
the mechanism's own masses, already accounted for by the inherent rotating
inertia, and including it overstates the flywheel fivefold at low speed.

### 3.2.4 From engine to distance

These cars are driven **burn and coast**: run hard from $v_{lo}$ to $v_{hi}$,
declutch, coast back. With $M v\, dv/dx = F$,

$$d = \int_{v_{lo}}^{v_{hi}} \frac{M v\, dv}{F(v)}, \qquad
  t = \int_{v_{lo}}^{v_{hi}} \frac{M\, dv}{F(v)}$$

with $F = P_w/v - F_{\text{res}}$ under power and $F = F_{\text{res}}$ coasting.

The bookkeeping is worth stating because the naive expectation is wrong. Over
one closed cycle the car starts and ends at the same speed, so the kinetic
energy nets to zero and the propulsive work equals the resistance work over the
**whole** distance: hard acceleration costs nothing in road load. What the
strategy buys is that the engine runs at high load; what it costs is
aerodynamic.

The minimum-average-speed rule is active at every optimum where the engine has
power to spare, which collapses the two-dimensional strategy search to one
dimension — and matters beyond speed, since a grid search would make the
objective a step function of the design variables and the optimizer downstream
would be differentiating quantisation noise.

## 3.3 Limitation: a quasi-static model has no parts

**The limitation.** Nothing in $X \to \lambda \to p \to \text{loads}$ determines
a cross-section. Without sections there is no mass, without mass no inertia, and
the objective above cannot be evaluated at all.

**The resolution.** Size the members against yield, fatigue and buckling. But
the sizing needs the loads and the loads need the masses, so with
$d$ the section diameters and $N, M$ the internal loads,

$$d = \mathcal{S}\big(N(d),\, M(d)\big)$$

is a fixed point rather than a sequence. It has one, and plain iteration reaches
it: a bending-critical member needs $d \sim F^{1/3}$, so $m \sim F^{2/3}$, and
composing with $F \sim m a$ gives $m \sim (Ca)\, m^{2/3}$ — a sub-linear loop
gain, hence a fixed point at $m = (Ca)^3$. The cubic dependence on acceleration
is why engine speed dominates the answer.

**Consequence.** The problem is multidisciplinary in the strict sense, and §3.5
must choose an architecture for it. Derivations: {doc}`theory` §9.

## 3.4 Limitation: the feasible set has measure zero

**The limitation.** Two requirements are stated as equalities,

$$STE(X) = 74\ \mathrm{mm}, \qquad \varepsilon(X) = 16$$

so the feasible set $\{X : g_{\text{eq}}(X) = 0\}$ is a nine-dimensional
manifold in $\mathbb{R}^{11}$. Its Lebesgue measure is zero, and a point drawn
from any continuous distribution lies on it with probability zero.

Measured consequences:

| | |
|---|---|
| COBYLA, 120 evaluations, 313 s | returns its starting point unchanged |
| uniform samples over the global box | 0 feasible in 4000 |
| uniform samples within 50 % of a feasible design | 0 feasible in 4000 |
| uniform samples within 10 % of a feasible design | 0 feasible in 4000 |

**Resolution (a): the optimizer.** Only a method that *moves along* the manifold
can be used, which selects SQP with exact gradients (§2.6). This is what makes
§3.5 necessary rather than merely desirable.

**Resolution (b): the specification.** A requirement stated as an equality
cannot be met by a manufactured part: no dimension is produced exactly, so the
stroke of a built engine is a random variable. The equality is a shorthand, and
restoring what it means gives

$$|STE - 74| \le \delta_{STE} = 0.05\ \mathrm{mm}, \qquad
  |\varepsilon - 16| \le \delta_\varepsilon = 0.05$$

The feasible set becomes full-dimensional and the problem well posed.

**Consequence.** The relaxation is a statement about tolerance, and §3.8 shows
it cannot be left unexamined.

## 3.5 Limitation: finite differences are wrong, not merely inaccurate

**The limitation.** Several constraints are extrema over the crank revolution,

$$\gamma(X) = \max_{\theta_1} \frac{|D(X,\theta_1)|}{\max_{\theta_1} |P(\theta_1)|}$$

and the sample attaining the maximum changes as $X$ moves. A difference quotient
taken across that switch does not approximate a derivative: measured at 25 %
error on $\gamma$ at a $10^{-4}$ mm step, and it does not improve with step size.

**The resolution.** The chain is closed form, so derivatives propagate forward
alongside each intermediate. For an extremum attained at $\theta^*$ the
envelope theorem gives

$$\frac{d}{dX}\Big[\max_{\theta} f(X,\theta)\Big]
  = \frac{\partial f}{\partial X}\Big|_{\theta^*}$$

because the term through $d\theta^*/dX$ carries $\partial f/\partial\theta = 0$
at the maximiser. The switching problem disappears rather than being mitigated.

Through the fixed point of §3.2 the same idea applies twice more. The
$18\times18$ equilibrium solve differentiates as

$$\frac{\partial x}{\partial p} = A^{-1}\Big(\frac{\partial b}{\partial p}
  - \frac{\partial A}{\partial p}\,x\Big)$$

reusing the factorisation already computed; and the sizing bisection is never
differentiated at all, because the diameter is defined implicitly by
$U(d, N, M) = 1$ and the implicit function theorem gives

$$\frac{\partial d}{\partial q} = -\frac{\partial U/\partial q}{\partial U/\partial d}$$

**Consequence.** SLSQP becomes applicable to a problem on which derivative-free
search returns its input. Derivations: {doc}`theory` §10.

## 3.6 Limitation: the coupling is a field, not a handful of scalars

**The limitation.** §3.3 leaves an MDA to be placed. §2.4 says the MDF/IDF trade
is decided by the dimension of the coupling, so it must be counted:

$$\dim(\text{coupling}) = \underbrace{2 \times 7 \times 360 \times 9}_{\text{load histories}}
   + \underbrace{7}_{\text{diameters}} + \underbrace{1}_{\text{piston mass}}
   = 45\,368$$

against **11** design variables. The couplings are not scalars but the internal
load history of every member at every crank angle at every station.

**The resolution.** MDF. IDF would carry 45 367 extra design variables and as
many consistency constraints to optimise eleven degrees of freedom; it is not
slower, it is unavailable.

**How coupled is it?** Gauss–Seidel converges linearly at a rate that *is* the
coupling strength, $\rho = \lim \|r_{k+1}\|/\|r_k\|$, so the claim is
measurable rather than asserted. At rest $\rho = 0$ exactly — with no inertia
there is no path from mass to load and the quasi-static problem is recovered —
rising to 0.68 at 1500 rpm. §6 tabulates it.

## 3.7 Limitation: the gear choice is discrete and pins a design variable

**The limitation.** A gear has an integer tooth count cut with a standard-module
hob. For the 2:1 pair, $r = mz/2$ and $z_1 = 2z_2$, so

$$I = \tfrac32\, m\, z_2, \qquad m \in \text{ISO 54},\quad z_2 \in \mathbb{Z},\ z_2 \ge 17$$

$I$ therefore lives on a lattice, not an interval — and $I$ is one of the
variables the equalities of §3.3 are satisfied *with*, so choosing the gears
throws the design off both. A 0.18 mm snap to the nearest lattice point moves
the top-dead-centre gap from 0.003 mm to 0.058 mm.

**The resolution.** State it as the MINLP it is and decompose. Of the two
candidates in §2.7, outer approximation is chosen because of what each master
requires:

- a **Benders** cut needs the optimal-value sensitivity $d\theta/dI$, which
  requires multipliers this problem's degenerate active set does not supply;
- an **outer-approximation** cut needs only $\nabla f$ and $\nabla g$ at the
  visited point, which §3.5 already provides exactly.

The master is then, with $y$ a one-hot selection over the lattice,

$$\min_{x,y,\eta}\ \eta \quad\text{s.t.}\quad
\begin{cases}
\eta \ge f_k + \nabla f_k^{\mathsf T}(x - x_k) \\
0 \ge g_k + \nabla g_k^{\mathsf T}(x - x_k) \\
I = \sum_j I_j y_j,\quad \textstyle\sum_j y_j = 1
\end{cases}$$

Infeasible sub-problems need no special machinery: their constraint
linearisations are added without an objective cut, which excludes that lattice
point on evidence.

## 3.8 Limitation: the relaxed bands are comparable with the scatter

**The limitation.** §3.4 relaxed the equalities into bands of half-width
$\delta$. That was a promise about tolerance, and it can be checked. With
$\Sigma$ the covariance of the dimensional errors and $\nabla g$ from §3.4,

$$\sigma_g = \sqrt{\nabla g^{\mathsf T} \Sigma\, \nabla g}$$

Evaluated at IT8 machining tolerances:

| | |
|---|---|
| band half-width $\delta_{STE}$ | 0.050 mm |
| scatter $\sigma_{STE}$ | 0.029 mm |
| best achievable index $\delta/\sigma$ | 1.7 |
| index of the reference design | 0.68 |

The band is under two standard deviations wide, so a perfectly centred design
misses it about 9 % of the time. **The relaxation that made the problem solvable
is the same quantity that makes it unreliable**, and a deterministic optimizer
satisfying $|STE-74| \le 0.05$ has no way to notice.

**Why a fixed margin is the wrong repair.** $g + k\sigma_g \le 0$ per constraint
is a reliability statement only under independence. Here every constraint is a
function of the same eleven dimensions and the measured correlations reach
$0.94$, with exactly $-1$ between the two sides of a relaxed band.

**The resolution.** Constrain a probability of failure, keeping the correlation:

$$P_f = 1 - \Phi_n(\beta;\, \rho) \le p_{\text{target}}, \qquad
\beta_i = -\frac{g_i}{\sigma_i}, \qquad
\rho_{ij} = \frac{\nabla g_i^{\mathsf T} \Sigma \nabla g_j}{\sigma_i \sigma_j}$$

a first-order (FORM) index per constraint combined through the multivariate
normal orthant. It costs one Jacobian evaluation, so it can sit inside the
optimization; sampling is the reference against which it is checked.

## 3.9 Limitation: the problem is nonconvex and the answer is one local solution

**The limitation.** Everything above yields a local solution. §2.9 lists the
global strategies, all of which require generating feasible starting points —
which §3.4 shows uniform sampling cannot do.

**The resolution.** Construct restarts *on* the manifold: perturb the incumbent,
project the perturbation back onto the two equalities by a minimum-norm Newton
step from the analytic Jacobians,

$$\Delta X = -J^{+} r, \qquad J = \begin{bmatrix}\nabla STE\\ \nabla \varepsilon\end{bmatrix}$$

and let the optimizer restore the inequalities from there.

**Consequence, stated in advance.** This is a local-search diversification, not
a global method. §6 reports what it settles and §7 what it does not.

## 3.10 The problem solved

$$
\begin{aligned}
\max_{X,\,y} \quad & R(X, y) && \text{range [km/L]} \\
\text{s.t.}\quad
& g_i(X) \le 0, \quad i \in \{mra,\, W,\, g,\, d,\, \gamma\}
  && \text{geometric} \\
& |STE(X) - 74| \le \delta_{STE}, \quad |\varepsilon(X) - 16| \le \delta_\varepsilon
  && \text{relaxed equalities (§3.4)} \\
& s(X, y) \le 0,\; \ell(X, y) \le 0,\; b(X, y) \le 0
  && \text{saturation, slenderness, bearing} \\
& r(X, y) \ge 0, \quad h(X, y) \ge 0
  && \text{engine runs, gears fit} \\
& I = \tfrac32 m z, \quad m \in \text{ISO 54},\ z \in \mathbb{Z},\ z \ge 17
  && \text{catalogue (§3.7)} \\
& X \in [X_{lb}, X_{ub}] \subset \mathbb{R}^{11}
\end{aligned}
$$

Twelve constraints, and $y$ is the converged MDA state of §3.6. Two details
the code makes visible and the notation above does not: each two-sided band is
attached as a pair of one-sided inequalities, so the scenario carries fourteen
constraint *functions* for these twelve constraints; and $I$ is not a free
variable but an output of the catalogue relation, leaving ten in the search.
The envelope bounds $H$ and $B$ that §2.3 discusses are *not* among the
constraints at all: once the objective prices size through mass, a separate
limit on it is redundant, and the range problem attaches none.

### Imposed, or merely checked

The statement above is the problem. Whether an optimizer is made to *hold* it
is a separate question, and the study answers it both ways:

| | the coupled and vehicle rows | reliability | results |
|---|---|---|---|
| {py:func}`~exlink.scenarios.build_range_scenario` | bind only at the end | audited after | §6.1 – §6.3 |
| {py:func}`~exlink.synthesis.maximise_range_from_target` | imposed at every step | audited after | §6.4 |
| the same, with ``beta_target`` | imposed at every step | **constrained** | **§6.4** |

The difference is not academic: the same SLSQP on the same problem reaches
4.9 % further under the second, and 1.7 % under the third, which also holds a
reliability target the other two only measure (§6.4). A constraint that binds
only at the end
lets the search spend its whole trajectory in a region it will later be told it
cannot use.

The second form needs one thing the first does not, and it is the subject of the
rest of this section.

### Why the objective needs a fallback

$R$ is not computable everywhere. A design whose kinematics closes can still
fail to size, fail to run, or fail to produce a four-stroke motion at all, and
at such a point the objective has no value for a line search to descend. A
constant penalty there leaves a flat region with no gradient — which is exactly
the failure the reliability constraint hits below, and the reason that search
could not move.

So the objective is a ladder, and a prescribed motion $\lambda^{\star}$ holds its
middle rung:

| the design | scores |
|---|---|
| range computable | $-R(X)$, the real objective |
| analysable, no range | a floor plus $\lVert\lambda(X) - \lambda^{\star}\rVert^2$ — the target takes over |
| motion is not a four-stroke cycle | a larger constant |

Each rung is strictly worse than the one above, so the search is pushed back
towards designs that run; and on the middle rung it still has something to
follow, because tracking $\lambda^{\star}$ is a proxy for returning to a working
cycle.

The target is a *fallback*, not a constraint. It is constructed to satisfy both
equality requirements exactly — they are functionals of $\lambda$ alone, so a
motion with the right two strokes meets them before any linkage exists — and it
is abandoned the moment the range becomes computable. §6.4 measures both halves
of that: the residual it ends at, and how often it fires.

§7.4 collects what the prescribed motion taught along the way, including the
attainability limit that decides whether a target is usable at all.

### Reliability: audited, not imposed

Under either form above this is a **deterministic** problem, and every design
reported in §6 came from one of them. The reliability analysis of §3.8 is
applied *to the solution*, so §6.2's 0.645 is an audit of a design obtained
without it rather than a target that was met or missed.

The loop can be closed:
{py:func}`~exlink.scenarios.build_range_scenario` takes a ``beta_target``, which
attaches {py:class}`~exlink.robustness.FailureProbabilityDiscipline` and adds

$$\beta_{\text{sys}}(X) \;=\; -\Phi^{-1}\bigl(P_f(X)\bigr) \;\ge\; \beta_{\text{target}}$$

as a thirteenth constraint. The index and not the probability: $P_f$ saturates
towards 0 and 1, where its gradient vanishes and the optimizer stalls, while
$\beta$ stays well scaled throughout. The discipline depends on the design
variables only, so it does not join the MDA and §3.6's coupling count is
unchanged.

**Running it does not produce a more reliable design, and the reason is
instructive.** From the coupled reference design at $\beta = -0.373$, SLSQP
returns its starting point unchanged — not only for a demanding target but for
$\beta \ge -0.2$, a step of 0.17 — reporting a positive directional derivative
for the line search. The reliability gradient is not at fault: central
differences on $\beta$ are smooth in all eleven variables at this point, with no
discontinuity. What defeats the line search is the same thinness §3.4 measures
elsewhere. A step of 0.05 mm along the normalised $\nabla\beta$ leaves the
geometric constraint set entirely, at which point $\beta$ is reported as $-8.2$
by the penalty fallback:

| step along $\nabla\beta$ [mm] | $\beta$ | nominally feasible |
|---|---|---|
| 0.00 | $-0.373$ | yes |
| 0.05 | $-2.570$ | no |
| 0.10 and beyond | $-8.210$ | no |

That better reliability is nevertheless available nearby is shown by sampling
rather than by descent, and the margin is not small. Drawing 2500 designs at
0.05 % to 1 % relative scatter about the reference, scoring each by $\beta$, and
then checking the twenty-five best against the *full* constraint set:

| | $\beta$ | $P_f$ | range |
|---|---|---|---|
| the reference design | $-0.373$ | 0.645 | 3338.3 km/L |
| best sampled, fully feasible | $\mathbf{+0.502}$ | **0.308** | **3341.7 km/L** |

All twenty-five are feasible. The probability of missing a requirement is more
than halved, and it costs **nothing**: the best of them is 0.10 % *better* in
range than the design it replaces.

**The deterministic optimum is therefore dominated, not merely unreliable.**
There is no reliability-against-range trade-off to negotiate here — a random
cloud finds points that beat the converged design on both objectives at once.
The reason is the one RBDO exists for: nothing in the deterministic problem
rewards backing away from a constraint boundary, so SLSQP converges *onto*
$g = 0$, where half the scatter falls outside. Retreating a few hundredths of a
millimetre costs almost no range and buys most of the reliability back, and the
deterministic formulation cannot see that because it prices only the nominal
point.

So the failure above is entirely algorithmic. It is not that the reliable
region is expensive to reach; it is adjacent, feasible, and free.

So closing the reliability loop is not a matter of attaching the discipline,
which is easy, but of reaching the reliable region, which a gradient method
starting from the deterministic optimum does not do. §7.3 lists what that needs.
The deterministic problem above therefore remains the one this study solves, and
the reliability numbers of §6.2 remain an audit — but now by a measured
limitation rather than an unexamined omission.

### Which constraints the probability covers

The uncertainty model is $\Sigma = \operatorname{diag}(\sigma^2)$ over the
eleven dimensions, from ISO 286 grades plus an angular clocking term
({py:func}`exlink.robustness.covariance`). It contains no material, load or
friction scatter. That fixes which constraints can honestly carry a probability:

| constraint | in $P_f$ | why |
|---|---|---|
| $STE$, $\varepsilon$ bands | **yes** | closed form in $X$; $\Sigma$ is their *complete* uncertainty |
| $mra$, $W$, $g$, $d$, $\gamma$ | **yes** | same |
| bearing | no | variance dominated by material and load scatter, which $\Sigma$ omits |
| saturation | no | a ceiling indicator: non-smooth, so the FORM linearisation carries no information at the kink |
| slenderness | no | a *model-validity* guard, not a limit state — see below |
| engine runs, gears fit | no | load-dependent, same omission as bearing |

Seven of the twelve, and the seven are exactly those that are functions of the
eleven dimensions alone. For them $\Sigma$ is the whole story and $\beta = -g/\sigma$
means what it says.

The exclusions are not one kind of thing. **Slenderness is a category error, not
a gap**: it fires when a link grows thicker than a third of its length, at which
point sizing it as a beam has stopped being credible. There is no probability
that beam theory applies; a limit state is a statement about a part failing, and
this is a statement about the model's domain. **Saturation is a genuine
limitation** — FORM linearises $g$ at the design point, and at a ceiling that
linearisation is meaningless; it needs sampling or a smoothed reformulation.
**Bearing is the weakest exclusion**, and the gradient it would need is already
available exactly
({py:attr}`exlink.dynamics_jacobian.CoupledJacobian.peak_bearing_load`); what is
missing is not the derivative but the uncertainty model behind it.

That last point generalises, and it is the reason the answer here is seven
rather than twelve. A probability of failure computed from an uncertainty model
that omits the dominant source is *worse* than a deterministic margin, because
it launders a partial variance into something that reads as a reliability
statement — and it would then enter the system union of §3.8 and corrupt a
number that is currently defensible. `bearing $\le 0$` with a safety factor is
at least honest about being a margin. Widening $\Sigma$ to carry material and
load scatter is the prerequisite for widening $P_f$, and §7.3 lists it as such.


---

Next: [4. Implementation framework](framework.md)
