# 2. State of the art

## 2.1 Scope of the review

The problem stated in §1 has five features that a solution method must
accommodate: it is a *mechanism-synthesis* problem, so the literature on how
such problems are posed governs the choice of objective; the objective is an
*application-level* quantity, so the vehicle- and engine-modelling literature
supplies the chain that prices geometry; the resulting model is *coupled*, so
the MDO-architecture literature governs how the disciplines are assembled; the
gear pair is *discrete*, so the MINLP literature governs the outer loop; and the
requirements are *tolerance-limited*, so the reliability literature governs how
they are enforced.

This section reviews each in turn. Every subsection states what the literature
established, what it left open for a problem of this shape, and which option is
therefore admissible here; the choices themselves are made in §3. §2.10 states
the position of this work relative to the whole.

## 2.2 Extended-expansion engines and the mechanisms that realise them

That an expansion ratio larger than the compression ratio raises the ideal cycle
efficiency has been known since Atkinson's patents, and the modern literature is
concerned with how to realise the asymmetry mechanically. {cite:t}`zhao2017`
surveys the alternatives and separates the two families: *valve-timed*
realisations (the Miller cycle, late or early inlet-valve closing), which shorten
the effective compression stroke and lose displacement, and *geometric*
realisations, in which a linkage gives the piston physically unequal strokes on
alternate revolutions.

The mechanism studied here belongs to the second family. Its canonical
description is {cite:t}`watanabe2006`, who analyse the multiple-linkage
general-purpose engine that became Honda's EXlink: an eccentric on a
half-speed shaft drives a swing link that modulates the connecting-rod
attachment, yielding a compression ratio of 8.5 against an expansion ratio of
12.3 and an indicated thermal efficiency raised from 27.3 % to 31.3 %.
{cite:t}`dumboeck2018` give an independent treatment of the same concept. The physical model used in §3.2 follows the standard
treatment of engine friction, heat release and mean effective pressures in
{cite:t}`heywood1988`.

Two things are settled by this literature and one is not. Settled: the
thermodynamic benefit is real, and it is measured at a *fixed* mechanism —
Watanabe et al. analyse a linkage whose dimensions are given, not chosen.
Also settled: the half-speed shaft that produces the asymmetry halves the
firing frequency. Not settled: what the linkage costs. None of these studies
carries the linkage dimensions through to the mass of the parts that must
carry the inertia loads at speed, so none can say whether the thermodynamic
gain survives the mechanism that delivers it. §6.3 carries it through, and
finds the advantage robust to how the two engines are matched but its
*attribution* fragile. At matched power strokes per minute the linkage keeps a
13–19 % range advantage over an optimised conventional engine; scoring the same
linkage as if it drove a two-revolution gas exchange turns that into a 4.3 %
deficit, which locates the benefit in the cycle rate rather than in the extended
expansion; and holding the conventional engine to the *same specification* as
the linkage widens the figure further. Which question is being asked has to be
stated with the number.

## 2.3 Formulating a mechanism-synthesis problem

The dominant formulation in mechanism synthesis is *dimensional synthesis to a
kinematic criterion*: choose link lengths to minimise a path or motion error, or
to maximise a transmission quality index, subject to bounds on the envelope.
{cite:t}`balli2002` review the transmission angle as such an index and its use
as both an objective and a constraint; {cite:t}`gosselin1990` give the
Jacobian-based characterisation of the configurations where force transmission
degenerates, which is what a transmission-angle criterion is a scalar proxy for.
{cite:t}`cabrera2002` are representative of how these problems are solved:
a genetic algorithm on a bounded box, the criterion evaluated kinematically.

The conventional treatment of *this* mechanism follows that pattern, maximising
a lever-arm quality measure

$$\eta = \frac{\int_0^{2\pi} M_r \, d\theta_1}{2 (STE + STC)\, \bar P}$$

subject to bounds on the two envelope dimensions $H$ and $B$, over the eleven
linkage variables. Three limitations follow, and they are limitations of the
formulation rather than of its solution.

It is **multi-objective without an exchange rate**. $\eta$, $H$ and $B$ compete
and nothing prices one against another, so the formulation yields a Pareto front
and never a design. Weighted sums, $\varepsilon$-constraint and moving limits
each produce *a* point, but the weights are the designer's, not the physics'.

Its central quantity is **not an efficiency**. With no friction in the model the
virtual-work identity makes $\int M_r \, d\theta_1 \equiv \int P \, d\lambda$ at
every crank angle, so $\eta$ is a ratio of two provably equal works: a kinematic
quality measure in which nothing is lost. This is a general property of
quasi-static transmission indices, not a defect peculiar to this one — the
literature treats them as *proxies* for goodness, and they are only ever as good
as the correlation between the proxy and the quantity of interest.

It **cannot see the parts**. Nothing in a quasi-static formulation determines a
cross-section, so the mechanism has no mass, no inertia loads, and no
speed dependence.

The alternative adopted here — carry the analysis through to the quantity the
application scores — is standard practice in vehicle-level MDO rather than a
novelty, and the ultra-efficiency-vehicle literature supplies the missing link
in the chain: {cite:t}`gechev2020` model an Eco-marathon prototype over its
actual track and show that the achievable consumption is set jointly by the
powertrain map and the driving strategy, so that neither can be scored without
the other. What that literature optimises is the *strategy* at fixed hardware;
what is optimised here is the *hardware*, with the strategy solved to optimality
inside each evaluation.

## 2.4 Architectures for a coupled problem

Once mass, inertia and structural sizing enter, the model is coupled: the
sizing depends on the loads and the loads depend on the sizing.
{cite:t}`cramer1994` established the formal vocabulary — multidisciplinary
feasible (MDF), individual discipline feasible (IDF), all-at-once — and
{cite:t}`martins2013architectures` give the modern survey, including the
bi-level family (of which {cite:t}`bliss2000` is the canonical member) that
decomposes system-level and discipline-level decisions.

**MDF** converges an MDA at every optimizer iteration, so every evaluated point
is physically consistent and the design space stays small; the cost is the inner
iteration and the need to differentiate through it.

**IDF** promotes the coupling variables to design variables with consistency
constraints, so no inner iteration runs; the cost is a design space that grows
by the dimension of the coupling.

The survey's trade favours IDF when the MDA is expensive relative to the
optimizer's handling of a larger space. That trade is decided here by a single
number — the dimension of the coupling — and §3.6 shows it is not close.

## 2.5 Obtaining derivatives

Whether a gradient method is admissible at all depends on whether accurate
derivatives are available. {cite:t}`martins2013derivatives` unify the available
routes; {cite:t}`sobieski1990` gives the global sensitivity equations that
extend them to coupled systems, which is what any gradient of an MDF
formulation must ultimately solve.

**Finite differences** are the default and are *wrong* on this problem, not
merely inaccurate: several constraints are maxima over the crank revolution, and
the sample attaining the maximum switches as the design moves. A difference
quotient taken across the switch is meaningless — measured at 25 % error on the
side-load ratio at a $10^{-4}$ mm step.

**Complex-step** {cite:p}`martins2003complexstep` removes subtractive
cancellation but not the switching problem, and requires the whole chain to be
complex-safe.

**Algorithmic differentiation** {cite:p}`griewank2008` is applicable and would
avoid hand derivation, at the cost of a dependency and of differentiating
through the fixed-point solver unless it is taught not to.

**Analytic derivatives with the envelope theorem** are used here. The chain is
closed form, so forward-mode propagation is direct; and for a maximum attained
at $\theta^{\star}$ the derivative is the partial derivative evaluated there,
because the term through the moving maximiser carries
$\partial f/\partial\theta = 0$. The result is classical —
{cite:t}`danskin1966` for max-functions and {cite:t}`milgrom2002` for the
general envelope theorem — and it removes the switching problem rather than
mitigating it. §3.5 develops it.

## 2.6 Solving a problem whose feasible set is thin

For a feasible set defined partly by equalities, the available families are:

| family | examples | applicability here |
|---|---|---|
| derivative-free direct search | COBYLA {cite:p}`powell1994`, Nelder–Mead | **no** — proceeds by sampling; see §3.4 |
| population methods | differential evolution, NSGA-II {cite:p}`deb2002` | **no** — same reason |
| penalty / augmented Lagrangian | external penalty, AL {cite:p}`nocedal2006` | possible, but inherits the conditioning |
| SQP with exact gradients | SLSQP {cite:p}`kraft1988` | **yes** — moves along the manifold |
| interior point | IPOPT | yes in principle; not evaluated here |

The distinction is not efficiency. A sampling method evaluates points that lie
off a measure-zero set with probability one, so it has no feasible point to
improve from; §3.4 gives the measurement. This is precisely why the
mechanism-synthesis literature reviewed in §2.3 can use genetic algorithms and
this study cannot: {cite:t}`cabrera2002` optimise over a box in which every
sampled point is feasible, whereas the two equality requirements here reduce
the feasible set to a nine-dimensional manifold in $\mathbb{R}^{11}$.

## 2.7 Mixed-integer nonlinear programming

The gear pair is discrete, so the problem is an MINLP.
{cite:t}`belotti2013` survey the field; the families relevant to a small
catalogue and an expensive continuous sub-problem are:

| family | note |
|---|---|
| exhaustive enumeration | exact on a small catalogue; no bound, no stopping criterion |
| relax-and-round | cheap; the rounded point generally leaves the feasible set |
| branch and bound | needs a relaxation whose bound is meaningful |
| generalized Benders (GBD) {cite:p}`benders1962,geoffrion1972` | master built from optimal-value sensitivities (Lagrange multipliers) |
| outer approximation (OA) {cite:p}`duran1986,fletcher1994` | master built from linearisations of $f$ and $g$ at the sub-problem solution |

GBD and OA share a decomposition — a continuous NLP for fixed integers, a
MILP master that accumulates cuts — and differ in what the master is built from.
{cite:t}`fletcher1994` also supply the treatment of *infeasible* sub-problems
that a thin feasible set makes routine here: a feasibility cut derived from the
sub-problem's own linearisation, rather than a heuristic. OA is chosen in §3.7,
for reasons that turn on which quantities this problem can supply reliably.

A separate strand treats discrete *catalogue* choices by continuous relaxation
with a penalty that drives the relaxation back to a vertex — DMO
{cite:p}`stegmann2005` and SFP {cite:p}`bruyneel2011` in composite design, both
descended from the SIMP penalty of {cite:t}`bendsoe1989`. That machinery is the
right tool when the catalogue is large and the sub-problem cheap; it is the
wrong one here, where the catalogue has a handful of entries and each
sub-problem solve is the expensive object, because it converts an exactly
solvable outer problem into an approximately solvable one.

## 2.8 Design under tolerance

The requirements in §5.2 are stated as equalities that no manufactured part can
meet exactly, so they must be relaxed to bands and the bands must then be
justified. Three levels of treatment are available, in increasing fidelity and
cost.

**Safety factors** applied to the nominal constraint. Simple, but the implied
reliability is unknown and varies constraint to constraint.

**Worst-case or fixed-margin robust design**, $g + k\sigma_g \le 0$.
{cite:t}`beyer2007` survey this family. It computes a margin, not a
probability, and it is a reliability statement only if the constraints are
independent — §3.8 shows they are not, with correlations reaching $\pm 1$
because the constraints are analytic functions of the same eleven dimensions.

**Reliability-based design optimization (RBDO)**, constraining a probability of
failure. The first-order reliability method of {cite:t}`hasofer1974`, made
algorithmic by {cite:t}`rackwitz1978`, gives the per-constraint reliability
index $\beta = -g/\sigma_g$ at a cost of one gradient. Its use *inside* an
optimization is the subject of {cite:t}`tu1999`, who separate the
reliability-index and performance-measure approaches, of {cite:t}`du2004`, who
decouple the nested loops, and of the survey by {cite:t}`valdebenito2010`.
Because the failure event here is a *system* event — any one of several
requirements violated — the union probability is needed rather than the
per-constraint one; {cite:t}`ditlevsen1979` gives the classical bounds, and
the correlated multivariate-normal orthant that the bounds approximate is
computed directly by the transformation of {cite:t}`genz1992`. Tolerance
classes are taken from {cite:p}`iso286`.

§3.8 uses FORM with the full correlation matrix, checked against sampling. It
is applied as an audit of the solved design rather than as a constraint of the
search, and over the seven constraints whose uncertainty the model actually
carries rather than over all twelve; §3.10 states both restrictions and why.

## 2.9 Escaping local optima

Multistart {cite:p}`rinnooykan1987`, basin hopping, and global methods (DIRECT,
simulated annealing, population methods) all require the ability to *generate a
feasible starting point*. The multi-level single-linkage theory of
{cite:t}`rinnooykan1987` in particular assumes sampling from the feasible
region. §3.4 shows that uniform generation cannot do so here, which restricts
the available options to restarts constructed *on* the feasible manifold.

## 2.10 Position of this work

Against that literature, this study is conventional in its machinery and
unconventional only in what it applies the machinery to.

*Standard practice, used as such.* The architecture comparison follows
{cite:t}`martins2013architectures`; the analytic sensitivities follow
{cite:t}`martins2013derivatives`; the outer loop is the algorithm of
{cite:t}`duran1986` and {cite:t}`fletcher1994` as implemented in GEMSEO
{cite:p}`gallard2018`; the reliability treatment is FORM
{cite:p}`hasofer1974,rackwitz1978` with a system probability
{cite:p}`ditlevsen1979,genz1992`. Nothing in §3 is a new algorithm.

*What the literature leaves open.* The extended-expansion literature
(§2.2) sizes no parts and therefore cannot price the mechanism against the cycle
it enables. The mechanism-synthesis literature (§2.3) optimises transmission
quality, which the virtual-work identity shows to be a proxy with no loss in it,
and does so on a box in which sampling methods are admissible. The
Eco-marathon literature (§2.3) optimises the strategy at fixed hardware. The
RBDO literature (§2.8) treats the bounds as given data rather than asking which
of them the mechanism is able to hold, and at what price in the objective.

*What this study contributes.* Three things, each a consequence of joining
those strands rather than of extending any one of them.

1. **A priced formulation.** Carrying the mechanism through mass, friction,
   structural sizing and the vehicle model to distance per unit fuel supplies
   the exchange rate that §2.3 shows the conventional formulation lacks, and
   makes the extended-expansion benefit separable from the firing-frequency
   effect that accompanies it in every geometric realisation. Priced against a
   conventional engine optimised under the same models rather than proportioned
   by hand, the topology's advantage survives matching the firing rate while the
   extended-expansion component alone reverses sign (§6.3).

2. **A demonstration that the geometry of the feasible set selects the
   algorithm.** The measure-zero feasible set is not a nuisance to be
   penalised: it excludes the sampling methods that the synthesis literature
   uses by default, forces multistart onto the manifold, and — through the
   band relaxation it obliges — makes reliability analysis a consequence of the
   formulation rather than an addition to it (§3.4, §3.8).

3. **A result about the requirements rather than about a design.** Asking not
   "what is the failure probability at this tolerance?" but "which of these
   bounds governs the reliability, and what does widening it cost in the
   objective?" inverts the usual RBDO question. Answering it here identifies
   one bound finer than the model's own resolution and one that governs
   everything else, and prices both in range (§6.2). Neither is visible in any
   nominal quantity, and optimization at a fixed specification cannot reveal
   them.

The problem is deliberately small enough that every claim above is checked
against a closed-form or sampled reference (§4.4), which is what makes it usable
as a demonstration rather than only as a design.

---

Next: [3. Methodology](methodology.md)
