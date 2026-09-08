# 2. Related work

Five features of the problem stated in §1 govern which methods may be applied to
it: it is a mechanism-synthesis problem, so the literature on how such problems
are posed governs the choice of objective; the objective is an application-level
quantity, so the vehicle- and engine-modelling literature supplies the chain that
prices geometry; the model is coupled, so the MDO-architecture literature governs
how the disciplines are assembled; the gear pair is discrete, so the MINLP
literature governs the outer loop; and the requirements are tolerance-limited, so
the reliability literature governs how they are enforced. Each is reviewed below
for what it settles and what it leaves open at this problem's shape.

## 2.1 Extended-expansion engines and the mechanisms that realise them

That an expansion ratio larger than the compression ratio raises the ideal cycle
efficiency has been known since Atkinson's patents, and the modern literature is
concerned with realising the asymmetry mechanically. {cite:t}`zhao2017` separates
the two families: *valve-timed* realisations (the Miller cycle, late or early
inlet-valve closing), which shorten the effective compression stroke and lose
displacement, and *geometric* realisations, in which a linkage gives the piston
physically unequal strokes on alternate revolutions. The mechanism studied here
belongs to the second family, whose canonical description is
{cite:t}`watanabe2006` — the multiple-linkage general-purpose engine that became
Honda's EXlink, with a compression ratio of 8.5 against an expansion ratio of
12.3 and an indicated thermal efficiency raised from 27.3 % to 31.3 %.
{cite:t}`dumboeck2018` give an independent treatment of the same concept, and the
physical model of §3.3 follows the standard treatment of engine friction, heat
release and mean effective pressures in {cite:t}`heywood1988`.

Two things are settled by this literature and one is not. The thermodynamic
benefit is real, and it is measured at a *fixed* mechanism: Watanabe et al.
analyse a linkage whose dimensions are given rather than chosen. The half-speed
shaft that produces the asymmetry does not change the firing frequency, power
being taken from the crankshaft as usual. What is not settled is what the linkage
*costs*: none of these studies carries the linkage dimensions through to the mass
of the parts that must carry the inertia loads at speed, so none can say whether
the thermodynamic gain survives the mechanism that delivers it. §5.4 carries it
through.

## 2.2 Formulating a mechanism-synthesis problem

Mechanism synthesis divides into two families by what the optimizer is allowed
to change. In the first the topology is given and only its **dimensions** are
chosen; in the second the **topology itself** is an outcome of the optimization.
The distinction matters here because the two families answer different questions
about this engine, and because only one of them can be applied to it.

The first family — *dimensional synthesis to a kinematic criterion* — is the
dominant formulation: choose link lengths to minimise a path or motion error, or
to maximise a transmission quality index, subject to bounds on the envelope.
{cite:t}`balli2002` review the transmission angle as such an index and its use as
both objective and constraint; {cite:t}`gosselin1990` give the Jacobian-based
characterisation of the configurations where force transmission degenerates,
which is what a transmission-angle criterion is a scalar proxy for.
{cite:t}`cabrera2002` are representative of how these problems are solved: a
genetic algorithm on a bounded box, the criterion evaluated kinematically.

The conventional treatment of *this* mechanism follows that pattern, maximising a
lever-arm quality measure

$$\eta = \frac{\int_0^{2\pi} M_r\,\mathrm{d}\theta_1}
              {2(\mathrm{STE} + \mathrm{STC})\,\bar P} \tag{2.1}$$

subject to bounds on the two envelope dimensions $H$ and $B$, over the eleven
linkage variables. Three limitations follow, and they belong to the formulation
rather than to its solution.

It is **multi-objective without an exchange rate**. $\eta$, $H$ and $B$ compete
and nothing prices one against another, so the formulation yields a Pareto front
and never a design. Weighted sums, $\varepsilon$-constraint and moving limits each
produce *a* point, but the weights are the designer's rather than the physics'.

Its central quantity is **not an efficiency**. With no friction in the model the
virtual-work identity makes $\int M_r\,\mathrm{d}\theta_1 \equiv \int
P\,\mathrm{d}\lambda$ at every crank angle, so (2.1) is a ratio of two provably
equal works: a kinematic quality measure in which nothing is lost. This is a
general property of quasi-static transmission indices rather than a defect
peculiar to this one — the literature treats them as *proxies*, and a proxy is
only as good as its correlation with the quantity of interest.

It **cannot see the parts**. Nothing in a quasi-static formulation determines a
cross-section, so the mechanism has no mass, no inertia loads and no speed
dependence.

The alternative adopted here — carry the analysis through to the quantity the
application scores — is standard practice in vehicle-level MDO rather than a
novelty, and the ultra-efficiency-vehicle literature supplies the missing link:
{cite:t}`gechev2020` model an Eco-marathon prototype over its actual track and
show that achievable consumption is set jointly by the powertrain map and the
driving strategy, so neither can be scored without the other. That literature
optimises the *strategy* at fixed hardware; what is optimised here is the
*hardware*, with the strategy solved to optimality inside each evaluation.

### Synthesising the topology, not only its dimensions

The second family generates the linkage itself. The line of work begun at Seoul
National University by {cite:t}`kim2007spring` is the one relevant here: the
design domain is filled with rigid blocks joined by **zero-length springs of
variable stiffness**, and the stiffnesses are the design variables. A spring
driven stiff welds two blocks into one rigid member; a spring driven soft becomes
a revolute joint; so a single continuous parameterisation covers every linkage the
domain can hold, and the topology, the joint positions and the link lengths are
determined in one solve rather than assumed. The relaxation is the same device
§2.4 describes for catalogue choices — a discrete present-or-absent decision made
continuous and then penalised back to its extremes — applied to joints instead of
to stock. The family has since been extended, in a form directly
relevant to a mechanism whose 2:1 relation is carried by a gear pair, to planar
**gear-linkage** mechanisms {cite:p}`yim2019gearlinkage`; and a recent
spring-connected *link* model {cite:p}`tran2024slm` reaches the same answers with
far fewer design variables.

Two things about it are worth stating precisely, because the springs invite a
misreading. **The result is a rigid-body linkage, not a compliant mechanism.**
The springs are a modelling device that the penalisation drives out; what is
delivered is a set of rigid members and revolute joints. That is the reason the
method belongs in this review at all. Compliant and soft mechanisms, which
deliver their motion by elastic deformation, are *not* adapted to this
application and are not considered anywhere in this study: the linkage carries
kilonewton gas loads and joint reactions of 6 to 12 kN (§5.1), it must survive
fatigue at a thousand cycles a minute, and an extended-expansion cycle depends on
a kinematically exact 2:1 relation between the shafts — a relation that a
deforming member does not hold, and whose error goes straight into the
top-dead-centre gap that §5.2 shows this mechanism can least afford.

This study belongs to the first family: the EXlink topology is given, and eleven
dimensions and a gear pair are chosen. That is a deliberate restriction rather
than an oversight, and it has a cost the results make explicit — §5.4 measures
what the topology is worth against a slider-crank and finds it cannot say
whether a *better* topology exists, because two topologies establish a contrast
and not a trend. The second family is the natural instrument for that question,
and §5.6 applies it: a spring-connected synthesis run on this engine's own
requirements, with no mechanism assumed.

## 2.3 Architectures, derivatives, and thin feasible sets

Once mass, inertia and structural sizing enter, the model is coupled.
{cite:t}`cramer1994` established the vocabulary — multidisciplinary feasible
(MDF), individual discipline feasible (IDF), all-at-once — and
{cite:t}`martins2013architectures` give the modern survey, including the bi-level
family of which {cite:t}`bliss2000` is canonical. MDF converges an MDA at every
optimizer iteration, so every evaluated point is physically consistent and the
design space stays small, at the cost of the inner iteration and of
differentiating through it; IDF promotes the coupling variables to design
variables with consistency constraints, at the cost of a design space that grows
by the dimension of the coupling. The survey's trade favours IDF when the MDA is
expensive relative to the optimizer's handling of a larger space; that trade is
decided here by a single number, and §4.3 shows it is not close.

Whether a gradient method is admissible at all depends on the derivatives.
{cite:t}`martins2013derivatives` unify the available routes and
{cite:t}`sobieski1990` gives the global sensitivity equations that extend them to
coupled systems. Finite differences are the default and are *wrong* here rather
than merely inaccurate, several constraints being maxima over the crank
revolution whose maximiser switches as the design moves. Complex-step
{cite:p}`martins2003complexstep` removes subtractive cancellation but not the
switching problem; algorithmic differentiation {cite:p}`griewank2008` is
applicable at the cost of a dependency and of differentiating through the
fixed-point solver unless taught not to. What is used here is analytic
propagation combined with the envelope theorem — classical, {cite:t}`danskin1966`
for max-functions and {cite:t}`milgrom2002` in general — which removes the
switching problem rather than mitigating it (§4.2).

The geometry of the feasible set then selects the optimizer. Derivative-free
direct search ({cite:t}`powell1994`, Nelder–Mead) and population methods
({cite:t}`deb2002`) proceed by sampling, and a sampling method evaluates points
that lie off a measure-zero set with probability one, so it has no feasible point
to improve from; §4.1 gives the measurement. Penalty and augmented-Lagrangian
methods {cite:p}`nocedal2006` are possible but inherit the conditioning. SQP with
exact gradients {cite:p}`kraft1988` moves along the manifold and is what is used.
This is precisely why the synthesis literature above can use genetic algorithms
and this study cannot: {cite:t}`cabrera2002` optimise over a box in which every
sampled point is feasible. The same property restricts the escape from local
optima, since multistart {cite:p}`rinnooykan1987`, basin hopping and the global
methods all assume the ability to generate a feasible starting point — leaving
restarts constructed *on* the manifold (§4.6).

## 2.4 Discrete catalogue choices

The gear pair is discrete, so the problem is an MINLP; {cite:t}`belotti2013`
survey the field. Exhaustive enumeration is exact on a small catalogue but
offers no bound and no stopping criterion; relax-and-round is cheap but its
rounded point generally leaves the feasible set; branch and bound needs a
relaxation whose bound is meaningful. The two decomposition methods share a
structure — a continuous NLP for fixed integers, a MILP master accumulating
cuts — and differ in what the master is built from: optimal-value sensitivities
for generalized Benders {cite:p}`benders1962,geoffrion1972`, linearisations of
$f$ and $g$ for outer approximation {cite:p}`duran1986,fletcher1994`.
{cite:t}`fletcher1994` also supply the treatment of *infeasible* sub-problems
that a thin feasible set makes routine here. §4.4 chooses between them on which
quantities this problem can supply reliably.

A separate strand treats catalogue choices by continuous relaxation with a
penalty driving the relaxation back to a vertex — DMO {cite:p}`stegmann2005` and
SFP {cite:p}`bruyneel2011`, both descended from the SIMP penalty of
{cite:t}`bendsoe1989`. That machinery suits a large catalogue and a cheap
sub-problem; here the catalogue has a handful of entries and each sub-problem
solve is the expensive object, so it would convert an exactly solvable outer
problem into an approximately solvable one.

## 2.5 Design under tolerance

The requirements of §3.1 are equalities that no manufactured part can meet
exactly, so they must be relaxed to bands and the bands justified. Safety factors
applied to the nominal constraint are simple, but the implied reliability is
unknown and varies from constraint to constraint. Worst-case or fixed-margin
robust design, $g + k\sigma_g \le 0$, surveyed by {cite:t}`beyer2007`, computes a
margin rather than a probability, and is a reliability statement only under
independence — which §4.5 shows fails here, correlations reaching $\pm 1$ because
the constraints are analytic functions of the same eleven dimensions.
Reliability-based design optimization constrains a probability of failure: the
first-order reliability method of {cite:t}`hasofer1974`, made algorithmic by
{cite:t}`rackwitz1978`, gives $\beta = -g/\sigma_g$ at the cost of one gradient,
and its use *inside* an optimization is treated by {cite:t}`tu1999`, by
{cite:t}`du2004`, and in the survey of {cite:t}`valdebenito2010`. Because the
failure event here is a *system* event, the union probability is needed;
{cite:t}`ditlevsen1979` gives the classical bounds and {cite:t}`genz1992` the
transformation that computes the correlated orthant directly. Tolerance classes
are taken from {cite:p}`iso286`.

## 2.6 Position of this work

This study is conventional in its machinery and unconventional only in what the
machinery is applied to. The architecture comparison follows
{cite:t}`martins2013architectures`, the analytic sensitivities
{cite:t}`martins2013derivatives`, the outer loop {cite:t}`duran1986` and
{cite:t}`fletcher1994` as implemented in GEMSEO {cite:p}`gallard2018`, and the
reliability treatment FORM {cite:p}`hasofer1974,rackwitz1978` with a system
probability {cite:p}`ditlevsen1979,genz1992`. Nothing in §4 is a new algorithm.

What the literature leaves open is fivefold. The extended-expansion literature
sizes no parts and so cannot price the mechanism against the cycle it enables.
The dimensional-synthesis literature optimises a transmission quality that the
virtual-work identity shows to be a proxy with no loss in it, and does so on a box
where sampling methods are admissible. The topology-synthesis literature
generates linkages against a *kinematic* target — a path, a motion, a
timing — and stops there, so it too determines no cross-section and prices no
part; a spring-connected block model asked for this engine would return a
mechanism, not a range. The Eco-marathon literature optimises the strategy at
fixed hardware. The RBDO literature treats the bounds as given data rather than
asking which of them the mechanism is able to hold, and at what price in the
objective. The contributions listed in §1.3 follow from joining these
strands rather than from extending any one of them, and the problem is
deliberately small enough that every claim is checked against a closed-form or
sampled reference (§4.9).
