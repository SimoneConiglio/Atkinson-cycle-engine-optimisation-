# 2. State of the art

Each subsection states a feature of the problem, reviews the methods available
for it, and says which is applicable here and why. The choices themselves are
made in §3.

## 2.1 Formulating a mechanism-design problem

The conventional treatment of this mechanism maximises a lever-arm quality
measure

$$\eta = \frac{\int_0^{2\pi} M_r \, d\theta_1}{2 (STE + STC)\, \bar P}$$

subject to bounds on the two envelope dimensions $H$ and $B$, over the eleven
linkage variables. Three things follow from that choice.

It is **multi-objective without an exchange rate**. $\eta$, $H$ and $B$ compete
and nothing prices one against another, so the formulation yields a Pareto front
and never a design. Weighted sums, $\varepsilon$-constraint and moving limits all
produce *a* point, but the weights are the designer's, not the physics'.

Its central quantity is **not an efficiency**. With no friction in the model the
virtual-work identity makes $\int M_r \, d\theta_1 \equiv \int P \, d\lambda$ at
every crank angle, so $\eta$ is a ratio of two provably equal works: a kinematic
quality measure in which nothing is lost.

It **cannot see the parts**. Nothing in a quasi-static formulation determines a
cross-section, so the mechanism has no mass and no inertia loads.

The alternative used here is to carry the analysis through to the quantity the
application scores. That is standard practice in vehicle-level MDO and removes
all three objections at once, at the cost of requiring the disciplines §3.2
introduces.

## 2.2 Solving a problem whose feasible set is thin

For a feasible set defined partly by equalities, the available families are:

| family | examples | applicability here |
|---|---|---|
| derivative-free direct search | COBYLA, Nelder–Mead | **no** — proceeds by sampling; see §3.3 |
| population methods | differential evolution, NSGA-II, GA | **no** — same reason |
| penalty / augmented Lagrangian | external penalty, AL | possible, but inherits the conditioning |
| SQP with exact gradients | SLSQP, NLOPT-SLSQP | **yes** — moves along the manifold |
| interior point | IPOPT | yes in principle; not evaluated here |

The distinction is not efficiency. A sampling method evaluates points that lie
off a measure-zero set with probability one, so it has no feasible point to
improve from; §3.3 gives the measurement.

## 2.3 Obtaining derivatives

Four routes, in descending order of cost per accuracy:

**Finite differences** are the default and are *wrong* on this problem, not
merely inaccurate: several constraints are maxima over the crank revolution, and
the sample attaining the maximum switches as the design moves. A difference
quotient taken across the switch is meaningless — measured at 25 % error on the
side-load ratio at a $10^{-4}$ mm step.

**Complex-step** removes subtractive cancellation but not the switching problem,
and requires the whole chain to be complex-safe.

**Algorithmic differentiation** is applicable and would avoid hand derivation, at
the cost of a dependency and of differentiating through the fixed-point solver
unless it is taught not to.

**Analytic derivatives with the envelope theorem** are used here. The chain is
closed form, so forward-mode propagation is direct; and for a maximum attained
at $\theta^*$ the derivative is the partial derivative evaluated there, because
the term through the moving maximiser carries $\partial f/\partial\theta = 0$.
That removes the switching problem rather than mitigating it. §3.4 develops it.

## 2.4 Architectures for a coupled problem

For two disciplines coupled through shared variables, the standard architectures
are MDF, IDF and their bi-level relatives.

**MDF** converges an MDA at every optimizer iteration, so every evaluated point
is physically consistent and the design space stays small; the cost is the inner
iteration.

**IDF** promotes the coupling variables to design variables with consistency
constraints, so no inner iteration runs; the cost is a design space that grows by
the dimension of the coupling.

The textbook trade favours IDF when the MDA is expensive. That trade is decided
here by a single number — the dimension of the coupling — and §3.5 shows it is
not close.

## 2.5 Mixed-integer nonlinear programming

The gear pair is discrete. The available families are:

| family | note |
|---|---|
| exhaustive enumeration | exact on a small catalogue; no bound, no stopping criterion |
| relax-and-round | cheap; the rounded point generally leaves the feasible set |
| branch and bound | needs a relaxation whose bound is meaningful |
| generalized Benders (GBD) | master built from optimal-value sensitivities |
| outer approximation (OA) | master built from linearisations of $f$ and $g$ |

GBD and OA share a decomposition and differ in what the master is built from.
OA is chosen in §3.6, for reasons that turn on which quantities this problem
can supply reliably.

## 2.6 Design under uncertainty

Three levels, in increasing fidelity and cost:

**Safety factors** applied to the nominal constraint. Simple, but the implied
reliability is unknown and varies constraint to constraint.

**Worst-case or fixed-margin robust design**, $g + k\sigma_g \le 0$. Computes a
margin, not a probability, and is a reliability statement only if the
constraints are independent — §3.7 shows they are not, with correlations
reaching $\pm 1$.

**Reliability-based design optimization (RBDO)**, constraining a probability of
failure. The probability may be estimated by sampling (accurate, expensive per
iteration) or by a first-order reliability method (cheap, approximate). §3.7
uses FORM inside the optimization with sampling as the reference.

## 2.7 Escaping local optima

Multistart, basin hopping, and global methods (DIRECT, simulated annealing,
population methods) all require the ability to *generate a feasible starting
point*. §3.3 shows that uniform generation cannot do so here, which restricts
the available options to restarts constructed on the feasible manifold.

---

Next: [3. Methodology](methodology.md)
