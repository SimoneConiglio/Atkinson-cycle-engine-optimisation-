# Tolerance and reliability

The central finding is about conditioning, so a design chosen for its nominal performance is exactly what a tolerance study exists to catch.

## What manufacturing does to a finished design

The central finding is about *conditioning* — the mechanism sits near a singularity — and a
design chosen for nominal performance in a badly conditioned region is exactly what a
tolerance study exists to catch. Presenting a deterministic optimum without one would be
negligent.

Tolerances are ISO 286 IT grades, not invented numbers: `i = 0.45·D^(1/3) + 0.001·D` µm, with
IT8 at 25i for a machined member. Errors propagate two ways — **first order from the exact
Jacobians**, so a full assessment costs one extra Jacobian evaluation, and **Monte Carlo** to
check the linearisation, which is precisely what should be distrusted near a singularity.

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

**`g ≤ 0.01 mm` cannot be held.** The dimensions producing the top-dead-centre gap are held to
±0.011–0.031 mm at IT8, and combine to give `g` a standard deviation of 0.013 mm — larger than
the constraint band itself. Process capability is **0.11** against an industrial target of
1.33, and two thirds of nominally conforming builds violate it.

Scanning the IT ladder settles what to do about it. Holding `g` would need a tolerance unit
multiple of **1.25i**, below the tightest grade in the table. *No machining grade fixes it.*
This is a defect in the specification, not in any design that meets it, and the remedy is a
shim at assembly or a relaxed bound — not a better optimizer. Every other constraint is
comfortable.

First order overestimates σ by up to 80 % here, so it is **conservative**, not optimistic. Worth
stating: the opposite would make first-order robust design unusable in this region.

---

## Reliability in the formulation

`robustness.tolerance_report` measures what manufacturing does to a *finished* design — the
wrong end of the process, since it can only tell the optimizer it was wrong. But the obvious
repair is also wrong, and it is worth recording why, because an earlier version of this
package shipped it.

A fixed margin `g + k·σ_g ≤ 0` applied to each constraint separately is a reliability
statement **only if the constraints are independent**. They are emphatically not: all eight
are functions of the same eleven dimensions. Measured, the correlation reaches 0.94, and
exactly −1 for the two sides of a relaxed equality:

```
             rod   comp   gap   side  st_up st_lo  r_up  r_lo
  rod_angle   1.0  -0.28 -0.56 -0.38 -0.18  0.18 -0.15  0.15
  tdc_gap   -0.56   0.62   1.0 -0.05  0.72 -0.72  0.71 -0.71
  stroke_up -0.18   0.78  0.72  0.02   1.0  -1.0  0.94 -0.94
```

Demanding all eight hold at 3σ *simultaneously* is far stronger than demanding the system be
reliable at 3σ, and it pays for that strength by rejecting designs that are actually
acceptable. It also never computes a probability of anything.

What is computed now is the quantity actually wanted — the probability that **any**
requirement is missed:

> `P_f = 1 − Φₙ(β; ρ)`,  `βᵢ = −gᵢ/σᵢ`,  `ρᵢⱼ = ∇gᵢᵀΣ∇gⱼ / (σᵢσⱼ)`

a first-order (FORM) reliability index per constraint, combined through the multivariate
normal orthant with the correlation kept. Every gradient is the exact one, so it costs one
Jacobian evaluation and can be evaluated at every optimizer iteration;
`FailureProbabilityDiscipline` exposes it as a GEMSEO constraint, so the requirement becomes
a single `P_f ≤ target` rather than eight separate margins.

The correlation does not always help, which is the point of keeping it. Here the two largest
contributors are *anti*-correlated, so the system probability (0.645) comes out **above** the
independent estimate (0.563), not below.

**Validated against sampling, and honest about where it is weak.** At system level FORM agrees
well — 0.645 against 0.653 from 4000 Monte Carlo samples. Per constraint it does not: the
top-dead-centre gap is strongly nonlinear and FORM under-predicts it, 0.42 against 0.54. So
the sampling estimate is the reference and the first-order one is what is cheap enough to sit
inside an optimization. `gemseo-umdo` supplies both routes properly — a `Sampling` formulation
with the `Probability` statistic, and a gradient-based `TaylorPolynomial` — with the
manufacturing scatter declared through `uncertain_design_variables` as `x = nominal + u`.

## What reliability says that the deterministic study could not

At the specified bounds, `COUPLED_DESIGN` has a **64.5 %** chance of missing at least one
requirement, essentially all of it the top-dead-centre gap. Inverting the reliability relation
answers *how much* the specification must give, with a number rather than "more": the gap
needs a bound of **0.054 mm** for a 10⁻³ target, against the 0.01 mm specified.

Relax it to that, and the problem moves rather than disappearing:

```
  constraint                g      sigma    beta     P(fail)
  tdc_gap            -0.04726    0.01531    3.09   1.010e-03    <- fixed
  stroke_lower       -0.01988    0.02910    0.68   2.473e-01    <- now binding
  system P(fail) 2.51e-01,  beta 0.67
```

The stroke band is *also* narrower than the scatter it is meant to represent, **and the design
sits off-centre inside it** — the residual is −0.030 of a ±0.05 band, so it spends most of its
margin on one side. Re-centring alone would take β from 0.68 to 1.7; reaching 10⁻³ needs a
band of ±0.09 mm.

Neither of those is something a deterministic optimum can see: it holds every constraint, at a
point from which half the builds fall out. That is what putting reliability in the formulation
buys, and it is a different answer from the fixed margin, not a rescaling of it.

Install it with `pip install exlink-opt[minlp]`; `exlink.minlp` is the only module that needs
the plugin and is deliberately not imported from the package root.

There is a second trap in the enumeration itself. Ranking candidate lattice points by
*distance* from the requested `I` is the obvious thing and the wrong thing: the nearest points
are reached with the **smallest** modules, and a small module needs a wide face to carry a
given tooth load. At the near-singular design — which puts 9 kN through the mesh, six times
what the backed-off one does — every point on the immediate lattice is unbuildable, and
reaching a workable pair means moving the centre distance by 13 %. `buildable_neighbours`
ranks by what can carry the load first and distance second.

That is one more way the singularity makes itself felt, and the geometric problem, where the
gears are two continuous radii, cannot see it at all.

Getting the module itself wrong is instructive, and it cost a debugging session worth
recording. Left free,
the module choice makes the objective a step function of `I`: the lightest workable module
changes at a threshold and the range jumps 40 km/L across it. A central difference straddling
that threshold returns a gradient of **3.7 × 10⁵** against a true gradient of order 10³. SLSQP,
handed a quadratic subproblem built from two constraint gradients that are both quantisation
noise, rejected it as *"inequality constraints incompatible"* and stopped at the starting
point having evaluated nothing. Pinning the pair — which also removes `I` from the design
space, leaving ten continuous variables — is both the fix and the correct algorithm.

---

Next: [Whether the answer is even a local optimum](search.md)
