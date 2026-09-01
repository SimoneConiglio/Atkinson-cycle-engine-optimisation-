# Local optima and multistart

Every gradient run here is a single solve from a single point. Closing that gap gives a mixed answer.

Every gradient run in this repository is a single SLSQP solve from a single starting point.
On a problem this nonconvex that is a real gap, and closing it gives a mixed answer worth
recording.

**GEMSEO's `MultiStart` cannot be used, and not for want of budget.** It draws restarts from an
LHS over the design box. The feasible set here contains two *equality* constraints, `STE = 74`
and `ε = 16`, so it is a codimension-two manifold of measure zero and uniform sampling hits it
with probability zero. Measured:

| sampled from | analysable | feasible |
|---|---|---|
| global box | 5.8 % | **0 / 4000** |
| ±50 % around a feasible design | 33.5 % | **0 / 4000** |
| ±10 % around a feasible design | 75.8 % | **0 / 4000** |

Run anyway, `MultiStart` returns its starting point unchanged.

**What works is restarting on the manifold.** `scenarios.multistart` perturbs the incumbent,
projects the perturbation back onto the two equalities with `project_onto_equalities` — a
deterministic Newton step from the analytic Jacobians — and lets SLSQP restore the
inequalities from there. Even so only 1–2 restarts in 10 reach feasibility, which is itself a
measure of how thin the set is.

## The efficiency optimum was local

`GRADIENT_DESIGN` came from one SLSQP run and reports **30.91 %**. Restarting finds **36.99 %** —
so it was a local optimum, by six points.

That is not a better engine, and the distinction matters. The better point is **443 mm** tall
against 320, and it sits exactly on the `g` bound: feasible at 720 crank samples, not at 1440.
The single-objective efficiency problem is unbounded in mechanism size — which is the stated
reason this package does not use it as the real objective — so a stronger search simply
exploits that harder. Multistart found a better answer to the wrong question.

## On the range problem the question is open

Range is bounded, so it is the one where a better local optimum would matter. There,
**0 of 6 restarts reached feasibility** in 40 minutes at 25 SLSQP iterations each. So this
does not show `COUPLED_DESIGN` is globally optimal; it shows the experiment was inconclusive
at an affordable budget, and that is how it is reported. Establishing it would need either a
much larger budget per restart or a feasibility-restoration phase ahead of each one.

---

Next: [What is checked against independent results](validation.md)
