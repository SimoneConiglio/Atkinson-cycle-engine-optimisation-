# Exact derivatives

The feasible set is a sliver, so a derivative-free search is not merely slow but inapplicable. Every derivative that matters is available in closed form.

The feasible set here is a **sliver**. At the reference design the two equality
constraints leave a band 0.1 mm wide on `STE` and 0.1 wide on `ε`, inside an 11-dimensional
box with sides of tens of millimetres, while `W` and `γ` sit within 0.4 % and 7 % of their
bounds. A derivative-free method cannot work in that — COBYLA returns its starting point
unchanged after 313 s and 120 evaluations, whatever the budget.

The whole analysis chain is closed form, so its derivatives are too. `exlink.jacobian`
propagates them forward through the same chain the kinematics evaluates, and gets the
extremum-based metrics from the **envelope theorem**: for `max_θ f(X, θ)` attained at `θ*`,
the derivative is `∂f/∂X` evaluated there, because the term through the moving maximiser
carries `∂f/∂θ = 0`.

That matters more than "faster". These metrics are maxima over the crank angle, so the
sample attaining them *switches* as the design moves, and a difference quotient taken across
the switch is simply wrong — on `γ` at a 1e-4 mm step, 25 % wrong. GEMSEO's own
`check_jacobian` passes on all of them; `tests/test_jacobian.py` pins both the agreement and
that failure mode.

| | COBYLA | SLSQP + finite differences | SLSQP + exact gradients |
|---|---|---|---|
| moved from start | 0 mm | 32 mm | 38 mm |
| time | 313 s | 15 s | **4 s** |

Started from the historical baseline design, the gradient-based search reaches
**η = 30.77 %**, feasible at every resolution from 720 samples up:

| | η | `H` | `W` | `g` | feasible |
|---|---|---|---|---|---|
| historical baseline (`PUBLISHED_DESIGN`) | 35.62 % | 283.2 mm | 0.9892 | 8.52 mm | **no** |
| augmented Lagrangian (`REFINED_DESIGN`) | 27.87 % | 238.5 mm | 0.9811 | 0.0060 | yes |
| SLSQP + exact gradients (`GRADIENT_DESIGN`) | **30.77 %** | 319.8 mm | 0.9850 | 0.0095 | yes |

Read that comparison carefully. The baseline's 35.62 % is not a real result: that design
violates five constraints when re-analysed, most glaringly `g = 8.5 mm` against a 0.01 mm
bound, and efficiency is unbounded above once the constraints are dropped. And nothing limits
the envelope in this single-objective form, so the gradient result's extra three points over
the augmented Lagrangian are bought partly with size — the mechanism grows to `H = 320 mm`.
What the comparison does show is that a derivative-free polish stops well short of the
efficiency available at comparable feasibility.

## Through the MDA, too

`exlink.dynamics_jacobian` differentiates the coupling itself, so GEMSEO assembles the
coupled derivative from each discipline's *local* Jacobian rather than differencing the whole
fixed point. Three ideas carry it:

- **The spectral operator is linear.** Accelerations are `Ω²·D²r`, so `da/dp = Ω²·D²(dr/dp)`
  — the same operator applied to the kinematic derivative arrays. No new mathematics.
- **The 18×18 solve** gives `dx/dp = A⁻¹(db/dp − dA/dp·x)`, reusing the factorisation the
  forward solve already made. `A`'s entries are moment arms, so `dA/dp` follows from the
  position derivatives.
- **The sizing bisection is never differentiated.** The diameter is defined implicitly by
  `U(d, N, M) = 1`, so the implicit function theorem gives `dd/dN = −(∂U/∂N)/(∂U/∂d)`.

Verified piece by piece against converged central differences — mass properties and
accelerations to round-off, the equilibrium solve to 5e-7, the internal loads to 4e-7, and
the sizing by directional derivatives — plus GEMSEO's `check_jacobian` on the dynamics
discipline.

The same finite-step trap appears here and is pinned in `tests/test_dynamics_jacobian.py`:
perturb the loads by 1e-4 and the crank angle attaining the fatigue extremum hops across a
near-flat minimum, putting the difference quotient tens of percent out. The implicit function
theorem doesn't care.

**What this buys.** Minimising total moving mass at 1000 rpm, from the augmented-Lagrangian
design, subject to every constraint and a 25 % efficiency floor:

| | COBYLA | SLSQP, FD through the MDA | SLSQP, analytic |
|---|---|---|---|
| result | did not move | did not finish | **1.039 → 0.234 kg** |
| cost | 120 evals, 313 s | timed out at 560 s | **40 evals, 148 s** |

| | refined (quasi-static optimum) | coupled optimum |
|---|---|---|
| total moving mass | 1.039 kg | **0.234 kg** |
| peak bearing load | 12 629 N | 7 504 N |
| `H` | 238.5 mm | 205.7 mm |
| `W` | 0.9811 | 0.9372 |
| `η` | 28.20 % | 25.00 % |

Four times lighter, a third off the bearing load, a smaller envelope, for three points of
efficiency — and it got there by moving off the transmission-angle singularity, which is
where the mass was going. Ships as `reference.COUPLED_DESIGN`, feasible at every crank-angle
resolution from 360 to 2880 samples.

**Still differenced:** `η`, `H`, `B` and the clearance in the analysis discipline. All are
smooth, none is tight, and `η` would additionally need the crank angle of top dead centre
differentiated, because the combustion pressure jump puts moving-boundary terms in its
integral.

---

Next: [The variables that are not continuous](discrete.md)
