# Theory

The mathematics behind `exlink`, and where this reconstruction departs from the
2015 report. Symbols follow the report; the code follows the symbols.

Units throughout: **mm**, **rad** (degrees only at the API surface), **N**,
**MPa** = N/mm², **N·mm**.

---

## 1. Parametrisation

The mechanism is a two-shaft linkage. The crankshaft `R1` carries the crank `q1`
ending at `Q`; the eccentric shaft `R2` carries `q2` ending at `D`. `R2` sits at
distance `I` from `R1` along direction `θ_r`. A pair of gears of primitive radii

    r₁ = 2I/3      r₂ = I/3      r₁/r₂ = 2

ties the two shafts. The swing rod `a` runs `Q → A`; the **trigonal link** is the
rigid triangle `A–D–E`; the piston rod `e` runs `E → P`; the piston crown `H` sits
`p = 16 mm` above `P`, on the cylinder axis `x = x₁`.

Describing the triangle by its three sides `b`, `c`, `d` would force the design
space to respect the triangle inequality, and would leave the sign of `θ_b`
undetermined. The report instead places `E` in the frame carried by `c = AD`:

    b = √(x_b² + y_b²)      θ_b = atan2(y_b, x_b)      d = √((x_b − c)² + y_b²)

Now `x_b` and `y_b` range freely over ℝ, `θ_b` carries its own sign, and the design
space is a plain box. `Design.b`, `.theta_b`, `.d` in `exlink/design.py`.

This is checked against the report's own Carnot expression
`θ_b = arccos((b² + c² − d²)/(2bc))` in `tests/test_design.py`.

---

## 2. Kinematics

Six degrees of freedom (`θ₁, θ₂, θ_T, θ_a, θ_e, λ`) and five constraints, so one
input `θ₁` fixes everything.

**Gear relation (1).** External gears turn opposite ways, at the inverse ratio of
their radii:

    θ₂ = −2θ₁ + θ_f

**Loop closures (2), (3).** Two vector chains close on themselves:

    R1 → Q → A → D → R2 → R1  = 0
    R1 → Q → A → E → P → H → R1  = 0

Projecting the first on the axes and isolating the terms in `a` and `c`:

    A = q₁ sin θ₁ − q₂ sin θ₂ + I cos θ_r
    B = −q₁ cos θ₁ + q₂ cos θ₂ + I sin θ_r        (3a)

Squaring and adding eliminates `θ_a` and `θ_T` separately, leaving their
difference `T = θ_a − θ_T`:

    A² + B² = a² + c² + 2ac cos T

    T = arccos( (A² + B² − a² − c²) / (2ac) )     (4)

**Compatibility condition (4a).** That arccosine argument must stay inside
(−1, 1). Define

    δ_c1 = max over θ₁ of |(A² + B² − a² − c²) / (2ac)|

If `δ_c1 ≥ 1` for even one crank angle, the four-bar cannot pass that angle: the
crankshaft rocks instead of turning. This is the Grashof condition for the
sub-mechanism `a, c, q₁, I, q₂`, written as something an optimizer can read.

Then, with `q = atan2(a sin T, a cos T + c)`:

    θ_T = atan2(B, A) − q          θ_a = θ_T + T          (5)

**Piston rod (6).** From the horizontal projection of the second chain:

    cos θ_e = (q₁ sin θ₁ − a cos θ_a − b cos(θ_b + θ_T) + x₁) / e

    δ_c2 = max over θ₁ of |cos θ_e|                (6a)

**Piston height (7).** From the vertical projection:

    λ = q₁ cos θ₁ + a sin θ_a + b sin(θ_b + θ_T) + e sin θ_e + p

> **Departure.** The report's equation (7) omits the constant `+ p`, so its `λ`
> is the height of the wrist pin `P`, not of the crown `H`. It is restored here
> so the two agree with the figure. Being constant it cancels out of every stroke
> and volume, so nothing downstream changes.

**Why the analytic inversion matters.** Newton–Raphson would solve the same
equations. But `δ_c1` and `δ_c2` only exist because the inversion is explicit, and
handing those two numbers to the optimizer — instead of letting the analysis
diverge — is what turns a problem full of hard failures into one with a usable
search landscape. The report makes this point, and it is the single most
important idea in the study.

**Critical configurations.** `|cos T| = 1` means `a` and `c` are parallel: the
swing rod stops working as a rod and the mechanism gains a degree of freedom.
The constraint is set at `W ≤ 0.985`, keeping `T` inside [10°, 170°]. The other
pair, `θ_e = 0` or `π`, is already excluded by the 10° rod-angle limit.

Verified in `tests/test_kinematics.py`: every link keeps its length to 1e-9 mm,
and `P` and `H` stay on `x = x₁`, over the whole revolution.

---

## 3. The Atkinson cycle

Over one crankshaft revolution `λ(θ₁)` must have **four monotone phases**: two
maxima (top dead centre, reached twice) and two *different* minima. The deeper
minimum ends expansion, the shallower one ends intake.

    TDC → deep BDC     expansion
    deep BDC → TDC     exhaust
    TDC → shallow BDC  intake
    shallow BDC → TDC  compression

    STE = λ_TDC − λ_deep        STC = λ_TDC − λ_shallow

    g = |λ_TDC,1 − λ_TDC,2|     the gap between the two top dead centres

`find_phases` counts sign changes of `dλ/dθ₁` and refines each extremum with a
parabola through its three neighbouring samples — necessary because `g` is
constrained to 0.01 mm on a 0.5° grid.

Volume and pressure, with `A_p = πΦ²/4` the bore area:

    V = V₀ + A_p (λ_TDC − λ)
    V₁ = V₀ + A_p·STC          ε = V₁/V₀

    intake, exhaust  P = P₀
    compression      P = P₀ (V₁/V)^γ
    combustion       P₃ = k·P₂,  P₂ = P₀ ε^γ     (instantaneous, at TDC)
    expansion        P = P₃ (V₀/V)^γ
    blow-down        P → P₀                      (instantaneous, at deep BDC)

The gas force on the crown uses the **gauge** pressure, `P_gas = (P − P₀)·A_p`, so
intake and exhaust are unloaded.

> **Departure.** The report writes compression as `P = P₀(V/V₁)^γ`, which would
> make pressure *fall* as the charge is compressed. The exponent sign is a typo:
> its own `P₂ = P₀ε^γ` requires `P = P₀(V₁/V)^γ`, which is what is implemented.

Given `Φ = 32 mm`, `V₀ = 3 cc` and `ε = 16`, the compression stroke is pinned at
`STC = 15·V₀/A_p ≈ 55.95 mm` against `STE = 74 mm` — the asymmetry the linkage
exists to produce.

**Designs that fail the phase test are penalised, not rejected:** `η = 0`,
`H = B = 1000`. A design whose piston goes up and down once per revolution is a
plain Otto engine; both failure modes are exercised in `tests/test_model.py`.

---

## 4. Quasi-static loads

Inertia is neglected — this is a first sizing iteration, and the masses are not
known until the parts have a shape.

**Piston.** With `P` the gas force and `θ_e` the rod angle:

    C = P / sin θ_e        rod load
    D = P cot θ_e          side load reacted by the liner

**Trigonal link.** Force balance at `A`, `D`, `E` plus the moment about `D`.
Writing `DE = b·u(θ_b + θ_T) − c·u(θ_T)` and `u_e = (cos θ_e, sin θ_e)`:

    DA × F_A + DE × F_E = 0
    ⇒ −c·A·sin(θ_a − θ_T) − C·(DE × u_e)_z = 0

    A = −C (DE × u_e)_z / (c sin(θ_a − θ_T))

    Q_x = C cos θ_e − A cos θ_a        Q_y = C sin θ_e − A sin θ_a

> **Departure — a sign correction.** The report prints this inversion **without
> the leading minus sign**. With its sign, the computed torque disagrees with the
> principle of virtual work by a factor of about −4, and the efficiency comes out
> negative. With the sign above, agreement is exact. See §6.

Note `sin(θ_a − θ_T) = sin T` in the denominator: the swing-rod load blows up at
exactly the critical configurations condition (4a) excludes — a second,
independent reason to enforce it.

**Shafts.** With `α` the gear pressure angle (20°, the standard involute value;
the report writes `α` symbolically but never gives a number):

    T_gear = −q₂(Q_y sin θ₂ + Q_x cos θ₂) / (r₂ cos α)

    M_r = q₁ A cos(θ_a − θ₁) + r₁ T_gear cos α

The two gear torques are `r₁ T cos α` and `r₂ T cos α` with the *same* sign, so
with `ω₂ = −2ω₁` and `r₁ = 2r₂` the pair transmits no net power — checked
explicitly in `tests/test_loads.py`, since a gear pair that generated power would
silently inflate the efficiency.

> **Departure.** The report's bearing reactions `R₁ₓ = A cos θ_a + T sin α` and
> `R₁ᵧ = A sin θ_a − T cos α` drop the `θ_r` dependence of the tooth-force
> direction. They should read `sin(θ_r + α)` and `cos(θ_r + α)`, as implemented.
> These reactions feed nothing downstream, so only their own values change.

---

## 5. Objectives and constraints

**Efficiency.** The report's average mechanical efficiency:

    η = ∮M_r dθ₁ / (2(STE + STC)·⟨P⟩) = (⟨M_r⟩/⟨P⟩) · π/(STE + STC)

a ratio of two works — the torque's on the crankshaft over the gas force's on the
piston. It measures the linkage's aptitude for turning force into torque, and it
grows without bound as the mechanism grows. That unboundedness is exactly why `H`
and `B` must enter the problem.

`φ = ⟨M_r⟩/⟨p⟩` is also reported. It has the dimension of a volume; the report
prints it as "94.46 %" with no stated normalisation, so `exlink` computes it in
mm³ and does not compare against that figure.

**Envelope.** `H` along the stroke, `B` across it: the bounding box of every body
over every configuration — joints, both gear primitives, and the piston over its
full travel.

**Clearance `d`.** The trigonal link must stay 10 mm clear of the cylinder.

> **Departure.** The report states this constraint but not the construction
> behind it. `cylinder_clearance` models the liner as the half-strip
> `x ∈ [x₁ ± Φ/2]`, `y ≥ y_bottom`, with `y_bottom` the lowest point the piston
> skirt reaches, and minimises the distance from the three triangle edges over
> the revolution. It is monotone in the right direction and vanishes on contact,
> which is what the constraint needs — but its numerical value is not expected to
> match the report's.

**Final formulation.**

    l_b ≤ X ≤ u_b
    min  f(X)    = (−η, H, B)ᵀ
    s.t. c(X)    = (mra−10, W−0.985, g−0.01, 10−d, γ−0.02)ᵀ ≤ 0
         c_eq(X) = (STE−74, ε−16)ᵀ = 0

with `W = max(δ_c1, δ_c2)` and `γ = max(D)/max(P)`.

---

## 6. Verification

The force chain is not merely transcribed — it is **pinned by an independent
identity**. In a massless, frictionless, quasi-static mechanism, instantaneous
power in equals power out:

    M_r(θ₁) · ω₁ = P(θ₁) · v_piston      ⇒      M_r = −P · dλ/dθ₁

Every step of the chain — piston, trigonal link, both shafts, the gear pair —
must conspire to satisfy this at *every* crank angle. `exlink` reproduces it to
machine precision (`tests/test_loads.py`), which is what exposed the sign slip in
§4: with the report's printed sign, the two disagree by a factor of about −4.

A second, independent route: the mean torque must equal the indicated p–V loop
area divided by 2π. That is also checked, and it never touches the force chain.

Together with the rigid-link check (§2), these leave very little room for the
model to be wrong in a way the tests would not see.

---

## 7. Numerical procedure

The report's own sequence, and its counterpart here:

1. **External penalty**, `F(X) = −η + r⁻²(c_eqᵀc_eq + ⟨c⟩ᵀ⟨c⟩)`, turning the
   constrained problem into an unconstrained one — `PenalisedExlinkDiscipline`.
   Accurate only for small `r`, badly conditioned when `r` is small; that
   trade-off is why the report finishes with an augmented Lagrangian.
2. **Local solvers** (conjugate gradient, simplex) — start-point dependent, and
   the problem is strongly non-convex. `NELDER-MEAD`, `SLSQP`, `NLOPT_COBYLA`.
3. **Global search.** The report needed ≥ 550 individuals over 11 variables.
   `DIFFERENTIAL_EVOLUTION`, via `maximise_efficiency`.
4. **MOEA.** Over the full box it returned solutions *worse* than the
   gradient-based ones. What worked was shrinking the box around an already-good
   design, then merging several such local fronts into one starting population —
   `local_pareto`, then `pareto_front`.
5. **Moving limits.** Treat `H` and `B` as constraints rather than objectives and
   walk the limit down, tracing the trade-off with ordinary single-objective
   solves — `sweep_moving_limits`.
6. **Augmented Lagrangian**, the final polish — `refine`, the default
   `Augmented_Lagrangian_order_0`.

Because NSGA-II takes no equality constraints, `relax_equalities` rewrites each
as a pair of one-sided inequalities, `|residual| ≤ tol` — the same pragmatism the
report applies to `g`.

Two practical notes on the multi-objective stage, both measured:

- **`g` is the binding constraint, not the equalities.** Sampled over a box
  around a good design, `g ≤ 0.01 mm` is met by 0.1 % of analysable designs
  against 4–8 % for the two real equalities. It is a third equality in disguise,
  and `moea_targets` relaxes it so the population has somewhere to live.
- **Generations matter more than population.** Seeded in a shrunk box, 20
  generations return a front of one point whatever the population size (60 or
  80, 1200 or 1600 evaluations); 35 generations give 6–33 points across seeds.
  The population needs *time* to work into the thin feasible sheet, not more
  parallel guesses. Spend extra budget on `max_gen` first.

The moving-limit sweep reaches the same trade-off without relaxing anything,
because each step is a warm-started local solve on the full problem:
`H ≤ 236 → η = 28.11 %`, `232 → 27.98 %`, `228 → 27.21 %`, `224 → 27.93 %`,
each meeting its limit exactly. Being a sequence of independent local solves,
the curve is not perfectly monotone and individual steps can come back
infeasible. It is not free either — each step must satisfy both equalities and
all five inequalities while pushed against a limit it did not previously meet,
which takes roughly 120 augmented Lagrangian outer iterations.

A note on judging feasibility: an augmented Lagrangian converges *onto* the
active bounds, so a converged design routinely lands a few parts in 1e5 outside
one of them (`g`, in practice). `is_feasible` therefore allows GEMSEO's own
`ineq_tolerance` of 1e-4, and 0.05 on the two equalities — 0.05 mm on a 74 mm
stroke being well inside anything the linkage could be built to.

GEMSEO's default convergence tests also stop NSGA-II after about five
generations here — the objectives barely move while the population is still
hunting for feasible designs — so `pareto_front` disables them and lets the
generation budget be the budget that applies.

---

## 8. On reproducing the published solution

Re-analysed exactly as printed, the report's design table gives `g = 8.5 mm`
against its reported 0.0069, and `STE = 79.6 mm` against 73.98. Since `g` is
precisely the quantity the optimizer drove to zero, the printed table cannot be
the design that produced the reported properties.

Rounding does not explain it: the table is printed to four significant figures,
and perturbing each variable by its rounding half-width moves `STE` by less than
0.2 mm. The design sits at `W = 0.982`, a hair from the singularity at `W = 1`,
where the piston motion is extremely sensitive to the link lengths — which is
also why a 2015 MATLAB script and a 2026 Python one would not agree to four
digits even given identical formulae.

What *does* reproduce is the physics and the workflow. Running the report's own
final step — the augmented Lagrangian, started from the published table — lands
on a fully feasible design at **η = 27.87 %** against the reported 27.76 %, with
`W`, `g`, `STE`, `ε` and `γ` all matching closely (see the table in the README).
That design ships as `exlink.reference.REFINED_DESIGN`.


---

## 9. Sizing the parts, and the coupling that creates

The report stops before this deliberately:

> to have the masses of the pieces we have to know their shape so we should have
> a first design, so those passages are for another iteration.

That iteration is `exlink.dynamics`, `exlink.sizing` and `exlink.coupled`. It
closes a loop the quasi-static study does not have:

```
    section diameters ──▶ member masses ──▶ inertia forces ──▶ internal loads
            ▲                                                        │
            └──────────── sizing against yield, ─────────────────────┘
                          fatigue and buckling
```

Neither half can go first. That is a genuine multidisciplinary coupling, and it
has to be **solved** rather than sequenced.

### 9.1 Why the report's solution method cannot just be extended

Without inertia every rod is a **two-force member**: the forces at its two ends
are equal, opposite and collinear with the rod. That is exactly what lets the
report eliminate unknowns one body at a time, piston to crankshaft.

Add mass and that collapses. A rod with a distributed d'Alembert load has end
forces that are neither collinear nor equal, so no body can be solved before its
neighbours. The whole mechanism must be solved at once.

Counting with Gruebler over 7 links (6 moving plus ground), 8 lower pairs (7
revolutes and the piston's guide) and 1 higher pair (the gear mesh):

```
M = 3(7 − 1) − 2(8) − 1 = 1
```

One degree of freedom, so the load problem is statically determinate: **18
unknowns against 6 bodies × 3 equilibrium equations**. `exlink.dynamics.solve`
assembles and solves that 18×18 system at every crank angle.

The determinant of that matrix is the same quantity condition (4a) protects: at
a critical configuration the mechanism gains a degree of freedom, the matrix goes
singular, and the internal forces blow up. The condition number is reported so
the connection is visible rather than implied.

### 9.2 Accelerations

Every history is smooth and periodic on a uniform grid, so the Fourier
derivative is exact where a finite difference is only `O(h²)`. That matters
because accelerations are *second* derivatives: on a 0.5° grid a
finite-difference second derivative loses about six digits, enough to pollute
the inertia forces and, through them, the sizing loop. Angles that accumulate
whole turns — `theta_2 = −2 theta_1 + theta_f` — are split into ramp plus
periodic part first (`exlink.derivatives`).

Constant crankshaft speed is inherited from the report, and simplifies the
bookkeeping in two ways worth stating:

- both shafts turn at constant rate, so neither has angular acceleration and
  neither contributes an inertia couple;
- shafts and gears are concentric with their own axes, so their centres of mass
  do not move and they exert no inertia force at all. Only the offset crank
  throws do, which is why the modelled bodies are the crank arms rather than
  whole shaft assemblies.

### 9.3 Failure modes

Each link is a solid round bar whose diameter is solved for — the smallest
section satisfying all three modes over the whole revolution:

| mode | criterion |
|---|---|
| static | peak fibre stress against `S_y / n_y`; uniaxial, so von Mises is `\|σ\|` |
| fatigue | Goodman, `σ_a/S_e + σ_m/S_u ≤ 1/n_f`, with Marin-corrected `S_e = k_a k_b k_c k_d k_e · 0.5 S_u` |
| buckling | Euler on the peak compressive load, `K = 1` for links, `2` for a crank throw |

Fatigue is evaluated **per extreme fibre**, so that `σ_a` and `σ_m` are taken at
a fixed material point rather than at whichever fibre happens to be worst at each
instant. A compressive mean stress is not credited as beneficial — the Goodman
line is truncated at `σ_m = 0` — which matters here because the connecting links
swing between tension and compression every revolution.

The required diameter comes from bisection, not a closed form: the fatigue size
factor `k_b` itself depends on the diameter being solved for.

**Internal loads.** For a member spanning two joints, the force and moment at a
section a fraction `s` along it follow from the free body of `[0, s]`. Because a
rigid body's acceleration varies *linearly* along any straight line through it,
that load is linear in `s` and both integrals close in form:

```
F(s) = m [a₁ s + (a₂ − a₁) s²/2] − F₁
M(s) = −m [(Δr × a₁)_z s²/2 + (Δr × Δa)_z s³/6] + s (Δr × F₁)_z
```

**One idealisation.** The trigonal link is a single rigid part, so as a frame it
is three times statically indeterminate. It is treated as a pin-jointed
triangle: each side takes an axial force from joint equilibrium (which *is*
determinate), and bends only under its own distributed inertia as a simply
supported beam. That captures the axial load exactly and under-estimates bending
at the corners — stated here rather than hidden.

Shafts are not sized. Their bearing span is not a design variable, so any
section would be arbitrary; and being concentric with their axes they do not
feed the inertia loop at all. The bearing *reactions* are computed and
constrained.

### 9.4 Does the loop converge?

Scaling answers it. A bending-critical member needs `d ~ F^(1/3)`, so its mass
goes as `m ~ d² ~ F^(2/3)`; the inertia force it then creates is `F ~ m a`.
Composing, `m ~ (C a) m^(2/3)`. The loop gain is **sub-linear**, so a fixed point
exists at

```
m = (C a)³
```

and plain Gauss-Seidel reaches it — but slowly, and with a mass that grows as the
*cube* of the acceleration level, hence the **sixth power of engine speed**.
That is why `MDAGaussSeidel` is the right MDA here (no coupled Jacobians needed,
and the bisection has no useful derivative to give a Newton method anyway), and
why speed — which never appears in the report — becomes one of the strongest
drivers in the problem.

### 9.5 What the dynamics changes

Two results, both checked in `tests/test_dynamics.py`:

**The mean torque does not move.** At constant speed the mechanism returns to its
starting state every revolution, so its kinetic energy is unchanged and the
inertia forces do no net work. Efficiency, being a mean-torque quantity, is
therefore *untouched* by speed. Only the peaks change — and the peaks are what
size the parts. With the gas load switched off, every reaction scales as exactly
`Ω²`.

**The quasi-static optimum is the wrong answer.** The report's design sits at
`W = 0.981`, a hair from the singularity, because that is where the quasi-static
lever arm is longest. But proximity to the singularity is also what amplifies
velocities and accelerations: joint `A` sees 75× the crank pin's acceleration.
Backing off costs nothing and saves almost everything:

| swing rod | `W` | `η` | `H` [mm] | mass [kg] | peak bearing [N] |
|---|---|---|---|---|---|
| ×1.00 | 0.9811 | 28.20 % | 238.5 | 8.43 | 244 800 |
| ×0.94 | 0.9670 | 27.79 % | 227.8 | 3.41 | 43 700 |
| ×0.85 | 0.9519 | 28.26 % | 215.6 | 2.20 | 13 900 |
| ×0.80 | 0.9475 | 28.62 % | 211.6 | 2.03 | 11 300 |

Four times lighter, twenty times lower bearing load, *better* efficiency and a
smaller envelope. The near-singular design was an artefact of leaving inertia
out.

### 9.6 The coupled formulation

```
l_b ≤ X ≤ u_b
min  f(X)    = (−η, H, B, M)ᵀ
s.t. c(X)    ≤ 0     the report's five, plus
                     saturation_margin   loop ran away, no section is thick enough
                     slenderness_margin  a "rod" thicker than a third of its length
                     bearing_margin      peak reaction against its limit
     c_eq(X) = 0
```

Structural mass cannot appear in the report's formulation at all — nothing in it
determines a cross-section. It only becomes an objective once the parts are
sized, and it is the objective the dynamics most affects.

Under GEMSEO this is the **MDF** formulation: `ExlinkDiscipline` is feed-forward
from the design variables, while `DynamicsDiscipline` and `StructureDiscipline`
are strongly coupled through `diameters` one way and the internal load histories
the other. MDF wraps them in an MDA so every point the optimizer evaluates is a
converged one. `warm_start` matters more than it looks: successive design points
are close together, so starting each MDA from the previous converged sections
turns a 70-sweep solve into a handful.
