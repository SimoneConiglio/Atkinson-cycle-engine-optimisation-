# Theory

The mathematics behind `exlink`: the kinematic inversion, the cycle model, the
load analysis, the sizing criteria, the optimization formulation, and the exact
derivatives of all of it. Symbols are defined where they first appear, and the
code follows them.

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
undetermined. `E` is therefore placed in the frame carried by `c = AD`:

    b = √(x_b² + y_b²)      θ_b = atan2(y_b, x_b)      d = √((x_b − c)² + y_b²)

Now `x_b` and `y_b` range freely over ℝ, `θ_b` carries its own sign, and the design
space is a plain box. `Design.b`, `.theta_b`, `.d` in `exlink/design.py`.

The two routes to `θ_b` — this one and the Carnot expression
`θ_b = arccos((b² + c² − d²)/(2bc))` — are checked against each other in
`tests/test_design.py`.

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

Note the constant `+ p`: `λ` is the height of the crown `H`, not of the wrist
pin `P`. Being a constant it cancels out of every stroke and volume, so only the
absolute datum depends on it.

**Why the analytic inversion matters.** Newton–Raphson would solve the same
equations. But `δ_c1` and `δ_c2` only exist because the inversion is explicit, and
handing those two numbers to the optimizer — instead of letting the analysis
diverge — is what turns a problem full of hard failures into one with a usable
search landscape. It is the single most important modelling decision here.

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

Note the sense of the compression exponent: `P = P₀(V₁/V)^γ`, so pressure rises
as the charge is compressed, and at `V = V₀` it gives `P₂ = P₀ ε^γ` as it must.

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

> **The leading minus sign matters.** Drop it — an easy slip when inverting the
> moment equation — and the computed torque disagrees with the principle of
> virtual work by a factor of about −4, with the efficiency coming out negative.
> With the sign above, agreement is exact. This is what the check in §6 is for.

Note `sin(θ_a − θ_T) = sin T` in the denominator: the swing-rod load blows up at
exactly the critical configurations condition (4a) excludes — a second,
independent reason to enforce it.

**Shafts.** With `α` the gear pressure angle (20°, the standard involute value;
20° is the standard involute value and is used throughout):

    T_gear = −q₂(Q_y sin θ₂ + Q_x cos θ₂) / (r₂ cos α)

    M_r = q₁ A cos(θ_a − θ₁) + r₁ T_gear cos α

The two gear torques are `r₁ T cos α` and `r₂ T cos α` with the *same* sign, so
with `ω₂ = −2ω₁` and `r₁ = 2r₂` the pair transmits no net power — checked
explicitly in `tests/test_loads.py`, since a gear pair that generated power would
silently inflate the efficiency.

The bearing reactions carry the tooth-force direction in full,
`R₁ₓ = A cos θ_a + T sin(θ_r + α)` and `R₁ᵧ = A sin θ_a − T cos(θ_r + α)`: the
line of action rotates with the shaft axis `θ_r`, not with `α` alone.

---

## 5. Objectives and constraints

**Efficiency.** The average mechanical efficiency is defined as:

    η = ∮M_r dθ₁ / (2(STE + STC)·⟨P⟩) = (⟨M_r⟩/⟨P⟩) · π/(STE + STC)

a ratio of two works — the torque's on the crankshaft over the gas force's on the
piston. It measures the linkage's aptitude for turning force into torque, and it
grows without bound as the mechanism grows. That unboundedness is exactly why `H`
and `B` must enter the problem.

A companion measure `φ = ⟨M_r⟩/⟨p⟩`, the torque-per-unit-pressure ratio, is also
computed. It has the dimension of a volume and `exlink` reports it in mm³; it is
proportional to `η` at fixed stroke, so it adds no independent information.

**Envelope.** `H` along the stroke, `B` across it: the bounding box of every body
over every configuration — joints, both gear primitives, and the piston over its
full travel.

**Clearance `d`.** The trigonal link must stay 10 mm clear of the cylinder.

> **An idealisation.** `cylinder_clearance` models the liner as the half-strip
> `x ∈ [x₁ ± Φ/2]`, `y ≥ y_bottom`, with `y_bottom` the lowest point the piston
> skirt reaches, and minimises the distance from the three triangle edges over
> the revolution. It is monotone in the right direction and vanishes on contact,
> which is what the constraint needs — but its numerical value is not expected to
> a detailed CAD clearance check.

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
§4: with the sign dropped, the two disagree by a factor of about −4.

A second, independent route: the mean torque must equal the indicated p–V loop
area divided by 2π. That is also checked, and it never touches the force chain.

Together with the rigid-link check (§2), these leave very little room for the
model to be wrong in a way the tests would not see.

---

## 7. Numerical procedure

A sequence that works on this problem, and what implements each stage:

1. **External penalty**, `F(X) = −η + r⁻²(c_eqᵀc_eq + ⟨c⟩ᵀ⟨c⟩)`, turning the
   constrained problem into an unconstrained one — `PenalisedExlinkDiscipline`.
   Accurate only for small `r`, badly conditioned when `r` is small; that
   trade-off is why a penalty method wants an augmented Lagrangian to finish.
2. **Local solvers** (conjugate gradient, simplex) — start-point dependent, and
   the problem is strongly non-convex. `NELDER-MEAD`, `SLSQP`, `NLOPT_COBYLA`.
3. **Global search.** Over 11 variables this needs a large population — ≥ 550
   individuals in the MATLAB study this problem originates from.
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
as a pair of one-sided inequalities, `|residual| ≤ tol` — the same pragmatism
already applied to `g`, which is constrained to a tolerance rather than to zero.

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

## 8. The historical baseline, and why the problem is ill-conditioned

`exlink.reference.PUBLISHED_DESIGN` is a design vector carried over from the
earlier MATLAB study this problem originates from (see the README's *Provenance*
section). It is kept as a **starting point**, not as a result, and re-analysing it
shows why.

As printed it gives `g = 8.5 mm` against a 0.01 mm bound and `STE = 79.6 mm`
against a required 74 mm — five constraints violated in all. Rounding does not
explain that: the vector is given to four significant figures, and perturbing
each variable by its rounding half-width moves `STE` by less than 0.2 mm.

The reason is conditioning, and it generalises well beyond this one design.
The vector sits at `W = 0.982`, a hair from the singularity at `W = 1`, and there
the piston motion is violently sensitive to the link lengths. A design in that
region is not reproducible to four digits across two independent
implementations — which is a property of the optimum's *location*, not of either
implementation, and is the first hint of the result in §9.5: the quasi-statically
optimal design sits exactly where the mechanism is worst conditioned.

Started from it, an augmented Lagrangian reaches a fully feasible design at
**η = 27.87 %**, shipped as `exlink.reference.REFINED_DESIGN`. That is the
quasi-static reference point the rest of this document measures against.


---

## 9. Sizing the parts, and the coupling that creates

A quasi-static study cannot size the parts, and the obstruction is circular: the
inertia loads need the part masses, the masses need the cross-sections, and the
cross-sections need the loads. `exlink.dynamics`, `exlink.sizing` and
`exlink.coupled` close that loop. It
closes a loop the quasi-static study does not have:

```
    section diameters ──▶ member masses ──▶ inertia forces ──▶ internal loads
            ▲                                                        │
            └──────────── sizing against yield, ─────────────────────┘
                          fatigue and buckling
```

Neither half can go first. That is a genuine multidisciplinary coupling, and it
has to be **solved** rather than sequenced.

### 9.1 Why sequential elimination cannot just be extended

Without inertia every rod is a **two-force member**: the forces at its two ends
are equal, opposite and collinear with the rod. That is exactly what lets the
loads be eliminated one body at a time, piston to crankshaft.

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

Constant crankshaft speed is assumed throughout, which simplifies the bookkeeping
in two ways worth stating:

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
why engine speed — which does not appear in the quasi-static problem at all —
becomes one of the strongest drivers here.

### 9.5 What the dynamics changes

Two results, both checked in `tests/test_dynamics.py`:

**The mean torque does not move.** At constant speed the mechanism returns to its
starting state every revolution, so its kinetic energy is unchanged and the
inertia forces do no net work. Efficiency, being a mean-torque quantity, is
therefore *untouched* by speed. Only the peaks change — and the peaks are what
size the parts. With the gas load switched off, every reaction scales as exactly
`Ω²`.

**The quasi-static optimum is the wrong answer.** Maximising efficiency without
inertia drives the design to `W = 0.981`, a hair from the singularity, because
that is where the quasi-static lever arm is longest. But proximity to the singularity is also what amplifies
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
s.t. c(X)    ≤ 0     the five of §5, plus
                     saturation_margin   loop ran away, no section is thick enough
                     slenderness_margin  a "rod" thicker than a third of its length
                     bearing_margin      peak reaction against its limit
     c_eq(X) = 0
```

Structural mass cannot appear in a quasi-static formulation at all — nothing in
one determines a cross-section. It only becomes an objective once the parts are
sized, and it is the objective the dynamics most affects.

Under GEMSEO this is the **MDF** formulation: `ExlinkDiscipline` is feed-forward
from the design variables, while `DynamicsDiscipline` and `StructureDiscipline`
are strongly coupled through `diameters` one way and the internal load histories
the other. MDF wraps them in an MDA so every point the optimizer evaluates is a
converged one. `warm_start` matters more than it looks: successive design points
are close together, so starting each MDA from the previous converged sections
turns a 70-sweep solve into a handful.


---

## 10. Derivatives

The feasible set is a **sliver**. At the reference design the two equalities leave
a band 0.1 mm wide on `STE` and 0.1 wide on `ε`, inside an eleven-dimensional box
with sides of tens of millimetres, while `W` and `γ` sit within 0.4 % and 7 % of
their bounds. A derivative-free method cannot work in that: COBYLA returns its
starting point unchanged after 313 s, and loosening the bands does not help.

The whole analysis chain is closed form, so its derivatives are too.

### 10.1 The two ideas

**Forward-mode chaining.** Each intermediate carries `d/dX` as an
`(n_angles, 11)` array, differentiated straight through the closed forms of §2.
One pass produces all eleven components.

**The envelope theorem.** Most metrics are extrema over the crank angle. For a
maximum attained at `θ*`,

```
d/dX max_θ f(X, θ) = ∂f/∂X |_{θ*}
```

because the term through the moving maximiser carries `∂f/∂θ = 0`. No derivative
of the *location* is needed.

That second point is not a convenience — it removes a failure mode. Because these
are maxima, the sample attaining them **switches** as the design moves, and a
difference quotient taken across the switch is wrong: on `γ` at a 1e-4 mm step,
25 % wrong. Both `tests/test_jacobian.py` and `tests/test_dynamics_jacobian.py`
pin that behaviour deliberately.

The extremum does not sit on a grid point, so the partial derivative is
interpolated to the parabolically refined location with a three-point Lagrange
quadratic — matching the order of the refinement that located it. Linear
interpolation leaves `STE` a factor of fifty worse.

### 10.2 Through the coupling

`exlink.dynamics_jacobian` differentiates the MDA itself, so GEMSEO assembles the
coupled derivative from local Jacobians instead of differencing the fixed point.

| piece | route |
|---|---|
| accelerations | `D²` is linear, so `da/dp = Ω² D²(dr/dp)` |
| mass properties | direct; `m = ρ(πd²/4)L`, centres of mass mass-weighted |
| 18×18 solve | `dx/dp = A⁻¹(db/dp − dA/dp·x)`, reusing the factorisation |
| internal loads | closed form; the trigonal truss via `dz = (MᵀM)⁻¹Mᵀ(dr − dM z)` |
| sizing bisection | implicit function theorem on `U(d, N, M) = 1` |

Verified against converged central differences: mass properties and accelerations
to round-off, the equilibrium solve to 5e-7, the internal loads to 4e-7, the
sizing by directional derivatives, and GEMSEO's own `check_jacobian` on both
disciplines.

### 10.3 What it changes

| problem | derivative-free | with gradients |
|---|---|---|
| maximise `η` from the published table | — | **η = 30.77 %** (vs 27.76 % reported) |
| minimise mass, coupled, 1000 rpm | COBYLA moved 0 in 313 s | **1.039 → 0.234 kg** in 148 s |

The coupled optimum reaches a quarter of the mass by moving `W` from 0.981 to
0.937 — off the singularity, which is exactly where the accelerations, and so the
sections, were coming from. It also comes out *smaller* (`H` 238.5 → 205.7 mm),
for three points of efficiency.

Left differenced, deliberately: `η`, `H`, `B` and the clearance. All smooth, none
tight, and `η` would additionally need `dθ_TDC/dX`, because the combustion
pressure jump puts moving-boundary terms in its integral.

---

## 11. Why `η` is not an efficiency

Section 5 defines

$$\eta = \frac{\int_0^{2\pi} M_r \, d\theta_1}{2 (STE + STC)\, \bar{P}}$$

and calls it an average mechanical efficiency. It is not one, and the reason
matters for everything that follows.

Section 6 verifies the force chain against the virtual-work identity
$M_r = -P \, d\lambda/d\theta_1$, which holds at *every* crank angle to machine
precision. With no friction in the model, the work the torque does on the
crankshaft is therefore identically the work the gas does on the piston. Nothing
is lost. A quantity that is a ratio of two provably equal works cannot be an
efficiency.

What $\eta$ actually measures is *kinematic quality*: how much **mean** torque a
given piston motion produces, normalised by the stroke. A motion that spends
its stroke where the gas pressure is high scores well; one that wastes it scores
badly. That is a useful figure of merit and it is the right one for the
geometric problem. But it grows without bound as the mechanism grows — a longer
lever arm always converts more — which is exactly why §5 has to add $H$ and $B$
as competing objectives, and why the resulting three-objective problem has no
solution, only a front.

### 11.1 The losses that are really there

Every quantity needed for a genuine mechanical efficiency is already solved for
in §9. Each of the seven revolute joints carries a known reaction $R_j$ and
rotates through a known relative angle $\Delta\theta_j$, so under boundary or
mixed lubrication

$$W_{\text{bearings}} = \sum_j \mu_j \, r_j \oint |R_j| \left|\frac{d(\Delta\theta_j)}{d\theta_1}\right| d\theta_1$$

The piston slides through $|d\lambda|$ against the liner reaction $D$ — the same
quantity the side-load constraint $\gamma$ bounds — plus a roughly constant ring
tension, which at the low loads of an economy engine is the larger term:

$$W_{\text{piston}} = \mu_p \oint \left(|D| + F_{\text{ring}}\right) \left|\frac{d\lambda}{d\theta_1}\right| d\theta_1$$

A spur mesh loses one to two per cent of what it transmits. Then

$$\eta_{\text{mech}} = \frac{W_{\text{ind}} - W_{\text{bearings}} - W_{\text{piston}} - W_{\text{mesh}}}{W_{\text{ind}}}$$

This is what makes the constraint set mean something. Without friction, the
bound $\gamma \le 0.02$ is an assertion about wear that never appears in any
objective, and the bearing-load bound is likewise inert. With it, both feed
directly into the quantity being maximised.

The model is Coulomb with constant coefficients, not hydrodynamic. It gets the
scaling right — loss proportional to load and to sliding distance — and responds
correctly to design changes, which is what an optimizer needs. It will not
predict absolute friction mean effective pressure better than about 30 %, and no
conclusion here rests on the absolute value, only on comparisons at equal
coefficients. The journal coefficient is the most uncertain number in it, so its
influence is swept rather than assumed.

### 11.2 The fuel, from the heat release

The idealised cycle of §3 has exactly one heat input, the constant-volume
combustion step, so the fuel follows from it rather than from an assumed
air-fuel ratio:

$$Q = \frac{V_0 (P_3 - P_2)}{\gamma - 1}$$

Deriving it this way keeps the range model consistent with the thermodynamics
instead of bolting a second, independent combustion model alongside it. The
indicated thermal efficiency $W_{\text{ind}}/Q$ then comes out at 47.7 %, which
is where an idealised Atkinson cycle at this expansion ratio belongs, and the
$p$–$V$ loop area agrees with the torque integral to 2 %.

One consequence should be stated rather than buried: with $k = 1.7$ fixed, the
implied charge is far leaner than any spark-ignition engine would burn. The
specification describes a heavily diluted, lightly loaded cycle. That is not
unreasonable for an economy engine, but the absolute IMEP of 2.7 bar is low and
the model should not be read as describing a conventional engine at full load.

---

## 12. Range, and why it is the right objective

The mechanism is for a Shell Eco-marathon car. That competition scores exactly
one thing: distance on a given quantity of fuel. Once that is written down, the
exchange rates between efficiency, envelope and mass stop being the designer's
taste and become physics.

### 12.1 The engine mass budget

The sizing loop of §9 returns the mass of the seven sized members, about 0.15 kg.
Reporting that as "the mass" and optimising it optimises a tail. The whole
reciprocating assembly is eight contributions, each sized from something the
analysis already knows: the sized members, the piston, the gear pair, the two
shafts, four bearings, the crankcase, the cylinder and head, and the flywheel.

Two of them change the shape of the problem.

**The crankcase converts $H$ and $B$ into kilograms.** A box has to enclose the
mechanism; its walls scale with the envelope area, at a thickness set by
castability and stiffness rather than stress. The main-bearing reactions do not
pass through those walls in bending — they go into local bosses — which is worth
stating because assuming otherwise asks for a 27 mm wall at these loads.

**The flywheel converts torque ripple into kilograms.** A single-cylinder engine
needs enough rotating inertia to carry it through compression:

$$J = \frac{\Delta E}{\delta \, \omega^2}, \qquad
\Delta E = \max_\theta E(\theta) - \min_\theta E(\theta), \qquad
E(\theta) = \int_0^\theta \left(M_r - \bar{M_r}\right) d\theta_1$$

The turning-moment diagram here must be the **gas** torque, not the total. The
inertia part of the full torque history is energy traded with the mechanism's
own moving masses, and that exchange is already accounted for by the inherent
rotating inertia of the cranks and gears; feeding the total in double-counts it
and overstates the flywheel by a factor of five at low speed.

So a linkage whose torque curve is flatter is lighter. No constraint in the
geometric problem expresses that, and it pushes against a long lever arm.

### 12.2 Burn and coast

These cars are not driven at constant speed. The engine runs hard from $v_{lo}$
to $v_{hi}$, is declutched, and the car coasts back down. With
$M v \, dv/dx = F$ the two phases integrate directly:

$$d = \int_{v_{lo}}^{v_{hi}} \frac{M v \, dv}{F(v)}, \qquad
t = \int_{v_{lo}}^{v_{hi}} \frac{M \, dv}{F(v)}$$

with $F = P_w/v - F_{\text{res}}$ during the burn and $F = F_{\text{res}}$ during
the coast.

The energy bookkeeping is worth stating because the naive expectation is wrong.
Over one closed cycle the car starts and ends at the same speed, so the kinetic
energy nets to zero and the propulsive work equals the resistance work over the
**whole** burn-plus-coast distance. Accelerating hard costs nothing in road
load. What the strategy buys is that the engine spends its running time at high
load; what it costs is aerodynamic, since drag goes as $v^2$ and swinging around
a mean burns more than cruising at it.

The minimum-average-speed rule is active at every optimum where the engine has
power to spare, because range always improves as the whole window comes down.
That collapses the two-dimensional strategy search to one dimension: for each
$v_{lo}$, the rule pins $v_{hi}$ exactly, by a bisection on a monotone function.
Solving the active constraint rather than sampling a grid matters for more than
speed — a grid makes the objective a step function of the design variables, and
the gradient-based optimizer downstream would be differentiating quantisation
noise.

### 12.3 The speed trade

Engine speed does not appear in the geometric problem at all. In the coupled
one it is the sharpest single variable, because it pulls two ways:

- **up**, because flywheel inertia for a given speed fluctuation goes as
  $1/\omega^2$, and the flywheel is frequently the heaviest item;
- **down**, because every inertia load goes as $\omega^2$, and §9.5 shows the
  structural mass following as $\omega^6$.

The optimum is interior and is not near either end. The near-singular geometry
of §8 cannot reach it: above 1000 rpm its friction exceeds its indicated work
and the engine will not turn itself over.

---

## 13. The discrete variable in the inter-axle distance

§1 treats the gear radii as continuous functions of $I$. That is a fiction. A
gear has an integer number of teeth cut with a standard-module hob, and for a
2:1 pair on a fixed centre distance

$$r = \frac{mz}{2}, \quad z_1 = 2 z_2 \quad \Longrightarrow \quad
I = \tfrac{3}{2} m z_2$$

so $I$ lives on a lattice. Undercutting sets $z_2 \ge 2/\sin^2\alpha = 17$, hence
$I \ge 25.5 m$.

This is not cosmetic, because $I$ is one of the variables the equalities
$STE = 74$ and $\varepsilon = 16$ are satisfied *with*. Snapping it to the
lattice moves the design off both, and the remaining continuous variables must
repair them. That is the mixed-integer structure of the problem: choose the
integers, repair the continuum.

The repair is substantial, not cosmetic. Moving $I$ onto the nearest buildable
lattice point — 0.18 mm — takes the top-dead-centre gap from 0.003 mm to
0.058 mm, five times its bound, and the expansion stroke 0.087 mm off target.
Since $\partial g/\partial I \approx 0.27$, the gap responds to a perturbation of
the geometry roughly thirty times faster than the perturbation itself.

That number is worth putting beside §14.1. The tolerance study finds
$\sigma_g \approx 0.013$ mm from IT8 machining tolerances, against a 0.01 mm
band; the lattice snap finds a 0.048 mm excursion from a 0.18 mm move. These are
the same fact reached by independent routes: $g$ is hypersensitive to the
geometry that produces it, and no geometric choice — optimised, toleranced or
gear-constrained — can hold it. It has to be an adjustment at assembly.

Enumerating the lattice needs care of its own. Ranking candidates by distance
from the requested $I$ is the obvious choice and the wrong one: the nearest
lattice points are reached with the smallest modules, and a small module needs a
wide face to carry a given tooth load, so ranking by distance offers the
structurally worst candidates first. At the near-singular design of §8 — which
puts 9 kN through the mesh, six times what the backed-off one does — *every*
point on the immediate lattice exceeds the face-width limit, and reaching a
buildable pair means moving the centre distance itself by 13 %. Candidates are
therefore ranked by whether they can carry the load, then by distance.

The practical consequence is that the discrete choice must be **enumerated
outside** the continuous solve, never left floating inside it. Left free, the
module choice makes the objective a step function of $I$: the lightest workable
module changes at a threshold and the range jumps across it. A central
difference straddling that threshold returns a gradient two orders of magnitude
too large, and SLSQP — handed a quadratic subproblem built from constraint
gradients that are quantisation noise — rejects it as *"inequality constraints
incompatible"* and terminates at the starting point without evaluating anything.
Pinning $(m, z)$ also removes $I$ from the design space, leaving a ten-variable
continuous problem per candidate pair.

The same discipline applies to section diameters. Preferred stock sizes are
applied **after** the fixed point of §9.4 converges, never inside it: a step
function inside a contraction mapping turns it into a limit cycle between two
adjacent stock sizes. Rounding is always upward, so it can never turn a section
the sizing certified as safe into one that is not.

---

## 14. Robustness

§8 argues that the published design is not reproducible because it sits where
the mechanism is ill-conditioned. That argument cuts both ways: a design chosen
for its nominal performance in such a region is exactly what a tolerance study
exists to catch, and presenting a deterministic optimum without one would be
negligent.

Tolerances are ISO 286 IT grades. The standard tolerance unit is
$i = 0.45\sqrt[3]{D} + 0.001 D$ µm for a nominal size $D$ in millimetres, and
grade $IT_n$ is a fixed multiple of it — 16i at IT7, 25i at IT8, 40i at IT9. A
machined linkage member is an IT8 part. The two clocking angles are assembly
rather than machining quantities and get a fixed angular band.

Propagation is done two ways. First order,

$$\sigma_g^2 = \nabla g^{\mathsf{T}} \Sigma \, \nabla g$$

costs one extra Jacobian evaluation, because §10 already supplies $\nabla g$
analytically — the exact derivatives built to make the optimizer work pay for
the uncertainty analysis as well. Monte Carlo then checks the linearisation,
which is precisely what should be distrusted near a singularity. On this problem
first order **overestimates** $\sigma$ by up to 80 %, so it is conservative
rather than optimistic; the opposite would make first-order robust design
unusable here.

### 14.1 The gap constraint cannot be manufactured

The result is a finding about the specification, not about any design.

The bound is $g \le 0.01$ mm. The dimensions that produce $g$ are held to
±0.011–0.031 mm each at IT8, and combine to give $g$ a standard deviation of
about 0.013 mm — larger than the constraint band itself. The one-sided process
capability is

$$C_{pk} = \frac{-g_{\text{nominal}}}{3\sigma} = 0.11$$

against an industrial target of 1.33, and roughly two thirds of nominally
conforming builds violate the constraint. Every other constraint is comfortable.

Since $\sigma$ is linear in the tolerance factor, the whole IT ladder can be
scanned from one Jacobian to ask which grade *would* hold it. The answer is
**1.25i**, below the tightest grade in the table. No machining grade fixes it.
The gap has to be taken up by adjustment at assembly — a shim — or the bound
relaxed and $g$ treated as a quantity to minimise rather than to constrain. No
optimizer can help, and one that reports a design satisfying $g \le 0.01$ mm is
reporting a nominal geometry, not a buildable engine.

---

## 15. How coupled is it, and does the formulation matter?

§9 asserts a strong coupling. Gauss–Seidel converges linearly at a rate that
*is* the coupling strength, so the assertion can be a measurement:

$$\rho = \lim_{k\to\infty} \frac{\|r_{k+1}\|}{\|r_k\|}$$

taken over the last few sweeps, where the iteration has reached its asymptotic
rate. At rest $\rho = 0$ exactly — with no inertia there is no path from mass to
load, the quasi-static problem of §4 is recovered, and there is nothing to
iterate. That is the sharpest available check that the measure reflects the
physics rather than the solver. The gain grows with $\omega^2$: 0.13 at 500 rpm,
0.65 at 1000, 0.68 at 1500, at which point each discipline rewrites two thirds
of the other's input on every sweep.

So "is this problem coupled?" has the answer *it depends on the operating
point*, and the dependence is steep.

MDF converges the MDA at every optimizer iteration, so every point evaluated is
physically consistent and the design space stays small. IDF hands the coupling
variables to the optimizer with consistency constraints, so no inner iteration
runs but the design space grows by the dimension of the coupling. The textbook
trade favours IDF when the MDA is expensive and penalises it when the coupling
variables are numerous.

Here the MDA is expensive, so the trade looks open until the coupling is
counted. It is not close. The strong couplings are `diameters` and
`piston_mass`, which are small, and `member_axial` and `member_bending`, which
are not scalars at all: each is the internal load history of every member, at
every crank angle, at every station along it, $7 \times 360 \times 9 = 22\,680$
values. The total is **45 367 coupling scalars against 11 design variables**.

IDF would therefore carry 45 367 extra design variables and 45 367 consistency
equalities in order to optimise eleven real degrees of freedom, and GEMSEO
refuses to build it for exactly that reason. On this problem IDF is not slower;
it is unavailable.

The lesson generalises away from this mechanism. IDF's cost scales with the
*dimension* of the coupling, not with the cost of the MDA, so any discipline
pair exchanging distributed fields -- load histories, pressure distributions,
temperature fields -- sits on the wrong side of the trade however expensive its
MDA. Reducing the coupling to a few scalars first, by a modal or spectral
representation of the histories, is what would make IDF a candidate at all.

One asymmetry is a property of this problem rather than of the formulations.
Intermediate IDF iterates are not physically consistent, so any quantity read
from one — the range, the mass budget, the friction — is meaningless until the
consistency constraints are satisfied. With MDF every iterate is a real engine.

---

## 16. Does the linkage earn its complexity?

The geometric problem cannot pose this question: it has no way to price a
member. With range as the objective it becomes a direct comparison against the
mechanism the EX-link replaces — a conventional slider-crank at the same bore,
clearance volume and compression ratio, whose piston follows
$\lambda = r\cos\theta + \sqrt{l^2 - r^2\sin^2\theta}$ and whose three bodies
give a 9×9 equilibrium system built exactly as §9.1 builds the 18×18 one, and
for the same reason.

The comparison is only meaningful under identical treatment, so the slider-crank
goes through the same sizing, friction, mass-budget and vehicle code; the member
list is a parameter of the sizing routine precisely so that the parity is real
rather than asserted. Only the kinematics, the equilibrium system and the cycle
differ.

Two independent checks: the Otto cycle reproduces $1 - \varepsilon^{1-\gamma}$
to four decimals, and its torque integral matches its $p$–$V$ loop to
$3\times10^{-5}$.

### 16.1 The advantage is firing frequency

The EX-link wins on range, 3338 against 2690 km/L. Decomposed, only about a
fifth of that is extended expansion: indicated efficiency 0.477 against 0.457.
Most of it is *mechanical* efficiency, 0.85 against 0.74 — and the EX-link has
seven journals to the slider-crank's three.

That inversion needs explaining rather than celebrating. As §1 notes, the roles
of the crankshaft and eccentric shaft are exchanged so that the whole four-stroke
cycle completes in **one** turn of the output shaft, while a conventional
four-stroke needs two. Per unit of work the EX-link therefore accumulates half
the journal rotation and half the piston sliding.

That assumption is too load-bearing to leave untested, so the comparison is
re-run with it removed — the same linkage driving a conventional four-stroke gas
exchange, hence twice the friction per cycle. The advantage falls from +24.1 %
to **+2.8 %**.

The honest statement of the result is therefore that the extended-expansion
linkage's benefit in this model comes overwhelmingly from firing frequency, and
that extended expansion alone very nearly fails to pay for the four extra
journals and the gear train it costs.

### 16.2 Conditioning decides the sign of the inertia effect

The second mechanism also settles whether §9.5's finding is a quirk of one
topology.

Speed *reduces* the slider-crank's peak main-bearing load — 4735 N at rest,
2985 N at 4000 rpm — and its structural mass with it. This is correct. The peak
gas force lands near top dead centre, where the reciprocating masses are
decelerating and their inertia force acts against the gas load. It is the
classic **inertia relief**, and the reason high-speed engines do not need
proportionally larger main bearings.

The near-singular EX-link does the opposite, to the point of having no feasible
structure above 1000 rpm. Same physics, opposite sign, and conditioning is what
decides which. Inertia relieves a well-conditioned linkage and amplifies an
ill-conditioned one — which is the general form of the claim that §9.5 could
only state for one mechanism.
