# Verification and validation

Which claims are checked against independent results, and which are not.

The model is checked against physical invariants rather than against another
implementation, so the checks stand on their own.

**Rigid-body kinematics.** Every link keeps its length to 1e-9 mm over the revolution, and
the piston pin and crown stay on the cylinder axis to the same tolerance
(`tests/test_kinematics.py`).

**The quasi-static force chain** is pinned by the principle of virtual work: in a massless,
frictionless mechanism the instantaneous power in equals the power out, so
`M_r(θ₁) = −P dλ/dθ₁` at *every* crank angle. The chain reproduces that to machine precision
(`tests/test_loads.py`), which constrains the whole chain — piston, trigonal link, both
shafts and the gear pair — not just its endpoints. A second, independent route agrees: the
mean torque equals the indicated p–V loop area over 2π.

That check earns its keep. It caught a **sign error in the trigonal-link moment inversion**
while the model was being built — the version that disagreed with virtual work by a factor
of about −4. Expanding `DA ∧ F_A + DE ∧ F_E = 0` gives
`−c A sin(θ_a − θ_T) − C (DE ∧ û_e)_z = 0`, so the swing-rod load carries a leading minus
sign; see `exlink/loads.py`.

**Dynamics.** At zero speed the 18×18 simultaneous solve reproduces the sequential
quasi-static elimination exactly — torque, joint forces and gear load to 1e-11. Mean torque
is provably independent of engine speed (inertia does no net work over a closed cycle) and
holds to 1e-6 from 0 to 3000 rpm. With the gas load removed every reaction scales as exactly
`Ω²`.

**Derivatives.** Every analytic derivative is compared against a *converged* central
difference, and independently by GEMSEO's own `check_jacobian`. Details in
[Gradients](#gradients).

**The gear pair** transmits no net power: with `ω₂ = −2ω₁` and `r₁ = 2r₂` the two gear
torques cancel exactly, checked explicitly — a pair that generated power would silently
inflate the efficiency.

---

## Does the model land in the right place?

Nothing here is calibrated to a target, so the agreements are checks rather than fits:

| quantity | model | expected |
|---|---|---|
| indicated thermal efficiency | 47.7 % | 40–55 % for an idealised Atkinson cycle |
| p–V loop area vs torque integral | agree to 2 % | identical, by virtual work |
| IMEP | 2.7 bar | low, and correctly so for `k = 1.7` |
| FMEP | 1.2 bar | 0.5–1.5 bar for a small single |
| range | 2100–3400 km/L | 2000–3500 for Prototype-class gasoline |

## Is the linkage worth it?

The geometric formulation could not pose this question: it has no way to price a member. With
range as the objective it is a straight comparison against the mechanism the EX-link
replaces — a conventional slider-crank at the same bore, same clearance volume, same
compression ratio.

The comparison is only meaningful if both are treated identically, so the slider-crank goes
through the *same code* wherever the code is not topology-specific: `size_from_arrays` takes
its member list as a parameter precisely so the parity is real rather than asserted. Same
material, same yield/fatigue/buckling checks, same friction coefficients, same crankcase,
bearing, shaft and flywheel models, same vehicle. What differs is only what must — the
kinematics, the equilibrium system, and the cycle.

Independently validated: the Otto cycle reproduces `1 − ε^(1−γ)` to four decimals, and the
torque integral matches the p–V loop to 3 × 10⁻⁵.

| | slider-crank (Otto) | EX-link (Atkinson) |
|---|---|---|
| members / journals | 2 / 3 | 7 / 7 |
| indicated efficiency | 0.457 | 0.477 |
| mechanical efficiency | 0.740 | 0.853 |
| brake efficiency | 0.338 | 0.407 |
| engine mass | 19.3 kg @ 2000 rpm | 12.2 kg @ 1000 rpm |
| **range** | **2690 km/L** | **3338 km/L** |

The EX-link wins by 24 %. But decomposing it, only about a fifth of that is extended
expansion — indicated efficiency 0.477 against 0.457. Most of it is *mechanical* efficiency,
0.85 against 0.74, and the EX-link has **seven** journals to the slider-crank's three.

That inversion needed explaining rather than celebrating. The cause is firing frequency: in
this model the EX-link completes four strokes in one crankshaft revolution (see [The design
problem](problem.md)), so per unit of work it accumulates half the journal rotation
and half the piston sliding of a four-stroke. That assumption was too load-bearing to leave
untested, so `firing_frequency_sensitivity` re-runs the comparison with it removed:

| | range | advantage |
|---|---|---|
| slider-crank | 2690 km/L | — |
| EX-link, as modelled | 3338 km/L | **+24.1 %** |
| EX-link, if it were a four-stroke | 2765 km/L | **+2.8 %** |

**The advantage is firing frequency, not extended expansion.** Reported the other way round it
would have been wrong.

## The generalisation

The second mechanism also settles whether the singularity finding is a quirk of one topology.
It is not, and the contrast is sharper than expected. Speed *reduces* the slider-crank's peak
main-bearing load:

| speed | peak main-bearing load | linkage mass |
|---|---|---|
| 0 rpm | 4735 N | 48.4 g |
| 2000 rpm | 4296 N | 45.1 g |
| 4000 rpm | 2985 N | 44.1 g |

That is correct, not a bug. The peak gas force lands near top dead centre, where the
reciprocating masses are decelerating and their inertia pulls the other way — the classic
**inertia relief** of the gas load, and the reason high-speed engines do not need
proportionally bigger main bearings.

The near-singular EX-link does the exact opposite: it has no feasible structure at all above
1000 rpm. Same physics, opposite sign, and **conditioning decides which**. Inertia relieves a
well-conditioned linkage and amplifies an ill-conditioned one. That is the central claim of
this study, now stated on two mechanisms instead of one.

---

Next: [What the optimizations produced](results.md)
