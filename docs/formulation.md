# Formulation

How the problem is posed for GEMSEO: what is optimised, what the disciplines are, and how they are coupled.

## Strategies available

The formulation is exposed through GEMSEO, so several strategies apply to the same problem:

| strategy | entry point |
|---|---|
| external penalty `F(X) = −η + r⁻²(c_eqᵀc_eq + ⟨c⟩ᵀ⟨c⟩)` | `PenalisedExlinkDiscipline` |
| gradient-based local search | `SLSQP` (`maximise_efficiency`, `minimise_mass`) |
| derivative-free local search | `NELDER-MEAD`, `NLOPT_COBYLA` |
| global search | `DIFFERENTIAL_EVOLUTION` (`maximise_efficiency`) |
| multi-objective front | `pareto_front` / `local_pareto` (`PYMOO_NSGA2`) |
| moving limits on `H` and `B` | `sweep_moving_limits` |
| ε-constraint on efficiency | `sweep_efficiency_floor` |
| constraint-accurate polish | `Augmented_Lagrangian_order_0` (`refine`) |

Which of these actually work on this problem is not a matter of taste, and the rest of this
section reports what was measured.

## Range: the objective the problem was always about

Everything above stops at the engine. But the engine is for a Shell Eco-marathon car, and
that competition scores exactly one thing — **distance on a given quantity of fuel** — so the
three-objective formulation was never the real problem. It could only ever produce a Pareto
front, because nothing in it says what a millimetre of height is worth in points of
efficiency.

Range prices all of them, and the prices are not adjustable:

```
eta_mech  --.
             +--> brake efficiency --.
cycle work -'                         +--> fuel per metre --> RANGE
                                     /
H, B  --> crankcase --.             /
torque ripple --> flywheel --> engine mass --> rolling resistance
```

Two pieces had to be built before that chain would close.

**A real mechanical efficiency.** The package's `η` is *not* an efficiency: with no friction
in the model, the virtual-work identity makes the torque's work exactly equal the gas force's
work at every crank angle. It is a kinematic quality measure — how much *mean* torque a
piston motion converts into — and nothing is lost in it. `friction.py` supplies the losses
that are really there, from quantities the dynamic solve already produces: seven journals
under known reactions turning through known relative angles, the ring pack and skirt sliding
against the liner reaction the side-load constraint already bounds, and the gear mesh. This
is what makes the constraint set *mean* something: a design that leans on the liner now burns
its fuel on the liner.

**A mass budget worth optimising.** The coupled sizing loop returns the mass of the seven
sized members: about 0.15 kg. No engine weighs 0.15 kg, and optimising that number optimises
a tail while a twelve-kilogram dog goes unmodelled.

## The mass budget

Eight contributions, each sized from something the analysis already knows. At the coupled
reference design, 1000 rpm:

| item | mass | share | set by |
|---|---|---|---|
| flywheel | 9.54 kg | 78 % | the cyclic torque fluctuation the linkage produces |
| crankcase | 0.90 kg | 7 % | **the envelope `H × B`**, at a stiffness-driven wall |
| shafts | 0.65 kg | 5 % | combined bending and torsion at the journals |
| bearings | 0.48 kg | 4 % | journal diameter, via catalogue proportions |
| cylinder + head | 0.26 kg | 2 % | bore and peak pressure |
| gears | 0.15 kg | 1 % | peak tooth load, and the chosen module |
| linkage | 0.15 kg | 1 % | the sizing fixed point |
| piston | 0.04 kg | < 1 % | peak gas pressure |
| **total** | **12.17 kg** | | |

Two of these change the *shape* of the problem rather than its scale.

**The crankcase makes `H` and `B` physical.** A box has to enclose the mechanism and its walls
scale with the envelope. The two envelope objectives convert to kilograms at a rate the
physics fixes, and kilograms convert to range through rolling resistance. That is what
collapses three objectives into one.

**The flywheel makes torque smoothness physical.** A single-cylinder engine needs enough
rotating inertia to carry it through compression, and the requirement follows from the
*fluctuation* of the turning-moment diagram. A linkage with a flatter torque curve is
lighter — a design driver no geometric constraint expresses, pushing against a long lever arm.
It is also, at these speeds, three quarters of the engine.

Sizing it correctly took two corrections. Feeding the *total* torque into the flywheel
calculation double-counts energy the mechanism's own masses already store, overstating it
fivefold at low speed; it uses the gas turning-moment diagram instead. And carrying the main
bearing reaction through a plate spanning the whole crankcase asked for a 27 mm wall — the
load goes into a local boss, and the wall is set by castability and stiffness.

## How the car is driven

Not at constant speed. Every serious team drives **burn and coast**: run the engine hard from
`v_lo` to `v_hi`, declutch, coast back down, repeat. This is where the accelerations and
decelerations enter, and the naive expectation about them is wrong. Over one closed
burn-and-coast cycle the car starts and ends at the same speed, so the kinetic energy nets to
zero and

> `W_burn = ∫ F_res(v) dx` over the **whole** burn *and* coast distance.

Accelerating hard costs **nothing** in resistance work. What burn-and-coast buys is that the
engine spends its running time at high load; what it costs is aerodynamic, since drag goes as
`v²` and swinging around a mean burns more than cruising at it. Both are in the model, and
there is a test pinning the energy balance.

The minimum-average-speed rule is active at every optimum where the engine has power to
spare, which collapses the two-dimensional strategy search to one dimension: for each `v_lo`,
the rule pins `v_hi` exactly. That is not only faster — a grid search would make the objective
a step function of the design variables, and the gradient-based optimizer downstream would be
differentiating quantisation noise.

## MDF or IDF?

`build_coupled_scenario` takes a `formulation_name`, so both can be run on the identical
problem. The comparison has a decisive answer, and it is structural rather than a matter of
timing:

| | objective | evaluations | seconds | feasible |
|---|---|---|---|---|
| MDF | 0.1444 kg | 25 | 41 | yes |
| IDF | — | 0 | — | **cannot be posed** |

IDF hands the coupling variables to the optimizer as design variables with consistency
constraints. Counting them settles it:

```
  design variables            11
  coupling variables       45368
    of which load histories  45360   (22680 each)
```

`member_axial` and `member_bending` are not scalars. They are the internal load history of
every member, at every crank angle, at every station along it — 7 × 360 × 9 apiece. IDF would
carry **45 367 variables and 45 367 equality constraints to optimise eleven real degrees of
freedom.**

So IDF is not slower here; it is unavailable. The general lesson separates cleanly from this
mechanism: IDF's cost scales with the *dimension* of the coupling, so a discipline pair that
exchanges distributed fields — load histories, pressure distributions, temperature fields —
rather than a handful of scalars will sit on the wrong side of that trade however expensive
its MDA.

---

Next: [The coupling at the centre of it](coupling.md)
