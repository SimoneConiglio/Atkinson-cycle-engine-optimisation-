# Sizing, dynamics and the coupling

Restoring inertia to the load path closes a loop between the section sizes and the loads they carry. This is the multidisciplinary core of the study.

A quasi-static study cannot do this, and the reason is circular: the inertia loads need the
part masses, and the masses need the sections, which need the loads. Restoring inertia
therefore closes a loop:

```
    section diameters ──▶ member masses ──▶ inertia forces ──▶ internal loads
            ▲                                                        │
            └──────────── sizing against yield, ─────────────────────┘
                          fatigue and buckling
```

Neither half can go first. Under GEMSEO that is the **MDF** formulation:
`DynamicsDiscipline` and `StructureDiscipline` are strongly coupled — `diameters` one way,
the internal load histories the other — and an MDA converges them before the optimizer sees
a number.

```bash
exlink size --rpm 1000                    # size every part, inertia included
exlink size --rpm 3000 --plot fig.png     # ...and watch the loop run away
```

```python
from exlink import solve_for_design
sized = solve_for_design(PUBLISHED_DESIGN, speed_rpm=1000)
print(sized.total_mass_kg, sized.diameters, sized.feasible)

from exlink.scenarios import minimise_mass
best = minimise_mass(speed_rpm=1000, min_efficiency=0.25)
```

**Why sequential elimination stops working.** Without inertia every rod is a two-force
member — the forces at its two ends are equal, opposite and collinear with the rod — which is
what lets the loads be eliminated one body at a time, from the piston down to the crankshaft.
Give a rod mass and its end forces are neither collinear nor equal, so nothing can be solved
before its neighbours. Gruebler gives the mechanism one degree of freedom, so the load
problem is statically determinate — **18 unknowns against 6 bodies × 3 equations** — and
`exlink.dynamics` assembles and solves that 18×18 system at every crank angle. Its
determinant is the same quantity condition (4a) protects: at a critical configuration the
matrix goes singular and the forces blow up.

Each member is then sized against **static yield**, **Goodman fatigue** (Marin-corrected
endurance limit, evaluated per extreme fibre) and **Euler buckling**, by bisection — needed
because the fatigue size factor `k_b` itself depends on the diameter being solved for.

## What it buys

![bearing load and torque against crank angle, at three speeds](figures/bearing_loads.png)

**Efficiency does not care about speed; the parts care enormously.** At constant speed the
mechanism returns to its starting state each revolution, so inertia does no net work and the
*mean* torque is provably unchanged — verified to 1e-6 across 0–3000 rpm. Only the peaks
move, and with the gas load switched off every reaction scales as exactly `Ω²`. Since mass
goes as the cube of the acceleration level, it goes as the **sixth power of engine speed**:

| speed | moving mass | peak bearing load |
|---|---|---|
| 0 rpm | 0.25 kg | 7.7 kN |
| 1000 rpm | 1.03 kg | 12.5 kN |
| 1500 rpm | 8.43 kg | 245 kN |
| 2000 rpm | *loop runs away — no section is thick enough* | |

![sized sections, mass and binding failure mode](figures/sizing.png)

Fatigue binds on six of the seven members; the long slender piston rod goes to Euler
buckling instead.

**The quasi-static optimum turns out to be the wrong answer.** Maximising efficiency without
inertia drives the design to `W = 0.981`, a hair from the transmission-angle singularity,
because that is where the quasi-static lever arm is longest. But that same proximity is what amplifies accelerations — joint `A` sees 75× the
crank pin's. Backing off costs nothing and saves almost everything:

| swing rod | `W` | `η` | `H` [mm] | mass [kg] | peak bearing [N] |
|---|---|---|---|---|---|
| ×1.00 | 0.9811 | 28.20 % | 238.5 | 1.039 | 12 629 |
| ×0.94 | 0.9670 | 27.79 % | 227.8 | 0.610 | 6 541 |
| ×0.88 | 0.9560 | 27.92 % | 218.0 | 0.498 | 6 647 |
| ×0.82 | 0.9488 | 28.56 % | 213.2 | 0.450 | 6 027 |

Half the bearing load, *better* efficiency, a smaller envelope — and less than half the mass,
at 1000 rpm. At 1500 rpm the same comparison spans 8.43 kg against 2.03 kg. The near-singular
design was an artefact of leaving inertia out.

Both tables above are printed by `python examples/05_sizing_and_dynamics.py` (about 30 s).

The coupled problem therefore gains an objective a quasi-static formulation cannot express —
structural mass, since nothing in it determines a cross-section — and three constraints that
only exist once the loads are dynamic: `saturation_margin` (the loop ran away),
`slenderness_margin` (a "rod" thicker than a third of its length is not a beam), and
`bearing_margin`.

Full derivation, including the one idealisation used for the trigonal link and why the loop
converges at all, is in [the theory notes](theory.md) §9.

---

## How strongly coupled is it?

"Strongly coupled" is an adjective. Gauss–Seidel converges linearly at a rate that *is* the
coupling strength, and the residual history already records it, so the claim can be a number:

```
      rpm       rho   sweeps   per decade   verdict
        0    0.0000        2          inf   weak
      250    0.0359        6          0.7   weak
      500    0.1307        9          1.1   moderate
      750    0.2558       12          1.7   moderate
     1000    0.6513       28          5.4   strong
     1250    0.6716       38          5.8   strong
     1500    0.6819       42          6.0   strong
```

At rest `ρ = 0` **exactly**: with no inertia there is no path from mass to load, so the
quasi-static problem is recovered and there is nothing to iterate. That is the sharpest check
that the measure reflects the physics rather than the solver. The gain grows with `ω²`, and by
1500 rpm each discipline rewrites two thirds of the other's input on every sweep.

So the answer to "does this problem need an MDA?" is *it depends on the operating point*, and
the dependence is steep.

---

Next: [The derivatives that make it tractable](gradients.md)
