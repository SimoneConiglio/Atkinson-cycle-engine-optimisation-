# The discrete variables

The gear pair is not continuous, and pretending otherwise quietly breaks both the geometry and the optimizer.

The 2:1 gear pair appears in the geometric problem only as two primitive radii, both
continuous functions of `I`. That is a fiction. A gear has an integer number of teeth cut with
a standard-module hob, so

> `r = mz/2` and `z₁ = 2z₂` ⟹ **`I = 1.5 · m · z₂`**

`I` lives on a **lattice**, not an interval. Asking for `I = 56.55 mm` gets you 56.40 (m=0.8,
z=47) or 56.25 (m=1.25, z=30), and nothing between. Undercutting sets the floor at `z ≥ 17`,
so `I ≥ 25.5m`.

This matters because `I` is one of the variables the equalities `STE = 74` and `ε = 16` are
satisfied *with*. Snapping it to the lattice breaks both, and the remaining continuous
variables have to repair them — the classic *choose the integers, repair the continuum*
structure.

And the repair is not a formality. Moving `I` onto the nearest buildable lattice point — a
shift of **0.18 mm** — blows the top-dead-centre gap from 0.003 mm to 0.058 mm, five times its
bound, and knocks the expansion stroke 0.087 mm off target. The remaining ten variables then
have a real optimization problem to solve just to get back to feasibility.

Which constraint does the repair struggle with? `g`, every time. That is the same constraint
the [tolerance study](reliability.md) shows cannot be manufactured, and the
two results are the same fact reached independently: `∂g/∂I ≈ 0.27` mm per mm, so `g` responds
to *any* perturbation of the geometry — a machining tolerance, a gear lattice — far faster than
its 0.01 mm band allows. It is not a constraint that geometry can hold.

## Solving it as a mixed-integer problem

Enumerating a handful of lattice points and keeping the best is not a method: it has no
bound, no stopping criterion, and no way to know whether an unvisited point would have won.
`exlink/minlp.py` states the gear choice as the MINLP it is and hands it to the
[`gemseo-bilevel-outer-approximation`](https://pypi.org/project/gemseo-bilevel-outer-approximation/)
plugin, which implements the Duran–Grossmann decomposition:

```
main       gear_choice, a one-hot selection over the lattice   (categorical)
           solved by BILEVEL_MASTER_OUTER_APPROXIMATION
                |                                    ^
                | I, gear_module, gear_teeth         | linearisations and
                v      (catalogue interpolation)     | feasibility
sub        the ten remaining linkage variables       (continuous)
           MDF over the coupled disciplines, solved by SLSQP
```

GEMSEO's `Benders` formulation performs the split itself — categorical variables to the main
problem, everything continuous to a sub-scenario it wraps in an `MDOScenarioAdapterBenders` —
so the main problem optimises the sub-problem's *optimum*. A `CatalogueDesignSpace`
categorical variable drives three catalogue interpolations (`I`, the module and the tooth
count) so that picking a lattice point sets all three; at unit SIMP penalty the interpolation
is exactly `I = Σⱼ yⱼ Iⱼ`, linear and analytically differentiable, which is what the
outer-approximation master linearises.

Two changes were needed to make the package fit the formulation, and both are the formulation
being right rather than convenient:

- **The gear pair became `RangeDiscipline` *inputs*.** With the module fixed at construction
  time there is nothing for a master to choose.
- **The range margins are published in both sign conventions.** `positive=True` makes GEMSEO
  rename a constraint to `-runs_margin`, which is not an output of anything, and the scenario
  adapter addresses constraints by output name. So `runs_violation` and `gear_violation` are
  emitted alongside the margins, and the bi-level formulation attaches those.

**Infeasible sub-problems need no special machinery**, which is the reason OA fits here.
Pinning `I` throws the design off the equalities `I` was one of the variables used to satisfy,
so several lattice points have no feasible continuous solution at all. Attaching the
constraints with `main_level=True` puts an `is_feasible` condition on the main problem, and
such a point is excluded on evidence — carrying exactly the information enumeration discarded.

Outer approximation's finite convergence to the *global* optimum needs the sub-problem convex
in the continuous variables for each fixed choice, which this problem violates comprehensively.
So the master's bound is a bound under an assumption that does not hold, and `minlp.exhaustive`
solves every candidate separately so the decomposition's answer can be checked against the
true best over the lattice.

## Convexification

Outer approximation's premature convergence on a nonconvex problem is exactly what the
plugin's convexification options exist for, and they are enabled by default: `posa` amplifies
the cut slopes so a linearisation is less likely to exclude a lattice point that is actually
better, and `adapt` corrects them by a secant method over the visited history so the cuts
behave like valid supports where the value function is not convex.

On this problem, measured over four candidates, they change nothing:

| setting | chosen | range | sub-solves |
|---|---|---|---|
| raw cuts (`posa=1`, `adapt=off`) | m=0.8, z=48 | 3366 km/L | 2 |
| `posa=2` only | m=0.8, z=48 | 3366 km/L | 2 |
| `adapt` only | m=0.8, z=48 | 3366 km/L | 2 |
| `posa=2` + `adapt` | m=0.8, z=48 | 3366 km/L | 2 |

The reason is structural rather than a failure of the options. The master terminates after two
sub-solves, and the adaptive secant correction needs more history than that before it can
adjust anything; and because the two equality constraints reach the master only through its
feasibility condition rather than as linearisations, there is little cut information for
`posa` to steepen. The options are left on because they are the right default — they can only
make the master more conservative — but on this problem they are a measured no-op, and saying
so is more useful than implying they fixed the shortfall.

## Why the equalities stay equalities

The plugin's source notes that the master does not linearise equality constraints and
suggests relaxing them into inequality bands first. Doing that here exhausts memory: four
extra sub-problem constraints, each carrying a post-optimal sensitivity through an MDF
sub-problem whose coupling variables are the **45 367** load-history entries counted in
[How coupled is it, really?](coupling.md). The process is OOM-killed past 15 GB.
With equalities the formulation adds a single `is_feasible` condition instead and the run
fits in **0.51 GB**.

That is the IDF result from earlier arriving in a second place: it is the *dimension* of the
coupling, not the cost of the MDA, that decides what a decomposition can afford. The cost of
the workaround is stated above — those two constraints inform the master only through
feasibility.

That check is worth reporting rather than skipping. Over four candidates at 1000 rpm, with a
25-iteration SLSQP budget per sub-problem:

| | chosen pair | range | sub-solves | seconds |
|---|---|---|---|---|
| outer approximation | m=0.8, z=48 | 3366 km/L | **2** | 575 |
| exhaustive | m=1.0, z=39 | 3385 km/L | 4 | 1056 |

The decomposition costs half the sub-solves and lands **0.6 % short** of the best point on the
lattice. That is exactly what nonconvexity buys you: the master stopped on a bound that is not
valid here, so it terminated before reaching the best candidate. The honest summary is that on
this problem the formulation's value is structural — a real mixed-integer statement, principled
handling of infeasible sub-problems, and a stopping criterion instead of a guessed budget —
rather than a better answer than enumeration. With a lattice too large to enumerate, that
structure is the only thing on offer.

## Manufacturability

`manufacturing.py` holds the R20 preferred bar diameters, the ISO 54 module series, and
minimum castable and machinable wall thicknesses. Rounding is applied *after* the fixed point
converges, never inside it: a step function inside a contraction turns it into a limit cycle
between two stock sizes. Rounding is always **up**, so it can never turn a certified section
unsafe, and `stock_premium` reports what buildability cost.

---

Next: [Whether the answer survives manufacturing](reliability.md)
