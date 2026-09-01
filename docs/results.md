# Results

What the optimizations actually produced, and what they did not.

## The reference designs

Each ships in {mod}`exlink.reference` and is reproducible from the examples. `eta` is the
lever-arm quality measure of {doc}`problem`, not a thermodynamic efficiency; feasibility is
against the geometric constraint set at 720 crank samples.

| design | `eta` | `H` mm | `B` mm | `W` | `g` mm | feasible |
|---|---|---|---|---|---|---|
| `PUBLISHED_DESIGN` | 35.62 % | 283 | 157 | 0.9892 | 8.5236 | no |
| `REFINED_DESIGN` | 27.80 % | 239 | 152 | 0.9811 | 0.0060 | yes |
| `GRADIENT_DESIGN` | 30.91 % | 320 | 159 | 0.9850 | 0.0095 | yes |
| `COUPLED_DESIGN` | 25.00 % | 198 | 131 | 0.9372 | 0.0070 | yes |
| `RANGE_DESIGN` | 25.46 % | 231 | 131 | 0.9319 | 0.0012 | no |

`PUBLISHED_DESIGN` is a historical baseline and does not reproduce; the reason is instructive
rather than alarming and is set out in {doc}`problem`. `GRADIENT_DESIGN` maximises `eta`
alone, which is unbounded in mechanism size — hence `H = 320` — and {doc}`search` shows it was
a local optimum besides. `COUPLED_DESIGN` is the one to compare against: it gives up five
points of `eta` to move off the transmission-angle singularity, and gets a lighter, faster,
longer-ranged engine for it.

## Vehicle-level performance

At 1000 rpm, through the whole chain:

| | `COUPLED_DESIGN` | slider-crank |
|---|---|---|
| indicated thermal efficiency | 47.7 % | 45.7 % |
| mechanical efficiency | 0.853 | 0.740 |
| brake thermal efficiency | 0.407 | 0.338 |
| engine mass | 12.17 kg | 19.3 kg |
| **range** | **3338 km/L** | **2690 km/L** |

The decomposition of that 24 % advantage — and why most of it is not extended expansion — is
in {doc}`validation`.

## The range optimization

SLSQP on `neg_range`, gear pair pinned, 1000 rpm, started from the coupled reference:

| | range | engine mass | `g` | strictly feasible |
|---|---|---|---|---|
| start (`COUPLED_DESIGN`) | 3338 km/L | 12.17 kg | 0.0067 mm | **yes** |
| best found (`RANGE_DESIGN`) | 3388 km/L | 12.47 kg | 0.0009 mm | no — see below |

The 1.5 % gain is modest, and the reason it is not simply banked is worth stating rather than
smoothing over.

`RANGE_DESIGN` satisfies every inequality, including the gap, at `g = 0.0009 mm`. It misses
the two *relaxed equalities* by 1.5 × 10⁻⁴ mm and 6.1 × 10⁻⁵ — SLSQP stopping within its own
convergence tolerance of the constraint it was handed. For scale, the tolerance study puts
the machining standard deviation of `STE` at 0.020 mm, **130 times larger**; no real part
would tell the two apart.

The obvious fix is to project it back onto the equality manifold, which
`project_onto_equalities` does exactly, by the minimum-norm Newton step from the analytic
Jacobians. That step is a few hundredths of a millimetre — and it moves `g` from 0.0009 to
0.0201 mm, twice its bound.

So the same wall appears from a fourth direction:

| perturbation | effect on `g` (bound: 0.01 mm) |
|---|---|
| IT8 machining tolerance on the members | `σ = 0.013 mm` |
| snapping `I` 0.18 mm onto the gear lattice | `0.003 → 0.058 mm` |
| minimum-norm equality projection | `0.0009 → 0.0201 mm` |
| tightest ISO grade that would hold it | 1.25i — off the ladder |

The honest reading is that **the specification is over-constrained**. The equality manifold
and the region `g ≤ 0.01` intersect in a sliver too thin for machining, gear selection or a
converged optimizer to land inside reliably. Treat `g` as an assembly adjustment — a shim on
the piston-rod length — or as a quantity to minimise, and `RANGE_DESIGN` is the answer at
3388 km/L. Under the specification as written, the best strictly feasible design is the one
we started from, and the 1.5 % is the price of a constraint that cannot be held.

That is not a result the optimizer could have delivered. It came out of the tolerance study,
and it is the single most useful thing in this repository for anyone who would actually build
the engine.

---

Next: [The full derivation](theory.md)
