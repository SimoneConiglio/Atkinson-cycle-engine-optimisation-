# 4. Implementation framework

§3 derived the methodology without reference to code. This section maps each
element onto the implementation and says where the approximations live. Every
module is documented in {doc}`api`.

## 4.1 Structure

The package separates the physics, which is plain NumPy and independently
testable, from the optimization, which is GEMSEO.

```
  physics                       optimization
  ────────────────────────      ───────────────────────────────
  design      kinematics        disciplines   GEMSEO wrappers
  constants   cycle             scenarios     design space, constraints,
  materials   loads                           workflows
  metrics     dynamics          formulations  coupling strength, MDF/IDF
  model       sizing            minlp         bi-level outer approximation
  derivatives coupled           robustness    tolerance and reliability
  friction    gears
  mass_budget manufacturing     jacobian              exact d/dX of the chain
  vehicle     performance       dynamics_jacobian     exact d/dX through the MDA
  slidercrank                   plots, animation, cli
```

## 4.2 Where each method lives

| §3 | method | module | notes |
|---|---|---|---|
| 3.3 | sizing/dynamics fixed point | {mod}`exlink.coupled` | reference Gauss–Seidel implementation |
| 3.3 | the $18\times18$ equilibrium solve | {mod}`exlink.dynamics` | d'Alembert, statically determinate |
| 3.3 | yield / fatigue / buckling | {mod}`exlink.sizing` | Goodman with Marin factors |
| 3.5 | forward-mode + envelope theorem | {mod}`exlink.jacobian` | exact, verified against converged differences |
| 3.5 | derivatives through the MDA | {mod}`exlink.dynamics_jacobian` | implicit function theorem on the bisection |
| 3.6 | MDF, and the coupling measurement | {mod}`exlink.formulations` | $\rho$ from the MDA residual history |
| 3.7 | bi-level outer approximation | {mod}`exlink.minlp` | `gemseo-bilevel-outer-approximation` |
| 3.7 | the gear lattice | {mod}`exlink.gears` | ISO 54 modules, undercut limit |
| 3.8 | tolerance and reliability | {mod}`exlink.robustness` | ISO 286, FORM, correlated orthant |
| 3.9 | restarts on the manifold | {mod}`exlink.scenarios` | `multistart`, `project_onto_equalities` |
| 3.10, 7.4 | prescribed motion, and range under every constraint | {mod}`exlink.synthesis` | the target as a fallback objective |

The section numbers in the first column were stale against §3's headings and
are corrected here; the mapping is worth checking rather than trusting, because
a limitation and its resolution live in different sections.

The objective chain of §3.1 is assembled in {mod}`exlink.performance`, which
composes {mod}`exlink.friction`, {mod}`exlink.mass_budget` and
{mod}`exlink.vehicle` and returns every intermediate so a result can be
interrogated rather than believed.

## 4.3 What is exact and what is not

The distinction matters because §3.5 argues that finite differences are wrong on
part of this problem. They are still used where they are safe, and the boundary
is deliberate:

| quantity | derivative | why |
|---|---|---|
| stroke, compression ratio, $W$, $mra$, $g$, $\gamma$ | **analytic** | extremum-based; differences are wrong (§3.5) |
| the $18\times18$ solve, the sizing bisection | **analytic** | the MDA would otherwise cost 11 fixed points per gradient |
| clearance $d$ | difference | a minimum over both crank angle and three edges; far from active |
| $\eta$, $H$, $B$ | difference | smooth, none tight |
| the range chain | difference | one load solve, no fixed point; ~0.3 s, so 18 columns is affordable |
| the reliability margin | difference | exact would need $\nabla^2 g$; see §7.2 |

## 4.4 Verification

Each discipline is checked against a result computed independently of it.

| check | agreement |
|---|---|
| force chain against virtual work, $M_r = -P\, d\lambda/d\theta_1$ | machine precision |
| torque integral against the p–V loop area | 2 % |
| mean gas torque against mean total torque (inertia does no net work) | $10^{-6}$ |
| rigid-link invariants over the revolution | machine precision |
| joint reactions scaling as $\Omega^2$ with gas load off | exact |
| slider-crank indicated efficiency against $1-\varepsilon^{1-\gamma}$ | 4 decimals |
| analytic Jacobians against converged central differences | $\sim 10^{-6}$ relative |
| GEMSEO `check_jacobian` on the disciplines | passes |
| FORM system probability against 4000-sample Monte Carlo | 0.645 vs 0.653 |

## 4.5 Running it

Install, command line, module map and a table mapping each result in §6 to the
example that produces it are in {doc}`implementation`.

---

Next: [5. The use case](use_case.md)
