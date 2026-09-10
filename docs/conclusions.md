# 6. Conclusions

## 6.1 What the study establishes

**The objective matters more than the algorithm.** The conventional formulation
prices nothing, its central quantity is not an efficiency, and it cannot see the
parts. Replacing it with range changes the answer qualitatively — not by a few
per cent, but from *the singularity is optimal* to *the singularity is the worst
place to be* (§5.1). Conditioning, not speed, decides the sign of the inertia
effect: a quasi-statically optimised linkage drifts to its transmission-angle
singularity and has no feasible structure above 2000 rev/min, while a
well-conditioned slider-crank shows peak bearing load *falling* with speed.

**Imposing a constraint and checking it are different searches.** The same SLSQP
on the same problem reaches 4.9 % further when the coupled and vehicle
constraints are held during the search rather than verified at the end, and
1.7 % further again when a reliability target is held too — a different and
better design rather than a smaller number (§5.5). Nothing about the algorithm
changed; the problem was posed better.

**A relaxation made for numerical reasons is a promise about tolerance.** The
equalities make the feasible set measure zero, forcing relaxation into bands;
those bands are 1.7 standard deviations wide against the scatter of the parts,
and the reference design has a 66.4 % chance of missing at least one requirement.
Which bounds are responsible is invisible in any nominal quantity: the
top-dead-centre bound is set finer than the model's own resolution, and once it
is widened the band on the expansion stroke governs everything over four orders
of magnitude of failure probability. Widening both costs 0.47 % of range. Most of
the rest is self-inflicted — designs sampled beside the deterministic optimum cut
the probability by 41 % at no cost in range, and the same happens to a
slider-crank held to an active constraint, so the effect belongs to deterministic
optimization under tolerance rather than to this mechanism (§5.2).

**A reliability model fails quietly, in two distinct ways.** It prices what its
uncertain vector contains and reports silence as safety on everything else:
widening the vector from eleven dimensions to seventeen leaves every geometric
constraint identical to the last figure and reveals that the binding constraint
of the whole design — the gear pair's face width, at $\beta = 0.00$ — was never in
the model. And it linearises whatever function it is handed: 150 000 sampled
builds put the study's result at $9.3\times10^{-3}$ against FORM's
$1.3\times10^{-3}$, because the expansion stroke is a maximum over two top dead
centres that the optimizer had driven to within 0.107 µm of each other. Neither
failure showed in any FORM output. Both fixes are free — a different catalogue
pair, and one extra row of a Jacobian already being computed (§5.3).

**The topology is worth 17.6 %, and the dimension count is not the reason.**
Against a conventional engine sized by identical models and optimised over its
own degrees of freedom, the linkage reaches 3395 km/L against 2888, both taking
720° of their crankshaft per cycle so that the comparison needs no correction.
Extended expansion — the feature the topology exists for — is the smallest of the
three contributions, five per cent of indicated efficiency; the larger two are a
lower side load and a lighter flywheel. Handing the conventional engine a third
freedom buys 0.51 %, which extrapolates to 4.6 % over nine rather than 17.6 %: the
dimensions are what let the linkage exploit its topology, not what give it the
advantage (§5.4).

**A schedule of speeds rewards smoothness.** Scored over four speeds as one
engine, both designs lose range and the conventional engine loses more, widening
the advantage to 19.3 %. The reason is structural: a wider speed range makes the
flywheel a larger share of the mass — 84 % of the linkage's and 96 % of the
baseline's — and the flywheel is the item the flat torque curve wins by a factor
of 1.79. The slider-crank is the lighter *mechanism* and the heavier *engine*.

**In a coupled problem the size of an effect is not the size of the thing it acts
on.** Solid round bars looked like the study's weakest assumption and the members
turn out to be 1.3 % of the engine, yet boring them buys 31 % of range at
3200 rev/min, because what a bore removes is inertia in the load path rather than
weight. Given to both engines the comparison moves from +17.6 % to +17.0 %.

**The trigonal link is the feature, not a detail.** Handed the design domain and
no topology, a spring-connected synthesis returns real linkages — discrete,
running at the stiffness floor — and none of them is an extended-expansion
mechanism. The one that reaches a clean four-stroke motion is an Otto engine,
with the two strokes equal to six figures, because it reaches the piston from the
geared shaft alone. That is an identity rather than an accident: a piston driven
from one shaft is periodic in that shaft's angle, so the asymmetry, which lives
entirely in the first harmonic, is *exactly* zero. Extended expansion requires
the piston train to be reached from both shafts, which is what EXlink's
three-cornered link does (§5.6). Handing the search a three-cornered body does
not by itself change the answer: §5.7 offers them, and the start that used two
of them built both on the geared shaft and returned an Otto engine again. The
body has to *bridge* the two shafts, not merely exist.

**A synthesis can be told what it is looking for and still be unable to look.**
The domain was adequate throughout — EXlink is three of its candidates and a gear
pair, reproduced to $2\times10^{-3}$ mm — but the objective scored a degenerate
answer ahead of it, because it charged for where the cycle happened to start.
Profiling that datum out, charging for the two structural conditions the identity
above implies, and drawing starts on the specification's own scale moves every
start from a one-shaft answer to a genuinely bridged one, and one of them to a
discrete mechanism with a two-to-one pair it chose from a catalogue. It is still
not the right motion. What that leaves is a well-posed search problem rather than
a formulation that could never have succeeded (§5.7).

**An under-constrained answer can report anything, and nothing standard catches
it.** Enriching the synthesis domain — three-cornered bodies as single
candidates, shaft angles as coordinates, gear pairs chosen from a catalogue —
produced a start that reported extended expansion with an asymmetry of 70.9 mm
and was not a mechanism at all: it kept no gear, so the geared shaft's angle was
a coordinate nothing resisted, and its motion was chosen by the solver rather
than by the linkage. Strain does not see this, because a loose mechanism strains
nothing. The reading that does is the piston's component in the null space of
the reduced Hessian: 0.71 there, against exactly zero at every genuine answer in
this study. Every enrichment of such a domain has the same exposure, so a
mobility test belongs beside the strain test in any run of this kind (§5.7).

**Decomposition buys structure, not speed.** Bi-level outer approximation halves
the sub-solves against enumeration and lands 0.6 % short, on a bound that is not
valid because the sub-problem is nonconvex. What it buys is a mixed-integer
statement, principled handling of infeasible lattice points, and a stopping
criterion in place of a guessed budget.

## 6.2 Limitations

**Modelling.** Coulomb friction with constant coefficients leaves absolute FMEP
uncertain by some 30 %, though rankings are robust since comparisons are at equal
coefficients. Instantaneous combustion with no heat transfer makes indicated
efficiency optimistic by several points, equally for both mechanisms. Constant
crankshaft speed is already priced by the flywheel sizing, and the pin-jointed
trigonal link is a stiff triangle either way. Two conservatisms are *not*
even-handed and both run in the linkage's favour, so §5.4's 17.6 % is a lower
bound on both accounts: neglecting gas exchange discards a recoverable loss two
and a half times larger for the conventional engine, and the baseline is allowed
to violate two limits the linkage is held to. One runs the other way — the
baseline's flywheel is sized for the worst turning-moment diagram in its class,
though by the same rule as the linkage's. Reliability is compared across
mechanisms of different dimensionality, which is a real difference rather than an
artefact but means §5.2's reliability figures are not like-for-like in the way
§5.4's range figures are.

**Method.** The reliability estimator is first order, and first order is enough
here once applied to the right function — but the branch check has been run at
two designs, and there is no automatic warning when a constraint is a maximum
evaluated near its tie. The widened uncertain vector takes its six added
parameters as independent, which is conservative for the two strengths, and costs
eighteen analyses rather than one Jacobian, so it remains a reporting tool rather
than something the optimizer can call. The mixed-integer bound is not a bound,
outer approximation's guarantee requiring a convex sub-problem this problem
violates comprehensively. Neither reported optimum is converged: both the
headline run and the re-solve stopped at their iteration caps, and the second is
0.13 % above the first, which bounds how much the cap hides without saying where
the formulation ends. The global optimum is not established — uniform multistart
is inapplicable and restoration makes six restarts feasible without answering
whether their solutions agree. The schedule's distance weights are stated rather
than derived, no survey of the track existing.

The topology synthesis of §5.6 searched a domain of twelve candidate members
from six random starts, and §5.7's enriched domain three of six inside a 90-minute cap;
both are small next to the published spring-connected studies. They
establish what extended expansion requires, not that no alternative topology
exists.

**Scope.** The results at $\beta \ge 3$ are stated at a widened specification: the
top-dead-centre gap at 0.1 mm and both equality bands at $\pm 0.15$. §5.2 prices
that widening at 0.47 % of range and shows what it buys, but whether those bounds
are acceptable is a question for the customer rather than for the optimizer. At
the bounds as written the mechanism reaches no reliable design at all. Two
topologies establish a contrast; three would establish a trend.

## 6.3 Further work

Five questions are open, in rough order of value per unit of effort.

*Re-optimise the linkage for a schedule rather than scoring it over one.* §5.4
scores the point optimum over four speeds; an eleven-variable solve against the
schedule would turn +19.3 % from a lower bound into an answer, and the schedule
makes the flywheel — the term the linkage wins on — a larger share of both
engines.

*Settle whether the range optimum is global.* Restoration takes six perturbed
starts from 0 of 6 feasible to 6 of 6, which makes the multistart answerable; the
answer needs a range solve from each, at six times the cost of the headline run.
It is open, and now cheaply askable.

*Automate the kink check.* Nothing in a FORM output signals that it is being
evaluated at a tie between branches of a maximum, and the study found its own
instance by sampling at two designs. A cheap detector — comparing the branch
values against the local scatter — would make the check routine rather than
fortunate. Second derivatives remain worth having for the exact
$\partial\beta/\partial x$ of the steered quantity, but must be taken away from
the tie.

*Widen the uncertainty model further and put it inside the search.* What remains
is a correlated model of the two strengths, second-order treatment of
`saturation` — a threshold on a fixed point, and the least linear constraint in
the set — and an implementation cheap enough for the optimizer to call rather than
for a report to quote.

*Widen the search; the formulation is now sound.* The conditions §5.7 identified
as missing are in it — bridging, mobility, a delivered stroke, and a datum that
is no longer charged for — and they work: a start now returns a discrete,
rigid, fully determined mechanism reaching the piston from both shafts with a
two-to-one pair it chose itself. What it does not yet return is the right motion,
and the reason is the one thing a graph condition cannot express: both shafts
must reach the piston with comparable *authority*, not merely be connected, and
the search commits to a short path in its first rung. The objective ranks EXlink
an order of magnitude ahead of everything four starts found, so what remains is a
search problem and nothing else — a wider multistart, a continuation that defers
the commitment, or restarts built near mechanisms that already balance the two
harmonics. Twelve bars, ten bodies and four starts is small next to the published
spring-connected studies.

*Then reach past the motion.* A kinematic target produces a candidate, and only
the range chain of §3.3 decides whether a candidate beats 3395 km/L. The
top-dead-centre gap that §5.2 shows governs this mechanism is a functional of the
synthesised motion, so it belongs inside the synthesis rather than in a check
afterwards. A functional IDF over the piston motion, sketched in {doc}`Appendix C
<supporting>`, is a larger question again and would need its own study.

## 6.4 Headline numbers

| what the linkage achieves | |
|---|---|
| best strictly feasible design, specification as written | 3338 km/L, 12.2 kg |
| best nominal design, all constraints imposed | 3501 km/L |
| best design that also holds a reliability target, bounds relaxed | **3395 km/L** |
| what the reliability requirement costs | −3 % |

| against a conventional engine | |
|---|---|
| optimised as a conventional engine, both at 720° per cycle | 2888 km/L |
| against the study's result, 3395 km/L | **+17.6 %** |
| against `COUPLED_DESIGN`, 3338 km/L | +15.6 % |
| over a four-point schedule, one engine each | 2733 vs 3260 km/L, **+19.3 %** |
| with tubular members on both, each re-optimised | 2929 vs 3425 km/L, **+17.0 %** |
| indicated efficiency | 0.457 → 0.480 |
| mechanical efficiency | 0.787 → 0.865 |
| engine mass | 16.9 → 12.9 kg |
| the same baseline given a third freedom | 2903 km/L, +0.51 % |

| what the bounds cost | |
|---|---|
| probability the reference design misses a requirement | 66.4 % |
| the same for the best design sampled beside it | 39.3 %, at +0.05 % range |
| gap bound above which the gap stops binding | 0.054 mm; 0.1 mm adopted |
| stroke band the *system* then needs | ±0.15 mm against ±0.05 |
| range given up by widening both | −0.47 % |
| the study's result, corrected constraint set | $P_f = 9.4\times10^{-3}$, $\beta = 2.35$ |
| the same against all thirteen constraints, as built | $5\times10^{-1}$, on the gear pair |
| with the gear pair the exhaustive search preferred | $9.4\times10^{-3}$, at +0.4 km/L |
