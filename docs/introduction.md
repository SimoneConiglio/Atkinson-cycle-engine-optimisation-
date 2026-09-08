# 1. Introduction

An engine built for a fuel-economy competition is scored on one quantity: the
distance it carries a car on a fixed quantity of fuel. Designing its mechanism
is therefore a problem with a single scalar objective, and the difficulty is
that almost nothing about the mechanism affects that objective directly.

The linkage studied here is an *extended-expansion* mechanism. It expands the
burnt gas through a larger volume ratio than it compressed the fresh charge,
recovering work that an Otto cycle sends out of the exhaust. Realising that
mechanically requires a piston that reaches top dead centre twice per cycle with
two *different* bottom dead centres — the short one setting the compression
stroke, the long one the expansion stroke. The design freedom is eleven linkage
dimensions and a gear pair.

## 1.1 Four routes to the objective, and they conflict

The design choices reach the fuel consumption of a car by four routes.
*Thermodynamically*, the expansion-to-compression ratio sets how much of the
heat release becomes indicated work; more expansion is better and needs a longer
stroke. *Mechanically*, every joint reaction and every unit of piston side load
is friction, which is indicated work that never reaches the crankshaft, and a
linkage with a long lever arm produces high joint loads. *Inertially*, the parts
have mass, and at speed their accelerations load the same joints — while the
masses depend on the sections, the sections on the loads, and the loads on the
masses, a closed loop. And at *vehicle level*, engine mass is carried for the
whole distance and rolling resistance is proportional to it, so a heavier engine
that is thermodynamically better may still lose.

A formulation that omits any of these gets the answer wrong in a specific and
predictable way. The formulation conventionally applied to this class of
mechanism — maximise a kinematic quality measure subject to envelope bounds,
with the performance requirements as equalities — omits three of the four, and
§2.2 sets out why.

## 1.2 Two structural properties decide what may be applied

Beyond the choice of objective, two properties of the problem determine which
methods are admissible at all.

Two of the seven requirements are **equalities**, which makes the feasible set a
measure-zero manifold in the design space: no method that proceeds by sampling
can find a feasible point in it. This is a property of the specification rather
than of the mechanism, and the repair belongs with the specification, so §3.5
relaxes the equalities into **tolerance bands** as part of stating the problem
rather than as a numerical expedient later. Those bands turn out to be comparable
in width with the manufacturing scatter of the parts, so a deterministic answer
inside the band is not the same thing as an engine that meets the requirement.

Neither property is a numerical inconvenience to be worked around. The first
dictates the class of optimizer; the second changes what "feasible" means and
obliges a probabilistic treatment. Reliability-based optimization arises here as
a *consequence* of the relaxation rather than as an addition to the problem.

## 1.3 Contribution

The study contributes, on a problem small enough to be verified throughout:

1. a formulation in which the objective is the application-level figure of merit
   and the competing geometric quantities are priced by the physical chain
   connecting them to it (§3);
2. a demonstration that each structural property of the resulting problem — a
   feasible set thin enough to exclude sampling, a field-valued coupling, a
   discrete gear choice, a tolerance comparable with the bands — determines an
   element of the admissible method rather than merely complicating it (§4);
3. quantitative results on what the formulation and the topology are each worth
   (§5.1, §5.4);
4. two failures of the reliability model, found by widening the uncertain vector
   and by sampling against the first-order estimate, together with the
   corrections they imply (§5.3).

## 1.4 Structure

§2 reviews the methods available for each feature of the problem and identifies
which are admissible. §3 states the mechanism, the specification, the
formulation, and the relaxation the specification turns out to require. §4 derives each element of the method from a property of that
formulation, and describes the implementation and its verification. §5 presents
the results, §6 the conclusions, limitations and further work. Derivations are
collected in {doc}`Appendix A <theory>`, the software in {doc}`Appendix B
<implementation>` and {doc}`api`, and the supporting measurements in
{doc}`Appendix C <supporting>`.

Every quantity reported is computed by the code described in §4.8 and pinned by
its test suite.
