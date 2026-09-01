# 1. Introduction

## 1.1 The problem

An engine for a fuel-economy competition is scored on one quantity: the distance
it carries a car on a fixed quantity of fuel. Designing its mechanism is
therefore a problem with a single scalar objective, and the difficulty is that
almost nothing about the mechanism affects that objective directly.

The linkage studied here is an **extended-expansion** mechanism. It expands the
burnt gas through a larger volume ratio than it compressed the fresh charge,
recovering work an Otto cycle sends out of the exhaust. Achieving that
mechanically requires a piston that reaches top dead centre twice per crankshaft
revolution with two *different* bottom dead centres — the short one setting the
compression stroke, the long one the expansion stroke.

The design freedom is eleven linkage dimensions plus a gear pair. What the
designer wants to know is how those choices reach the fuel consumption of a car.

## 1.2 The stakes

They reach it by four routes, and the routes conflict:

**Thermodynamic.** The expansion-to-compression ratio sets how much of the heat
release becomes indicated work. More expansion is better, and it requires a
longer stroke.

**Mechanical.** Every joint reaction and every unit of piston side load is
friction, which is indicated work that never reaches the crankshaft. A linkage
with a long lever arm produces high joint loads.

**Inertial.** The parts have mass, and at speed their accelerations load the
same joints. The masses depend on the sections, the sections depend on the
loads, and the loads depend on the masses: a closed loop.

**Vehicle-level.** Engine mass is carried for the whole distance, and rolling
resistance is proportional to it. A heavier engine that is thermodynamically
better may still lose.

A formulation that omits any of these gets the answer wrong in a specific,
predictable way, and §2 shows that the formulation conventionally used for this
mechanism omits three of the four.

## 1.3 What decides whether the problem is solvable

Two structural properties, both established in §3, determine which methods can
be applied at all:

- Two of the seven requirements are **equalities**, which makes the feasible set
  a measure-zero manifold in the design space. No method that proceeds by
  sampling can find a feasible point.
- The equalities must therefore be relaxed into **tolerance bands**, and those
  bands turn out to be comparable in width with the manufacturing scatter of the
  parts. A deterministic answer inside the band is not the same as an engine
  that meets the requirement.

Neither property is a numerical inconvenience to be worked around. The first
dictates the choice of optimizer; the second changes what "feasible" means and
obliges a probabilistic treatment.

## 1.4 Scope and structure

This document states the problem (§3), reviews the methods available for each of
its features (§2), derives the methodology (§3), describes the implementation
(§4), presents the use case (§5), reports and discusses the results (§6), and
sets out what the framework cannot yet do (§7). Derivations are collected in the
{doc}`theory` appendix; every module is documented in {doc}`api`.

Every quantity reported is computed by the code described in §4 and pinned by
its test suite.

---

Next: [2. State of the art](state_of_the_art.md)
