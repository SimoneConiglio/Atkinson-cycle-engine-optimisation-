"""Coupling strength, measured rather than asserted."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.formulations import (
    CouplingStrength,
    coupling_curve,
    coupling_strength,
    format_coupling,
)
from exlink.reference import COUPLED_DESIGN


def test_at_rest_the_disciplines_are_decoupled() -> None:
    """With no inertia there is no feedback, and the loop is not a loop.

    This is the sharpest available check that the coupling measure means what
    it claims: at zero speed the sections cannot influence the loads at all,
    because the only path from mass to load is through the inertia forces.  The
    quasi-static problem of the geometric formulation has to be recovered
    exactly, and a measure that reported a non-zero coupling there would be
    measuring its own solver instead of the physics.
    """
    measured = coupling_strength(COUPLED_DESIGN, speed_rpm=0.0)
    assert measured.rho < 1.0e-6
    assert measured.descriptor == "weak"


def test_coupling_grows_with_speed() -> None:
    """The loop gain scales with the inertia forces, so with ``omega^2``."""
    curve = [item for item in coupling_curve(COUPLED_DESIGN) if item.converged]
    values = [item.rho for item in curve]
    assert values == sorted(values)
    assert values[-1] > 0.4


def test_a_strong_coupling_is_reported_as_one() -> None:
    """The claim that this problem needs an MDA has to be falsifiable.

    A contraction factor above 0.5 means each discipline rewrites a large part
    of the other's input on every sweep, and a single sequential pass would be
    badly wrong.  Below 0.1 the honest thing would be to say so and drop the
    MDA.
    """
    measured = coupling_strength(COUPLED_DESIGN, speed_rpm=1500.0)
    assert measured.descriptor == "strong"
    assert measured.sweeps > 10


def test_sweeps_per_decade_matches_the_measured_rate() -> None:
    measured = coupling_strength(COUPLED_DESIGN, speed_rpm=1000.0)
    assert measured.sweeps_per_decade == pytest.approx(
        -1.0 / np.log10(measured.rho), rel=1.0e-9
    )


def test_an_unanalysable_design_is_reported_as_divergent() -> None:
    broken = COUPLED_DESIGN.replace(a=0.4 * COUPLED_DESIGN.a)
    measured = coupling_strength(broken, speed_rpm=1000.0)
    assert not measured.converged
    assert measured.descriptor == "divergent"
    assert not np.isfinite(measured.sweeps_per_decade)


def test_format_coupling_renders_every_row() -> None:
    curve = coupling_curve(COUPLED_DESIGN, speeds=(0.0, 1000.0))
    rendered = format_coupling(curve)
    assert "rho" in rendered
    assert len(rendered.splitlines()) == 6


def test_descriptor_thresholds_are_monotone() -> None:
    words = [CouplingStrength(0.0, rho, 1, True).descriptor for rho in (0.01, 0.3, 0.8, 1.2)]
    assert words == ["weak", "moderate", "strong", "divergent"]


def test_the_coupling_is_load_histories_not_scalars() -> None:
    """Why the formulation question has a structural answer here.

    IDF's cost scales with the dimension of the coupling, and the coupling in
    this problem is not a handful of scalars: ``member_axial`` and
    ``member_bending`` are the internal load history of every member, at every
    crank angle, at every station.  Four orders of magnitude more coupling
    variables than design variables means IDF cannot be posed, never mind won.
    """
    from exlink.disciplines import COUPLED_SAMPLES
    from exlink.dynamics import MEMBER_NAMES
    from exlink.formulations import coupling_dimension
    from exlink.sizing import STATIONS

    sizes = coupling_dimension()
    history = len(MEMBER_NAMES) * COUPLED_SAMPLES * STATIONS
    assert sizes["member_axial"] == history
    assert sizes["member_bending"] == history
    assert sizes["total"] > 1000 * sizes["design_variables"]


def test_idf_reports_its_own_refusal_rather_than_raising() -> None:
    """A formulation that cannot be built is a comparison outcome, not an error.

    The run has to come back with the reason attached, so the comparison can
    state *why* IDF is unavailable instead of simply omitting it.
    """
    from exlink.formulations import compare_formulations, format_formulations
    from exlink.reference import COUPLED_DESIGN

    rows = compare_formulations(
        initial=COUPLED_DESIGN, speed_rpm=0.0, max_iter=2, names=("IDF",)
    )
    assert len(rows) == 1
    assert not rows[0].feasible
    assert rows[0].error
    rendered = format_formulations(rows)
    assert "coupling variables" in rendered
