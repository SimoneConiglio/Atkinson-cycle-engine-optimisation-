"""The closed-form inversion of the loop-closure equations."""

from __future__ import annotations

import numpy as np
import pytest

from exlink import PUBLISHED_DESIGN
from exlink.kinematics import crank_angles, solve


@pytest.fixture(scope="module")
def kin():
    return solve(PUBLISHED_DESIGN, samples=720)


def test_the_crank_angles_are_measured_from_the_y_axis() -> None:
    """``theta_1`` and ``theta_2`` are referred to ``+y``, ``theta_r`` to ``+x``.

    The loop closure carries ``q_1 sin(theta_1)`` in its x-projection and
    ``-q_1 cos(theta_1)`` in its y-projection, which is exactly this convention
    written down; nothing else in the module states it, and a figure that
    assumes ``+x`` puts every crank arc a right angle away from the member it
    measures.  Pinned here so the drawings and the derivations cannot drift
    apart from the code.
    """
    design = PUBLISHED_DESIGN
    for degrees in (0.0, 45.0, 90.0):
        angle = np.radians(degrees)
        frame = solve(design, theta_1=np.array([angle]))
        crank = frame.Q[0] - frame.R1[0]
        assert crank == pytest.approx(
            design.q_1 * np.array([-np.sin(angle), np.cos(angle)]), abs=1e-9
        )
        throw = frame.D[0] - frame.R2[0]
        second = float(frame.theta_2[0])
        assert throw == pytest.approx(
            design.q_2 * np.array([-np.sin(second), np.cos(second)]), abs=1e-9
        )
        # theta_r, by contrast, is the direction of R1 -> R2 from +x.
        assert frame.R2[0] == pytest.approx(
            design.I * np.array([np.cos(design.theta_r_rad), np.sin(design.theta_r_rad)]),
            abs=1e-9,
        )


def test_every_link_keeps_its_length(kin) -> None:
    """The whole point of a rigid linkage: no bar may stretch.

    This is the strongest single check on the inversion -- it fails for almost
    any sign or indexing error in the chain.
    """
    design = kin.design
    for start, end, expected, label in (
        (kin.R1, kin.Q, design.q_1, "q_1"),
        (kin.Q, kin.A, design.a, "a"),
        (kin.A, kin.D, design.c, "c"),
        (kin.A, kin.E, design.b, "b"),
        (kin.D, kin.E, design.d, "d"),
        (kin.E, kin.P, design.e, "e"),
        (kin.R2, kin.D, design.q_2, "q_2"),
    ):
        lengths = np.linalg.norm(end - start, axis=1)
        assert lengths == pytest.approx(expected, abs=1e-9), label


def test_the_piston_stays_on_the_cylinder_axis(kin) -> None:
    """``P`` and ``H`` must sit at ``x = x_1`` for every crank angle."""
    assert kin.P[:, 0] == pytest.approx(kin.design.x_1, abs=1e-9)
    assert kin.H[:, 0] == pytest.approx(kin.design.x_1, abs=1e-9)


def test_the_piston_length_is_respected(kin) -> None:
    assert kin.H[:, 1] - kin.P[:, 1] == pytest.approx(kin.spec.piston_length, abs=1e-9)


def test_the_shafts_are_fixed(kin) -> None:
    assert pytest.approx(0.0) == kin.R1
    assert np.ptp(kin.R2, axis=0) == pytest.approx(0.0)
    assert np.linalg.norm(kin.R2[0]) == pytest.approx(kin.design.I)


def test_the_gear_relation_holds(kin) -> None:
    """``theta_2 = -2 theta_1 + theta_f``, equation (1)."""
    expected = -2.0 * kin.theta_1 + kin.design.theta_f_rad
    assert kin.theta_2 == pytest.approx(expected)


def test_the_transmission_angle_matches_its_definition(kin) -> None:
    assert kin.transmission_angle == pytest.approx(kin.theta_a - kin.theta_T)


def test_compatibility_measures_are_the_arccosine_arguments(kin) -> None:
    """``delta_c1`` and ``delta_c2`` are conditions (4a) and (6a)."""
    design = kin.design
    a, c = design.a, design.c
    proj_a = (
        design.q_1 * np.sin(kin.theta_1)
        - design.q_2 * np.sin(kin.theta_2)
        + design.I * np.cos(design.theta_r_rad)
    )
    proj_b = (
        -design.q_1 * np.cos(kin.theta_1)
        + design.q_2 * np.cos(kin.theta_2)
        + design.I * np.sin(design.theta_r_rad)
    )
    cos_t = (proj_a**2 + proj_b**2 - a**2 - c**2) / (2.0 * a * c)
    assert kin.delta_c1 == pytest.approx(np.max(np.abs(cos_t)))
    assert kin.delta_c2 == pytest.approx(np.max(np.abs(np.cos(kin.theta_e))))
    assert kin.compatibility == pytest.approx(max(kin.delta_c1, kin.delta_c2))


def test_the_published_design_is_kinematically_compatible(kin) -> None:
    assert kin.feasible
    assert kin.compatibility < 1.0


def test_an_impossible_linkage_is_flagged_not_raised() -> None:
    """A swing rod far too short to close the loop must not raise.

    The optimizer relies on being handed ``W >= 1`` as a number it can follow
    downhill, rather than an exception it has to catch.
    """
    broken = PUBLISHED_DESIGN.replace(a=25.0, c=25.0)
    result = solve(broken, samples=180)
    assert not result.feasible
    assert result.compatibility >= 1.0


def test_solving_at_explicit_angles(kin) -> None:
    angles = np.radians([0.0, 90.0, 180.0])
    result = solve(PUBLISHED_DESIGN, theta_1=angles)
    assert result.theta_1 == pytest.approx(angles)
    assert result.lam.shape == (3,)


def test_crank_angles_span_one_revolution() -> None:
    angles = crank_angles(360)
    assert angles[0] == 0.0
    assert angles[-1] < 2.0 * np.pi
    assert angles.size == 360


def test_crank_angles_rejects_a_coarse_grid() -> None:
    with pytest.raises(ValueError, match="at least 8"):
        crank_angles(4)


def test_bodies_are_polylines_over_every_angle(kin) -> None:
    bodies = kin.bodies
    n = kin.theta_1.size
    assert bodies["trigonal"].shape == (n, 4, 2)
    assert bodies["piston"].shape == (n, 2, 2)
    # The trigonal polyline closes back on A.
    assert bodies["trigonal"][:, 0] == pytest.approx(bodies["trigonal"][:, 3])
