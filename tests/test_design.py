"""The design vector and its derived trigonal-link geometry."""

from __future__ import annotations

import math

import numpy as np
import pytest

from exlink import GLOBAL_BOUNDS, VARIABLE_NAMES, Bounds, Design


def test_round_trips_through_an_array(published: Design) -> None:
    assert Design.from_array(published.to_array()) == published


def test_round_trips_through_a_mapping(published: Design) -> None:
    assert Design.from_mapping(published.to_mapping()) == published


def test_array_order_matches_the_declared_names(published: Design) -> None:
    values = published.to_array()
    for index, name in enumerate(VARIABLE_NAMES):
        assert values[index] == pytest.approx(getattr(published, name))


def test_rejects_a_vector_of_the_wrong_length() -> None:
    with pytest.raises(ValueError, match="11 design variables"):
        Design.from_array([1.0, 2.0, 3.0])


def test_replace_rejects_unknown_variables(published: Design) -> None:
    with pytest.raises(ValueError, match="unknown design variables"):
        published.replace(nonsense=1.0)


def test_trigonal_reparametrisation_satisfies_carnot(published: Design) -> None:
    """``theta_b`` from ``atan2(y_b, x_b)`` must agree with the Carnot theorem.

    ``(x_b, y_b)`` is introduced precisely so that ``b``, ``c`` and ``d``
    always close a triangle; this checks the two routes to ``theta_b`` agree.
    """
    b, c, d = published.b, published.c, published.d
    carnot = math.acos((b**2 + c**2 - d**2) / (2.0 * b * c))
    assert abs(published.theta_b) == pytest.approx(carnot, abs=1e-12)


def test_trigonal_sides_are_consistent(published: Design) -> None:
    assert published.b == pytest.approx(math.hypot(published.x_b, published.y_b))
    assert published.d == pytest.approx(math.hypot(published.x_b - published.c, published.y_b))


def test_gear_radii_follow_the_transmission_ratio(published: Design) -> None:
    assert published.r_1 + published.r_2 == pytest.approx(published.I)
    assert published.r_1 / published.r_2 == pytest.approx(2.0)


def test_bounds_around_a_design_contain_it(published: Design) -> None:
    box = Bounds.around(published, relative=0.1)
    assert box.contains(published)


def test_bounds_around_gives_angles_an_absolute_window(published: Design) -> None:
    """A multiplicative window collapses on angles that may be near zero."""
    box = Bounds.around(published, relative=0.1, absolute_angle=20.0)
    index = VARIABLE_NAMES.index("theta_r")
    assert box.upper[index] - box.lower[index] == pytest.approx(40.0)


def test_bounds_reject_a_crossed_box() -> None:
    with pytest.raises(ValueError, match="lower bound above upper bound"):
        Bounds(lower=np.ones(11), upper=np.zeros(11))


def test_clip_projects_into_the_box(published: Design) -> None:
    box = Bounds.around(published, relative=0.01)
    outside = published.replace(a=published.a * 3.0)
    assert box.contains(box.clip(outside))


def test_the_published_design_lies_inside_the_global_bounds(published: Design) -> None:
    assert GLOBAL_BOUNDS.contains(published)
