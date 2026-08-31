"""Manufacturability rules: preferred sizes, floors, and what they cost."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.gears import lattice_inter_axle, lattice_neighbours, size_pair, tooth_count
from exlink.manufacturing import (
    MIN_TEETH,
    STANDARD_MODULES,
    STOCK_DIAMETERS,
    round_to_module,
    round_up_to_stock,
    stock_premium,
)


def test_rounding_never_weakens_a_member() -> None:
    """Stock rounding must go up, always.

    Rounding down would silently turn a section the sizing loop certified as
    safe into one that is not, which is the one failure mode a manufacturing
    layer must never introduce.
    """
    required = np.linspace(3.0, 95.0, 400)
    supplied = round_up_to_stock(required)
    assert np.all(supplied >= required - 1.0e-12)


def test_rounding_lands_on_the_catalogue() -> None:
    """Every rounded value is a size a supplier actually stocks."""
    supplied = round_up_to_stock(np.linspace(3.0, 99.0, 200))
    assert set(np.unique(supplied)).issubset(set(STOCK_DIAMETERS.tolist()))


def test_oversize_requirement_is_not_silently_clipped() -> None:
    """A runaway section stays visible as a runaway.

    Clipping to the largest stock size would report a buildable design for a
    mechanism that cannot be built at all.
    """
    assert round_up_to_stock(500.0) == pytest.approx(500.0)


def test_exact_stock_size_is_left_alone() -> None:
    for size in STOCK_DIAMETERS:
        assert round_up_to_stock(float(size)) == pytest.approx(size)


def test_stock_premium_is_zero_on_the_catalogue_and_positive_off_it() -> None:
    exact = [8.0, 12.0, 20.0]
    assert stock_premium(exact, exact) == pytest.approx(0.0)
    assert stock_premium([8.1], round_up_to_stock([8.1])) > 0.0


def test_module_rounds_to_the_standard_series() -> None:
    for value in (0.55, 1.1, 1.4, 2.3, 7.0):
        assert round_to_module(value) in set(STANDARD_MODULES.tolist())


def test_inter_axle_distance_lives_on_a_lattice() -> None:
    """``I = 1.5 m z`` is the whole point of treating the module as discrete.

    A continuous ``I`` that is not on the lattice describes a mechanism whose
    gears cannot be cut, so every candidate returned must reproduce its own
    module and tooth count exactly.
    """
    for module, teeth, value in lattice_neighbours(56.5)[:10]:
        assert value == pytest.approx(1.5 * module * teeth)
        assert teeth >= MIN_TEETH
        assert tooth_count(value, module) == teeth


def test_lattice_neighbours_bracket_the_request() -> None:
    """The candidate list must contain a distance on each side."""
    target = 56.5
    values = [value for _m, _z, value in lattice_neighbours(target)]
    assert min(values) < target < max(values)


def test_lattice_neighbours_are_ordered_by_distance() -> None:
    target = 61.37
    errors = [abs(value - target) for _m, _z, value in lattice_neighbours(target)]
    assert errors == sorted(errors)


def test_undercut_limit_is_respected() -> None:
    """A 20 deg pinion below 17 teeth is undercut and must not be offered."""
    assert tooth_count(10.0, 3.0) == MIN_TEETH
    assert all(teeth >= MIN_TEETH for _m, teeth, _v in lattice_neighbours(20.0))


def test_gear_face_width_carries_the_load() -> None:
    """Doubling the tooth load must not leave the pair overstressed."""
    light = size_pair(56.5, 400.0, module=1.5)
    heavy = size_pair(56.5, 800.0, module=1.5)
    assert heavy.face_width >= light.face_width
    assert max(heavy.bending_utilisation, heavy.contact_utilisation) <= 1.0 + 1.0e-9


def test_gear_geometry_is_self_consistent() -> None:
    pair = size_pair(56.5, 600.0, module=2.0)
    assert pair.teeth_large == 2 * pair.teeth_small
    assert pair.radius_large == pytest.approx(2.0 * pair.radius_small)
    assert pair.inter_axle == pytest.approx(pair.radius_small + pair.radius_large)
    assert pair.inter_axle == pytest.approx(lattice_inter_axle(pair.module, pair.teeth_small))


def test_automatic_module_choice_is_feasible_and_light() -> None:
    """Left to choose, the sizer must return a workable pair.

    And it must not simply take the largest module: that would always work
    and would always be heavy.
    """
    chosen = size_pair(56.5, 600.0)
    assert chosen.feasible
    assert chosen.module < STANDARD_MODULES[-1]
    forced = size_pair(56.5, 600.0, module=float(STANDARD_MODULES[-1]))
    assert chosen.mass <= forced.mass
