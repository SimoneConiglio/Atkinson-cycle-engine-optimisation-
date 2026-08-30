"""Spectral differentiation of the periodic histories."""

from __future__ import annotations

import numpy as np
import pytest

from exlink.derivatives import ramp_derivative, spectral_derivative, winding_number

N = 360
THETA = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)


def test_derivative_of_a_trigonometric_polynomial_is_exact() -> None:
    """Spectral differentiation is exact on band-limited data, unlike a difference."""
    values = np.cos(3.0 * THETA) + 0.4 * np.sin(5.0 * THETA)
    expected = -3.0 * np.sin(3.0 * THETA) + 2.0 * np.cos(5.0 * THETA)
    assert spectral_derivative(values, 1) == pytest.approx(expected, abs=1e-10)


def test_second_derivative_is_exact() -> None:
    values = np.cos(3.0 * THETA) + 0.4 * np.sin(5.0 * THETA)
    expected = -9.0 * np.cos(3.0 * THETA) - 10.0 * np.sin(5.0 * THETA)
    assert spectral_derivative(values, 2) == pytest.approx(expected, abs=1e-8)


def test_zeroth_derivative_returns_a_copy() -> None:
    values = np.cos(THETA)
    result = spectral_derivative(values, 0)
    assert result == pytest.approx(values)
    assert result is not values


def test_negative_order_is_rejected() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        spectral_derivative(np.cos(THETA), -1)


def test_winding_number_of_the_gear_relation() -> None:
    """``theta_2 = -2 theta_1 + theta_f`` turns twice backwards per revolution."""
    wrapped = np.angle(np.exp(1j * (-2.0 * THETA + 0.7)))
    assert winding_number(wrapped) == pytest.approx(-2.0, abs=1e-9)


def test_winding_number_of_a_rocking_link_is_zero() -> None:
    assert winding_number(0.4 * np.sin(THETA)) == pytest.approx(0.0, abs=1e-9)


def test_ramp_derivative_handles_a_wrapped_winding_angle() -> None:
    """A wrapped angle that accumulates turns still differentiates exactly."""
    raw = -2.0 * THETA + 0.3 * np.sin(THETA) + 0.7
    wrapped = np.angle(np.exp(1j * raw))
    assert ramp_derivative(wrapped, 1) == pytest.approx(-2.0 + 0.3 * np.cos(THETA), abs=1e-10)
    assert ramp_derivative(wrapped, 2) == pytest.approx(-0.3 * np.sin(THETA), abs=1e-8)


def test_derivative_works_along_the_last_axis() -> None:
    stacked = np.stack([np.cos(THETA), np.sin(THETA)])
    result = spectral_derivative(stacked, 1)
    assert result[0] == pytest.approx(-np.sin(THETA), abs=1e-10)
    assert result[1] == pytest.approx(np.cos(THETA), abs=1e-10)
