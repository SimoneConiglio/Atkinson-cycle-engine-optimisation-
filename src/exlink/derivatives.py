"""Derivatives of the linkage's periodic histories with respect to crank angle.

Every history the analysis produces -- joint positions, link angles, piston
height -- is a smooth function of ``theta_1`` sampled uniformly over exactly one
revolution.  For such data the Fourier (spectral) derivative is exact to machine
precision, where a finite difference on the same grid would be accurate only to
``O(h^2)``.  That matters because accelerations are *second* derivatives, and a
finite-difference second derivative on a 0.5 deg grid loses about six digits --
enough noise to pollute the inertia forces and, through them, the sizing loop.

Two cases arise:

* **Periodic** histories -- joint coordinates, ``lambda`` -- return to their
  starting value after one revolution and can be differentiated directly.
* **Ramped** histories -- link angles -- may accumulate a whole number of turns
  per revolution, as ``theta_2 = -2 theta_1 + theta_f`` does.  Such a history is
  a linear ramp plus a periodic part; :func:`ramp_derivative` splits the two,
  differentiates the periodic part spectrally, and adds the ramp's slope back.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def spectral_derivative(values: FloatArray, order: int = 1) -> FloatArray:
    """Differentiate a ``2 pi``-periodic history sampled uniformly.

    Args:
        values: Samples over ``[0, 2 pi)``, endpoint excluded, shaped ``(n,)``
            or ``(..., n)`` -- the last axis is the crank angle.
        order: Derivative order; 1 for ``d/dtheta_1``, 2 for ``d^2/dtheta_1^2``.

    Returns:
        The derivative, same shape as ``values``.

    Raises:
        ValueError: If ``order`` is negative.
    """
    if order < 0:
        msg = "derivative order must be non-negative"
        raise ValueError(msg)
    if order == 0:
        return np.asarray(values, dtype=float).copy()

    array = np.asarray(values, dtype=float)
    n = array.shape[-1]
    wavenumbers = np.fft.rfftfreq(n, d=1.0 / n)
    factor = (1j * wavenumbers) ** order
    if order % 2 == 1 and n % 2 == 0:
        # The Nyquist mode carries no reliable sign for an odd derivative; the
        # usual remedy is to drop it rather than let it alias.
        factor[-1] = 0.0
    return np.fft.irfft(np.fft.rfft(array, axis=-1) * factor, n=n, axis=-1)


def winding_number(angles: FloatArray) -> float:
    """Whole turns an angle history accumulates over one revolution.

    Args:
        angles: An angle history over ``[0, 2 pi)``, possibly wrapped.

    Returns:
        The slope ``k`` of the linear ramp, so that ``angles - k theta_1`` is
        periodic.  ``-2`` for the eccentric shaft, ``0`` for a link that rocks.
    """
    unwrapped = np.unwrap(np.asarray(angles, dtype=float))
    # Continue the unwrap one sample past the end, back onto the first sample.
    step = np.angle(np.exp(1j * (unwrapped[0] - unwrapped[-1])))
    total = unwrapped[-1] + step - unwrapped[0]
    return float(total / (2.0 * np.pi))


def ramp_derivative(angles: FloatArray, order: int = 1) -> FloatArray:
    """Differentiate an angle history that may accumulate whole turns.

    Args:
        angles: An angle history over ``[0, 2 pi)``, possibly wrapped.
        order: Derivative order.

    Returns:
        ``d^order angles / d theta_1^order``.
    """
    array = np.asarray(angles, dtype=float)
    n = array.shape[-1]
    theta = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    slope = winding_number(array)
    periodic = np.unwrap(array) - slope * theta
    derivative = spectral_derivative(periodic, order=order)
    if order == 1:
        derivative = derivative + slope
    return derivative


def rate_and_acceleration(values: FloatArray, speed: float) -> tuple[FloatArray, FloatArray]:
    """Convert crank-angle derivatives into time derivatives.

    The mechanism is analysed at constant crankshaft speed, so
    ``d/dt = Omega d/dtheta_1`` and ``d^2/dt^2 = Omega^2 d^2/dtheta_1^2``
    exactly -- no ``theta_1``-dot-dot term.

    Args:
        values: A periodic history over one revolution.
        speed: Crankshaft speed ``Omega`` [rad/s].

    Returns:
        ``(first, second)`` time derivatives.
    """
    return (
        speed * spectral_derivative(values, order=1),
        speed**2 * spectral_derivative(values, order=2),
    )
