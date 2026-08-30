"""Material data and fatigue-strength corrections for component sizing.

Unit system
-----------
The structural code uses the consistent set **mm, N, MPa, tonne, s**, so that
``F = m a`` needs no conversion factor: one tonne accelerated at one mm/s^2 is
one newton.  Steel density is therefore ``7.85e-9 t/mm^3``.  Masses are exposed
in kilograms at the API surface, where a reader expects them.

Fatigue model
-------------
Rotating-bending endurance limit corrected by Marin factors, then a Goodman
line, which is the standard textbook treatment (Shigley) and the right level of
fidelity for a first sizing iteration:

.. math::
    \\frac{\\sigma_a}{S_e} + \\frac{\\sigma_m}{S_u} \\le \\frac{1}{n_f}

with ``S_e = k_a k_b k_c k_d k_e S_e'`` and ``S_e' = 0.5 S_u``.  The size factor
``k_b`` depends on the diameter being solved for, which is one reason the sizing
of each member is a small root-find rather than a closed form.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]

STEEL_DENSITY = 7.85e-9
"""Density of steel in the consistent unit system [tonne/mm^3]."""


@dataclass(frozen=True)
class Material:
    """An isotropic metal, described well enough to size a part against it."""

    name: str = "42CrMo4 quenched and tempered"
    density: float = STEEL_DENSITY
    """[tonne/mm^3]"""

    youngs_modulus: float = 210_000.0
    """``E`` [MPa]."""

    yield_strength: float = 700.0
    """``S_y`` [MPa]."""

    ultimate_strength: float = 900.0
    """``S_u`` [MPa]."""

    surface_factor_a: float = 4.51
    surface_factor_b: float = -0.265
    """``k_a = a S_u^b`` for a machined surface (Shigley, MPa units)."""

    reliability_factor: float = 0.814
    """``k_e`` for 99 % reliability."""

    temperature_factor: float = 1.0
    """``k_d``; 1.0 at ambient."""

    @property
    def density_kg_per_mm3(self) -> float:
        """Density in kg/mm^3, for reporting masses in kilograms."""
        return self.density * 1000.0

    @property
    def endurance_limit_uncorrected(self) -> float:
        """``S_e' = 0.5 S_u`` [MPa], capped at 700 MPa as usual for steels."""
        return min(0.5 * self.ultimate_strength, 700.0)

    @property
    def surface_factor(self) -> float:
        """``k_a``, the machined-surface correction."""
        return self.surface_factor_a * self.ultimate_strength**self.surface_factor_b

    def size_factor(self, diameter: FloatArray | float) -> FloatArray:
        """``k_b``, the size correction for a solid round section.

        Args:
            diameter: Section diameter [mm].

        Returns:
            ``k_b``, elementwise.
        """
        d = np.atleast_1d(np.asarray(diameter, dtype=float))
        factor = np.ones_like(d)
        small = (d >= 2.79) & (d <= 51.0)
        large = d > 51.0
        factor[small] = 1.24 * d[small] ** -0.107
        factor[large] = 1.51 * d[large] ** -0.157
        return factor

    def endurance_limit(
        self, diameter: FloatArray | float, load_factor: float = 1.0
    ) -> FloatArray:
        """Corrected endurance limit ``S_e`` [MPa].

        Args:
            diameter: Section diameter [mm], entering through ``k_b``.
            load_factor: ``k_c``; 1.0 bending, 0.85 axial, 0.59 torsion.  Use
                1.0 where bending dominates, which it does for every link here.

        Returns:
            ``S_e``, elementwise in ``diameter``.
        """
        return (
            self.surface_factor
            * self.size_factor(diameter)
            * load_factor
            * self.temperature_factor
            * self.reliability_factor
            * self.endurance_limit_uncorrected
        )


@dataclass(frozen=True)
class SafetyFactors:
    """Design factors applied to each failure mode."""

    static: float = 1.5
    """``n_y`` against first yield."""

    fatigue: float = 2.0
    """``n_f`` on the Goodman line."""

    buckling: float = 3.0
    """``n_b`` against Euler buckling of members in compression.

    Higher than the others because Euler buckling is sudden, and because a
    real link has imperfections the ideal formula ignores.
    """


DEFAULT_MATERIAL = Material()
DEFAULT_SAFETY = SafetyFactors()


def goodman_utilisation(
    alternating: FloatArray,
    mean: FloatArray,
    endurance: FloatArray,
    ultimate: float,
) -> FloatArray:
    """Goodman utilisation; ``<= 1 / n_f`` is safe.

    Args:
        alternating: ``sigma_a`` [MPa], non-negative.
        mean: ``sigma_m`` [MPa], signed.
        endurance: ``S_e`` [MPa].
        ultimate: ``S_u`` [MPa].

    Returns:
        ``sigma_a / S_e + max(sigma_m, 0) / S_u``.

    A compressive mean stress is not credited as beneficial: the Goodman line is
    truncated at ``sigma_m = 0``, which is the conservative convention and
    matters here because the connecting members swing between tension and
    compression every revolution.
    """
    return alternating / endurance + np.maximum(mean, 0.0) / ultimate
