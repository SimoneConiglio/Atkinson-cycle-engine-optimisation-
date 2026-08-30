"""The 11-dimensional design vector ``X`` of the EX-link mechanism.

The report parametrises the linkage by

.. math:: X = (a, c, I, x_b, y_b, x_1, e, q_1, q_2, \\theta_f, \\theta_r)^T

Rather than describing the trigonal link by its three sides ``b``, ``c``, ``d``
-- which must satisfy the triangle inequality and leave the sign of
``theta_b`` undetermined -- the report places ``E`` in the frame carried by
``c = AD``.  ``x_b`` and ``y_b`` are then free to take any value, positive or
negative, and ``theta_b = atan2(y_b, x_b)`` carries its own sign.  That
reparametrisation is what makes the design space a plain box.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import astuple, dataclass, fields
from typing import ClassVar

import numpy as np
from numpy.typing import ArrayLike, NDArray

#: Order of the design variables in the flat vector representation.
VARIABLE_NAMES: tuple[str, ...] = (
    "a",
    "c",
    "I",
    "x_b",
    "y_b",
    "x_1",
    "e",
    "q_1",
    "q_2",
    "theta_f",
    "theta_r",
)

#: Human-readable description of each design variable.
VARIABLE_DESCRIPTIONS: dict[str, str] = {
    "a": "length of the swing rod QA [mm]",
    "c": "side AD of the trigonal link [mm]",
    "I": "distance between crankshaft and eccentric shaft [mm]",
    "x_b": "abscissa of E in the frame carried by AD [mm]",
    "y_b": "ordinate of E in the frame carried by AD [mm]",
    "x_1": "lateral offset of the cylinder axis from the crankshaft [mm]",
    "e": "length of the piston rod EP [mm]",
    "q_1": "length of the crank fixed to the crankshaft [mm]",
    "q_2": "length of the crank fixed to the eccentric shaft [mm]",
    "theta_f": "dephasing between the two cranks at theta_1 = 0 [deg]",
    "theta_r": "orientation of the crankshaft-to-eccentric-shaft axis [deg]",
}

#: The two entries of the design vector that are angles, stored in degrees.
ANGULAR_VARIABLES: frozenset[str] = frozenset({"theta_f", "theta_r"})


@dataclass(frozen=True)
class Design:
    """A single set of mechanism dimensions.

    Angles are stored in **degrees** because that is how they are reported and
    bounded; every consumer converts to radians through
    :attr:`theta_f_rad` / :attr:`theta_r_rad`.
    """

    a: float
    c: float
    I: float
    x_b: float
    y_b: float
    x_1: float
    e: float
    q_1: float
    q_2: float
    theta_f: float
    theta_r: float

    names: ClassVar[tuple[str, ...]] = VARIABLE_NAMES

    # -- alternative constructors -------------------------------------------------

    @classmethod
    def from_array(cls, values: ArrayLike) -> Design:
        """Build a design from a flat vector ordered as :data:`VARIABLE_NAMES`."""
        array = np.asarray(values, dtype=float).ravel()
        if array.size != len(VARIABLE_NAMES):
            msg = f"expected {len(VARIABLE_NAMES)} design variables, got {array.size}"
            raise ValueError(msg)
        return cls(*array.tolist())

    @classmethod
    def from_mapping(cls, mapping: dict[str, ArrayLike]) -> Design:
        """Build a design from a ``{name: value}`` mapping (GEMSEO input data)."""
        values = [float(np.ravel(mapping[name])[0]) for name in VARIABLE_NAMES]
        return cls(*values)

    # -- conversions --------------------------------------------------------------

    def to_array(self) -> NDArray[np.float64]:
        """Return the design as a flat vector ordered as :data:`VARIABLE_NAMES`."""
        return np.array(astuple(self), dtype=float)

    def to_mapping(self) -> dict[str, NDArray[np.float64]]:
        """Return the design as GEMSEO-style ``{name: array([value])}`` data."""
        return {name: np.array([getattr(self, name)], dtype=float) for name in VARIABLE_NAMES}

    def replace(self, **changes: float) -> Design:
        """Return a copy with some variables overridden."""
        current = {field.name: getattr(self, field.name) for field in fields(self)}
        unknown = set(changes) - set(current)
        if unknown:
            msg = f"unknown design variables: {sorted(unknown)}"
            raise ValueError(msg)
        current.update(changes)
        return Design(**current)

    # -- derived geometry ---------------------------------------------------------

    @property
    def theta_f_rad(self) -> float:
        """Crank dephasing ``theta_f`` in radians."""
        return math.radians(self.theta_f)

    @property
    def theta_r_rad(self) -> float:
        """Shaft-axis orientation ``theta_r`` in radians."""
        return math.radians(self.theta_r)

    @property
    def b(self) -> float:
        """Side ``AE`` of the trigonal link, ``b = hypot(x_b, y_b)`` [mm]."""
        return math.hypot(self.x_b, self.y_b)

    @property
    def theta_b(self) -> float:
        """Signed angle ``atan2(y_b, x_b)`` between ``AE`` and ``AD`` [rad].

        Equivalent to the Carnot expression
        ``arccos((b^2 + c^2 - d^2) / (2 b c))`` of the report, but with a sign.
        """
        return math.atan2(self.y_b, self.x_b)

    @property
    def d(self) -> float:
        """Side ``DE`` of the trigonal link [mm]."""
        return math.hypot(self.x_b - self.c, self.y_b)

    @property
    def r_1(self) -> float:
        """Primitive radius of the crankshaft gear, ``r_1 = 2 I / 3`` [mm]."""
        return 2.0 * self.I / 3.0

    @property
    def r_2(self) -> float:
        """Primitive radius of the eccentric-shaft gear, ``r_2 = I / 3`` [mm]."""
        return self.I / 3.0

    def __str__(self) -> str:  # pragma: no cover - presentation only
        parts = []
        for name in VARIABLE_NAMES:
            unit = "deg" if name in ANGULAR_VARIABLES else "mm"
            parts.append(f"{name}={getattr(self, name):.4g} {unit}")
        return "Design(" + ", ".join(parts) + ")"


@dataclass(frozen=True)
class Bounds:
    """Box bounds ``l_b <= X <= u_b`` on the design vector."""

    lower: NDArray[np.float64]
    upper: NDArray[np.float64]

    def __post_init__(self) -> None:
        lower = np.asarray(self.lower, dtype=float).ravel()
        upper = np.asarray(self.upper, dtype=float).ravel()
        if lower.size != len(VARIABLE_NAMES) or upper.size != len(VARIABLE_NAMES):
            msg = f"bounds must have {len(VARIABLE_NAMES)} entries"
            raise ValueError(msg)
        if np.any(lower > upper):
            bad = [n for n, lo, up in zip(VARIABLE_NAMES, lower, upper, strict=True) if lo > up]
            msg = f"lower bound above upper bound for {bad}"
            raise ValueError(msg)
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)

    @classmethod
    def around(
        cls,
        design: Design,
        relative: float = 0.1,
        absolute_angle: float = 20.0,
    ) -> Bounds:
        """Bounds centred on ``design``, as used for the local Pareto fronts.

        The report shrinks the box to ``0.9 X_0 <= X <= 1.1 X_0`` before running
        a MOEA seeded near a known solution.  A purely multiplicative rule
        collapses on variables that may be negative or near zero (``x_1``,
        ``y_b``, the angles), so angles get a symmetric absolute window instead.

        Args:
            design: The centre of the box.
            relative: Half-width as a fraction of ``|X_0|`` for length variables.
            absolute_angle: Half-width in degrees for ``theta_f`` and ``theta_r``.
        """
        centre = design.to_array()
        half = np.abs(centre) * relative
        for index, name in enumerate(VARIABLE_NAMES):
            if name in ANGULAR_VARIABLES:
                half[index] = absolute_angle
        return cls(lower=centre - half, upper=centre + half)

    def to_array(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return ``(lower, upper)`` as a pair of arrays."""
        return self.lower, self.upper

    def clip(self, design: Design) -> Design:
        """Project a design back into the box."""
        return Design.from_array(np.clip(design.to_array(), self.lower, self.upper))

    def contains(self, design: Design) -> bool:
        """Whether a design lies inside the box."""
        values = design.to_array()
        return bool(np.all(values >= self.lower) and np.all(values <= self.upper))

    def items(self) -> Iterable[tuple[str, float, float]]:
        """Iterate over ``(name, lower, upper)`` triples."""
        return zip(VARIABLE_NAMES, self.lower, self.upper, strict=True)


#: Global bounds used by the global searches.
#:
#: The report never tabulates ``l_b`` and ``u_b``; these are chosen to comfortably
#: contain its published solution while keeping every length physically sensible
#: for a single-cylinder engine of 74 mm expansion stroke.
GLOBAL_BOUNDS = Bounds(
    lower=np.array([20.0, 20.0, 20.0, -200.0, -200.0, -100.0, 40.0, 2.0, 2.0, -180.0, -180.0]),
    upper=np.array([250.0, 250.0, 150.0, 250.0, 200.0, 100.0, 300.0, 60.0, 80.0, 180.0, 180.0]),
)
