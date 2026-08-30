"""Closed-form kinematics of the EX-link linkage.

The report inverts the loop-closure equations analytically rather than solving
them with Newton-Raphson.  That choice matters twice over: the evaluation is
fast enough to sit inside an evolutionary algorithm, and the two arccosine
arguments it exposes -- ``delta_c1`` (4a) and ``delta_c2`` (6a) -- become
explicit *compatibility* measures.  A design whose ``|argument|`` reaches 1 for
some crank angle cannot turn continuously: the crankshaft only rocks.  Feeding
those two numbers to the optimizer as a constraint, instead of letting the
analysis fail, is what makes the problem globally solvable.

Solution chain for a given crank angle ``theta_1``:

1. ``theta_2 = -2 theta_1 + theta_f``                                   (1)
2. ``(A, B)``, the closure of the four-bar ``R1 Q A D R2``            (3a)
3. ``T = arccos((A^2 + B^2 - a^2 - c^2) / (2 a c))``                    (4)
4. ``q = atan2(a sin T, a cos T + c)``                                 (2a)
5. ``theta_T = atan2(B, A) - q`` and ``theta_a = theta_T + T``          (5)
6. ``theta_e = arccos((q_1 sin theta_1 - a cos theta_a
   - b cos(theta_b + theta_T) + x_1) / e)``                            (6)
7. ``lambda = q_1 cos theta_1 + a sin theta_a
   + b sin(theta_b + theta_T) + e sin theta_e + p``                     (7)

Every inverted cosine keeps the positive root, as the report prescribes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .constants import DEFAULT_SPEC, EngineSpec
from .design import Design

FloatArray = NDArray[np.float64]

DEFAULT_SAMPLES = 720
"""Crank angles per revolution used by default.

720 samples (0.5 deg) resolve the two top dead centres finely enough that the
parabolic refinement in :mod:`exlink.metrics` lands within a micrometre.
"""


def crank_angles(samples: int = DEFAULT_SAMPLES) -> FloatArray:
    """Return ``samples`` crank angles spanning one revolution, endpoint excluded."""
    if samples < 8:
        msg = "at least 8 crank angles are needed to detect the four phases"
        raise ValueError(msg)
    return np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False)


@dataclass(frozen=True)
class Kinematics:
    """Pose of every body over one crankshaft revolution.

    All array attributes are shaped ``(n_angles,)`` unless stated otherwise;
    point attributes are shaped ``(n_angles, 2)``.
    """

    design: Design
    spec: EngineSpec
    theta_1: FloatArray
    theta_2: FloatArray
    theta_T: FloatArray
    theta_a: FloatArray
    theta_e: FloatArray
    lam: FloatArray
    """Piston crown height ``lambda`` above the crankshaft axis [mm]."""

    transmission_angle: FloatArray
    """``T = theta_a - theta_T``, the swing-rod / trigonal transmission angle."""

    delta_c1: float
    """``max |cos T|`` over the revolution -- compatibility condition (4a)."""

    delta_c2: float
    """``max |cos theta_e|`` over the revolution -- compatibility condition (6a)."""

    feasible: bool
    """Whether the crankshaft completes a full revolution.

    ``False`` as soon as either compatibility argument reaches 1 somewhere, in
    which case the angles below are computed from clipped arguments and are
    physically meaningless.
    """

    # -- joint trajectories -------------------------------------------------------
    R1: FloatArray
    R2: FloatArray
    Q: FloatArray
    A: FloatArray
    D: FloatArray
    E: FloatArray
    P: FloatArray
    H: FloatArray

    @property
    def compatibility(self) -> float:
        """``W = max(delta_c1, delta_c2)`` of the report."""
        return max(self.delta_c1, self.delta_c2)

    @property
    def bodies(self) -> dict[str, FloatArray]:
        """Named polylines, each shaped ``(n_angles, n_points, 2)``.

        Used by the plotting and animation helpers and by the bounding-box
        computation, so that both see exactly the same geometry.
        """
        return {
            "crank_1": np.stack([self.R1, self.Q], axis=1),
            "swing_rod": np.stack([self.Q, self.A], axis=1),
            "trigonal": np.stack([self.A, self.D, self.E, self.A], axis=1),
            "crank_2": np.stack([self.R2, self.D], axis=1),
            "piston_rod": np.stack([self.E, self.P], axis=1),
            "piston": np.stack([self.P, self.H], axis=1),
        }


def solve(
    design: Design,
    samples: int = DEFAULT_SAMPLES,
    spec: EngineSpec = DEFAULT_SPEC,
    theta_1: FloatArray | None = None,
) -> Kinematics:
    """Solve the linkage over one crankshaft revolution.

    Args:
        design: The mechanism dimensions.
        samples: Number of crank angles, ignored when ``theta_1`` is given.
        spec: Fixed engine data (only ``piston_length`` is used here).
        theta_1: Explicit crank angles [rad]; defaults to a uniform revolution.

    Returns:
        The pose of every body at every crank angle.  Check
        :attr:`Kinematics.feasible` before trusting the angles: an incompatible
        design is still returned, with clipped arccosine arguments, so that the
        optimizer can read ``delta_c1`` / ``delta_c2`` and be steered by them.
    """
    angles = crank_angles(samples) if theta_1 is None else np.asarray(theta_1, float)

    a, c, e = design.a, design.c, design.e
    q_1, q_2, I = design.q_1, design.q_2, design.I
    b, theta_b = design.b, design.theta_b
    theta_r = design.theta_r_rad

    # (1) the 1:2 gear pair ties the eccentric shaft to the crankshaft.
    theta_2 = -2.0 * angles + design.theta_f_rad

    # (3a) closure of the four-bar R1-Q-A-D-R2 projected on the axes.
    A_proj = q_1 * np.sin(angles) - q_2 * np.sin(theta_2) + I * np.cos(theta_r)
    B_proj = -q_1 * np.cos(angles) + q_2 * np.cos(theta_2) + I * np.sin(theta_r)

    # (4)/(4a) transmission angle between the swing rod and the trigonal link.
    cos_T = (A_proj**2 + B_proj**2 - a**2 - c**2) / (2.0 * a * c)
    delta_c1 = float(np.max(np.abs(cos_T)))
    T = np.arccos(np.clip(cos_T, -1.0, 1.0))

    # (2a)/(5) orientation of the trigonal link, then of the swing rod.
    q = np.arctan2(a * np.sin(T), a * np.cos(T) + c)
    theta_T = np.arctan2(B_proj, A_proj) - q
    theta_a = theta_T + T

    # (6)/(6a) piston-rod angle from the horizontal projection of chain (3).
    cos_e = (
        q_1 * np.sin(angles) - a * np.cos(theta_a) - b * np.cos(theta_b + theta_T) + design.x_1
    ) / e
    delta_c2 = float(np.max(np.abs(cos_e)))
    theta_e = np.arccos(np.clip(cos_e, -1.0, 1.0))

    # (7) piston crown height.  The report's equation (7) omits the constant
    # piston length p; it is restored here so that ``lam`` really is the height
    # of H.  Being a constant it cancels out of every stroke and volume anyway.
    lam = (
        q_1 * np.cos(angles)
        + a * np.sin(theta_a)
        + b * np.sin(theta_b + theta_T)
        + e * np.sin(theta_e)
        + spec.piston_length
    )

    feasible = bool(delta_c1 < 1.0 and delta_c2 < 1.0)

    # -- joint positions, rebuilt from the same chains ---------------------------
    zeros = np.zeros_like(angles)
    R1 = np.stack([zeros, zeros], axis=-1)
    R2 = np.stack(
        [np.full_like(angles, I * np.cos(theta_r)), np.full_like(angles, I * np.sin(theta_r))],
        axis=-1,
    )
    Q = np.stack([-q_1 * np.sin(angles), q_1 * np.cos(angles)], axis=-1)
    A_pt = Q + a * np.stack([np.cos(theta_a), np.sin(theta_a)], axis=-1)
    D_pt = A_pt + c * np.stack([np.cos(theta_T), np.sin(theta_T)], axis=-1)
    E_pt = A_pt + b * np.stack([np.cos(theta_b + theta_T), np.sin(theta_b + theta_T)], axis=-1)
    P_pt = E_pt + e * np.stack([np.cos(theta_e), np.sin(theta_e)], axis=-1)
    H_pt = np.stack([np.full_like(angles, design.x_1), lam], axis=-1)

    return Kinematics(
        design=design,
        spec=spec,
        theta_1=angles,
        theta_2=theta_2,
        theta_T=theta_T,
        theta_a=theta_a,
        theta_e=theta_e,
        lam=lam,
        transmission_angle=T,
        delta_c1=delta_c1,
        delta_c2=delta_c2,
        feasible=feasible,
        R1=R1,
        R2=R2,
        Q=Q,
        A=A_pt,
        D=D_pt,
        E=E_pt,
        P=P_pt,
        H=H_pt,
    )
