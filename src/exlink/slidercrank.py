"""The mechanism the EX-link replaces, modelled on identical terms.

Every conclusion in this package so far is about one linkage.  A reader is
entitled to ask two questions that one linkage cannot answer.

**Is the singularity finding general, or a quirk of this topology?**  The claim
is that a mechanism optimised quasi-statically drifts to its transmission-angle
singularity, because that is where the lever arm is longest -- and that the
same proximity is what amplifies the inertia loads, so the quasi-static optimum
is the worst possible starting point once the parts have mass.  If that is a
real mechanism-design principle it should appear in the simplest possible
linkage too, and here it does: a slider-crank's lever arm also peaks as its
connecting rod approaches alignment with the crank.

**Is the extra complexity worth it at all?**  The EX-link buys extended
expansion -- the gas expands through a larger volume ratio than it was
compressed through -- and pays for it with seven members instead of two, seven
journals instead of three, and a gear pair.  The geometric formulation could
not even pose that question: it has no way to price a member.  With range as
the objective it becomes a straight comparison, and the answer is not obvious
in advance.  Extended expansion is worth a few points of indicated efficiency;
four extra journals and a gear train are worth a few points of mechanical
efficiency in the other direction, plus their mass.

Both sides have to be optimised
-------------------------------
A comparison between an optimised EX-link and a slider-crank whose proportions
were written down from a textbook measures the optimization, not the topology,
and it flatters the side that was searched.  :func:`optimise_slidercrank`
removes that asymmetry by maximising the conventional engine's range over the
two degrees of freedom it has -- the rod length and the speed -- under the same
models.  It is worth 7.4 % to the baseline, and it changes the *sign* of the
study's headline comparison once the firing-frequency difference below is taken
out.

Modelling parity
----------------
The comparison is only worth making if both mechanisms are treated identically,
so the slider-crank goes through the *same* code wherever the code is not
topology-specific: the same material and safety factors, the same yield,
fatigue and buckling sizing in :mod:`exlink.sizing`, the same Coulomb friction
coefficients from :mod:`exlink.friction`, the same crankcase, bearing, shaft and
flywheel models from :mod:`exlink.mass_budget`, and the same vehicle and
burn-and-coast strategy from :mod:`exlink.vehicle`.  What differs is only what
must: the kinematics, the equilibrium system, and the cycle.

One asymmetry is *not* a modelling choice and has to be stated plainly.  In
this model the EX-link completes all four strokes in one crankshaft revolution
-- that is exactly what :func:`exlink.cycle.find_phases` requires of it -- while
a conventional slider-crank needs two.  The EX-link therefore fires twice as
often for the same displacement and speed.  That is a power-density advantage,
not an efficiency one, and since range is fuel per unit of *work* it does not
flatter the EX-link's range directly; it shows up only through the operating
point the vehicle can use.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .constants import DEFAULT_SPEC, DEFAULT_TARGETS, DesignTargets, EngineSpec
from .derivatives import ramp_derivative, spectral_derivative
from .friction import JOURNAL_FRICTION, PISTON_FRICTION, RING_TENSION
from .mass_budget import (
    CASE_CLEARANCE,
    FLYWHEEL_WEB_FACTOR,
    SPEED_FLUCTUATION,
    MassBudget,
    bearing_mass,
    crankcase_mass,
    cylinder_mass,
    flywheel_requirement,
    shaft_diameter,
)
from .materials import (
    DEFAULT_MATERIAL,
    DEFAULT_SAFETY,
    FloatArray,
    Material,
    SafetyFactors,
)
from .sizing import (
    MEMBER_FIXITY,
    internal_loads,
    piston_mass_from_pressure,
    size_from_arrays,
)
from .vehicle import Vehicle, best_strategy, brake_efficiency

STATIONS = 9
"""Sections evaluated along each member, matching the EX-link sizing."""


@dataclass(frozen=True)
class SliderCrank:
    """A conventional slider-crank, sized by its two lengths."""

    crank: float
    """Crank radius ``r`` [mm]."""

    rod: float
    """Connecting-rod length ``l`` [mm]."""

    @property
    def stroke(self) -> float:
        """Piston travel ``2r`` [mm]."""
        return 2.0 * self.crank

    @property
    def obliquity(self) -> float:
        """``r / l``, the rod obliquity ratio.

        The slider-crank's own conditioning number.  As it approaches 1 the rod
        aligns with the crank at some point in the revolution, the mechanism
        loses its ability to turn the gas force into torque there, and the
        joint forces diverge -- the same phenomenon the EX-link's ``W``
        constraint guards against.  Practical engines sit at 0.25 to 0.33.
        """
        return self.crank / self.rod

    @classmethod
    def for_compression_ratio(
        cls,
        ratio: float,
        obliquity: float = 0.30,
        spec: EngineSpec = DEFAULT_SPEC,
    ) -> SliderCrank:
        """Build the slider-crank that realises a compression ratio.

        With a fixed clearance volume, the compression ratio fixes the swept
        volume and hence the stroke; the rod length then follows from the
        chosen obliquity.

        Args:
            ratio: Required compression ratio ``epsilon``.
            obliquity: ``r / l``.
            spec: Fixed engine data.

        Returns:
            The mechanism.
        """
        swept = spec.dead_volume * (ratio - 1.0)
        crank = 0.5 * swept / spec.piston_area
        return cls(crank=crank, rod=crank / obliquity)


def kinematics(mechanism: SliderCrank, samples: int = 720) -> dict[str, FloatArray]:
    """Piston motion and rod angle over one crankshaft revolution.

    .. math:: \\lambda(\\theta) = r \\cos\\theta + \\sqrt{l^2 - r^2 \\sin^2\\theta}

    Args:
        mechanism: The slider-crank.
        samples: Crank angles per revolution.

    Returns:
        ``theta``, ``lam``, ``rod_angle``, and the two joint trajectories.
    """
    theta = np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False)
    r, length = mechanism.crank, mechanism.rod
    sin_phi = r * np.sin(theta) / length
    phi = np.arcsin(np.clip(sin_phi, -1.0, 1.0))
    lam = r * np.cos(theta) + length * np.cos(phi)
    crank_pin = np.stack([r * np.sin(theta), r * np.cos(theta)], axis=-1)
    wrist = np.stack([np.zeros_like(theta), lam], axis=-1)
    return {
        "theta": theta,
        "lam": lam,
        "rod_angle": phi,
        "crank_pin": crank_pin,
        "wrist": wrist,
    }


def otto_cycle(
    lam: FloatArray, spec: EngineSpec = DEFAULT_SPEC
) -> dict[str, FloatArray | float]:
    """The equivalent-Otto cycle on a slider-crank motion.

    Compression and expansion through the *same* volume ratio -- which is the
    whole difference from the Atkinson cycle, and the reason extended expansion
    is worth building a linkage for.  The four strokes occupy two crankshaft
    revolutions, so the pressure is assigned over a doubled angle base and the
    work per revolution is half the work per cycle.

    Args:
        lam: Piston height over one crank revolution [mm].
        spec: Fixed engine data.

    Returns:
        ``pressure``, ``gauge_pressure``, ``piston_force`` over two
        revolutions, plus ``p_compression_end``, ``p_combustion``,
        ``compression_ratio`` and ``indicated_work`` per cycle [N.mm].
    """
    top = float(np.max(lam))
    bottom = float(np.min(lam))
    swept = (top - bottom) * spec.piston_area
    volume_top = spec.dead_volume
    v_bottom = volume_top + swept
    ratio = v_bottom / volume_top

    gamma = spec.heat_capacity_ratio
    p_intake = spec.p_intake
    p_two = p_intake * ratio**gamma
    p_three = spec.explosion_ratio * p_two

    # The cycle occupies two crankshaft revolutions, so the motion is repeated
    # and the volume follows it.
    doubled = np.concatenate([lam, lam])
    volume = volume_top + (top - doubled) * spec.piston_area

    n = lam.size
    bottom_index = int(np.argmin(lam))
    index = np.arange(n)
    descending = index < bottom_index

    pressure = np.empty_like(volume)
    # First revolution: intake down at plenum pressure, then adiabatic
    # compression back up from bottom dead centre.
    pressure[:n] = np.where(descending, p_intake, p_intake * (v_bottom / volume[:n]) ** gamma)
    # Second: adiabatic expansion down from the post-combustion state, then
    # exhaust up at plenum pressure.  Expansion and compression share the same
    # volume ratio -- that is exactly what extended expansion changes.
    pressure[n:] = np.where(descending, p_three * (volume_top / volume[n:]) ** gamma, p_intake)

    gauge = pressure - p_intake
    work = float(np.trapezoid(pressure, volume))
    return {
        "pressure": pressure,
        "gauge_pressure": gauge,
        "piston_force": gauge * spec.piston_area,
        "volume": volume,
        "p_compression_end": p_two,
        "p_combustion": p_three,
        "compression_ratio": ratio,
        "indicated_work": work,
    }


#: The two sized members, and their Euler end-fixity kind.
SC_MEMBERS: tuple[tuple[str, str], ...] = (("crank", "cantilever"), ("rod", "link"))


@dataclass(frozen=True)
class SliderCrankResult:
    """A sized, loaded slider-crank at one operating point."""

    mechanism: SliderCrank
    speed: float
    """Crankshaft speed [rad/s]."""

    diameters: dict[str, float]
    member_mass: dict[str, float]
    """[tonne]"""

    piston_mass: float
    """[tonne]"""

    torque: FloatArray
    """Output torque over two crankshaft revolutions [N.mm]."""

    reaction: dict[str, FloatArray]
    liner_force: FloatArray
    lam: FloatArray
    """Piston height over the two revolutions [mm]."""

    gas_force: FloatArray
    """Gas force on the crown over the two revolutions [N]."""

    indicated_work: float
    """Work per cycle, from the p-V loop [N.mm]."""

    heat_release: float
    """Heat added per cycle [N.mm]."""

    peak_pressure: float
    """Peak gauge pressure [MPa]."""

    conditioning: float
    """Worst condition number of the 9x9 equilibrium matrix."""

    converged: bool

    @property
    def mean_torque(self) -> float:
        """Mean output torque over the full cycle [N.mm]."""
        return float(np.mean(self.torque))

    @property
    def peak_bearing_load(self) -> float:
        """Largest main-journal reaction [N]."""
        return float(np.max(np.linalg.norm(self.reaction["O"], axis=1)))

    @property
    def height(self) -> float:
        """Envelope along the stroke [mm]."""
        return float(np.max(self.lam) + self.mechanism.crank)

    @property
    def width(self) -> float:
        """Envelope across the stroke [mm]."""
        return 2.0 * self.mechanism.crank


def _solve_equilibrium(
    mechanism: SliderCrank,
    motion: dict[str, FloatArray],
    gas_force: FloatArray,
    masses: dict[str, float],
    piston: float,
    speed: float,
    spec: EngineSpec,
) -> tuple[dict[str, FloatArray], FloatArray, FloatArray, float, dict[str, FloatArray]]:
    """Simultaneous equilibrium of the three bodies, with inertia.

    Nine unknowns -- three joint force pairs, the liner normal force, the liner
    reaction moment and the output torque -- and nine equations.  Exactly the
    construction :mod:`exlink.dynamics` uses on the EX-link, at a third of the
    size, and for the same reason: with mass in the rod it stops being a
    two-force member and sequential elimination no longer closes.
    """
    theta = motion["theta"]
    n = theta.size
    doubled = np.concatenate([theta, theta + 2.0 * np.pi])
    del doubled

    crank_pin = np.concatenate([motion["crank_pin"], motion["crank_pin"]])
    wrist = np.concatenate([motion["wrist"], motion["wrist"]])
    origin = np.zeros_like(wrist)

    def second(points: FloatArray) -> FloatArray:
        return np.stack(
            [spectral_derivative(points[:, 0], 2), spectral_derivative(points[:, 1], 2)],
            axis=-1,
        )

    scale = speed**2
    a_pin = scale * second(crank_pin)
    a_wrist = scale * second(wrist)
    com_crank = 0.5 * crank_pin
    com_rod = 0.5 * (crank_pin + wrist)
    com_piston = wrist + np.array([0.0, 0.5 * spec.piston_length])
    a_crank = scale * second(com_crank)
    a_rod = scale * second(com_rod)
    a_piston = scale * second(com_piston)

    rod_angle = np.concatenate([motion["rod_angle"], motion["rod_angle"]])
    alpha_rod = scale * ramp_derivative(rod_angle, 2)

    m_crank, m_rod = masses["crank"], masses["rod"]
    length_rod = mechanism.rod
    i_rod = m_rod * length_rod**2 / 12.0

    total = 2 * n
    matrix = np.zeros((total, 9, 9))
    rhs = np.zeros((total, 9))
    # Unknown order: R_Ox R_Oy R_Ax R_Ay R_Bx R_By N M_liner M_torque
    o_x, o_y, a_x, a_y, b_x, b_y, i_n, i_m, i_t = range(9)

    def cross(u: FloatArray, v_index: tuple[int, int], sign: float, row: int) -> None:
        matrix[:, row, v_index[0]] += -sign * u[:, 1]
        matrix[:, row, v_index[1]] += sign * u[:, 0]

    # -- crank: R_O - R_A = m a, moments about its own centre of mass ------------
    matrix[:, 0, o_x] += 1.0
    matrix[:, 0, a_x] += -1.0
    matrix[:, 1, o_y] += 1.0
    matrix[:, 1, a_y] += -1.0
    cross(origin - com_crank, (o_x, o_y), +1.0, 2)
    cross(crank_pin - com_crank, (a_x, a_y), -1.0, 2)
    # Signed so that a positive unknown is torque delivered *out* of the shaft.
    # This crank turns the opposite way to the EX-link's, so the sign is the
    # opposite of the one in :mod:`exlink.dynamics`; the check that fixes it is
    # that the torque integral must equal the p-V loop area, not a convention.
    matrix[:, 2, i_t] += 1.0
    rhs[:, 0] += m_crank * a_crank[:, 0]
    rhs[:, 1] += m_crank * a_crank[:, 1]

    # -- rod: R_A - R_B = m a -----------------------------------------------------
    matrix[:, 3, a_x] += 1.0
    matrix[:, 3, b_x] += -1.0
    matrix[:, 4, a_y] += 1.0
    matrix[:, 4, b_y] += -1.0
    cross(crank_pin - com_rod, (a_x, a_y), +1.0, 5)
    cross(wrist - com_rod, (b_x, b_y), -1.0, 5)
    rhs[:, 3] += m_rod * a_rod[:, 0]
    rhs[:, 4] += m_rod * a_rod[:, 1]
    rhs[:, 5] += i_rod * alpha_rod

    # -- piston: R_B + N x + gas = m a -------------------------------------------
    matrix[:, 6, b_x] += 1.0
    matrix[:, 6, i_n] += 1.0
    matrix[:, 7, b_y] += 1.0
    cross(wrist - com_piston, (b_x, b_y), +1.0, 8)
    matrix[:, 8, i_m] += 1.0
    rhs[:, 6] += piston * a_piston[:, 0]
    # The gas resultant is (0, -F) on the crown, so it moves to the right-hand
    # side as +F, exactly as in the EX-link model.
    rhs[:, 7] += piston * a_piston[:, 1] + gas_force

    solution: FloatArray = np.linalg.solve(matrix, rhs[..., None])[..., 0].astype(
        np.float64, copy=False
    )
    conditioning = float(np.max(np.linalg.cond(matrix)))
    reaction = {
        "O": solution[:, 0:2],
        "A": solution[:, 2:4],
        "B": solution[:, 4:6],
    }
    accelerations = {
        "pin": a_pin,
        "wrist": a_wrist,
        "origin": np.zeros_like(a_pin),
    }
    return reaction, solution[:, i_n], solution[:, i_t], conditioning, accelerations


def solve(
    mechanism: SliderCrank,
    speed_rpm: float,
    samples: int = 360,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
    spec: EngineSpec = DEFAULT_SPEC,
    tolerance: float = 1.0e-6,
    max_iterations: int = 200,
) -> SliderCrankResult:
    """Solve the slider-crank's sizing/dynamics fixed point.

    The same coupling as the EX-link's and solved the same way -- sections set
    masses, masses set inertia loads, loads set sections -- so that the two
    mechanisms are compared under identical treatment rather than one being
    given a converged structure and the other an assumed one.

    Args:
        mechanism: The slider-crank.
        speed_rpm: Crankshaft speed [rev/min].
        samples: Crank angles per revolution.
        material: The material.
        safety: The design factors.
        spec: Fixed engine data.
        tolerance: Convergence tolerance on the diameter change [mm].
        max_iterations: Sweep limit.

    Returns:
        The converged result.
    """
    motion = kinematics(mechanism, samples)
    cycle = otto_cycle(motion["lam"], spec)
    gas = np.asarray(cycle["piston_force"], dtype=float)
    peak_pressure = float(np.max(np.asarray(cycle["gauge_pressure"])))
    speed = speed_rpm * 2.0 * math.pi / 60.0

    crown, piston = piston_mass_from_pressure(peak_pressure, material, safety, spec)
    del crown

    lengths = np.array([mechanism.crank, mechanism.rod])
    diameters = np.array([8.0, 8.0])
    doubled_lam = np.concatenate([motion["lam"], motion["lam"]])

    reaction: dict[str, FloatArray] = {}
    liner = np.zeros(2 * samples)
    torque = np.zeros(2 * samples)
    conditioning = 0.0
    converged = False

    for _ in range(max_iterations):
        area = math.pi * diameters**2 / 4.0
        masses = {
            "crank": float(material.density * area[0] * lengths[0]),
            "rod": float(material.density * area[1] * lengths[1]),
        }
        reaction, liner, torque, conditioning, acceleration = _solve_equilibrium(
            mechanism, motion, gas, masses, piston, speed, spec
        )

        crank_pin = np.concatenate([motion["crank_pin"], motion["crank_pin"]])
        wrist = np.concatenate([motion["wrist"], motion["wrist"]])
        origin = np.zeros_like(wrist)

        axial_crank, bending_crank = internal_loads(
            crank_pin,
            origin,
            -reaction["A"],
            masses["crank"],
            acceleration["pin"],
            acceleration["origin"],
            stations=STATIONS,
        )
        axial_rod, bending_rod = internal_loads(
            crank_pin,
            wrist,
            reaction["A"],
            masses["rod"],
            acceleration["pin"],
            acceleration["wrist"],
            stations=STATIONS,
        )
        sized = size_from_arrays(
            np.stack([axial_crank, axial_rod]),
            np.stack([bending_crank, bending_rod]),
            lengths,
            material,
            safety,
            fixity=np.array([MEMBER_FIXITY[0], MEMBER_FIXITY[1]]),
            names=("crank", "rod"),
        )
        updated = np.array([sized["crank"].diameter, sized["rod"].diameter])
        residual = float(np.max(np.abs(updated - diameters)))
        diameters = updated
        if residual <= tolerance:
            converged = True
            break

    area = math.pi * diameters**2 / 4.0
    member_mass = {
        "crank": float(material.density * area[0] * lengths[0]),
        "rod": float(material.density * area[1] * lengths[1]),
    }
    quantity = (
        spec.dead_volume
        * (float(cycle["p_combustion"]) - float(cycle["p_compression_end"]))
        / (spec.heat_capacity_ratio - 1.0)
    )
    return SliderCrankResult(
        mechanism=mechanism,
        speed=speed,
        diameters={"crank": float(diameters[0]), "rod": float(diameters[1])},
        member_mass=member_mass,
        piston_mass=piston,
        torque=torque,
        reaction=reaction,
        liner_force=liner,
        lam=doubled_lam,
        gas_force=gas,
        indicated_work=float(cycle["indicated_work"]),
        heat_release=float(quantity),
        peak_pressure=peak_pressure,
        conditioning=conditioning,
        converged=converged,
    )


def friction_work(
    result: SliderCrankResult,
    journal_friction: float = JOURNAL_FRICTION,
    piston_friction: float = PISTON_FRICTION,
    ring_tension: float = RING_TENSION,
) -> float:
    """Mechanical loss over one full cycle [N.mm].

    Three journals instead of the EX-link's seven, and no gear mesh.  That is
    the slider-crank's structural advantage, and it is exactly what has to be
    weighed against extended expansion.

    Args:
        result: A solved slider-crank.
        journal_friction: Coulomb coefficient at the journals.
        piston_friction: Coulomb coefficient at the liner.
        ring_tension: Ring-pack radial load [N].

    Returns:
        Friction work per cycle [N.mm].
    """
    total = result.torque.size
    step = 4.0 * math.pi / total
    motion = kinematics(result.mechanism, total // 2)
    rod_angle = np.concatenate([motion["rod_angle"], motion["rod_angle"]])
    crank_angle = np.linspace(0.0, 4.0 * math.pi, total, endpoint=False)

    # Relative rotation: the main journal turns through the full crank angle,
    # the crank pin through crank minus rod angle, the wrist pin through the
    # rod angle alone.
    rates = {
        "O": np.ones(total),
        "A": np.abs(1.0 - ramp_derivative(rod_angle, 1)),
        "B": np.abs(ramp_derivative(rod_angle, 1)),
    }
    del crank_angle
    radii = {
        "O": 0.5 * result.diameters["crank"],
        "A": 0.5 * max(result.diameters["crank"], result.diameters["rod"]),
        "B": 0.5 * result.diameters["rod"],
    }
    bearings = sum(
        journal_friction
        * radii[joint]
        * float(np.sum(np.linalg.norm(result.reaction[joint], axis=1) * rates[joint]))
        * step
        for joint in ("O", "A", "B")
    )
    slide = np.abs(spectral_derivative(result.lam, 1))
    normal = np.abs(result.liner_force) + ring_tension
    piston = piston_friction * float(np.sum(normal * slide)) * step
    return float(bearings + piston)


def mass_budget(
    result: SliderCrankResult,
    fluctuation: float = SPEED_FLUCTUATION,
    material: Material = DEFAULT_MATERIAL,
    safety: SafetyFactors = DEFAULT_SAFETY,
    spec: EngineSpec = DEFAULT_SPEC,
) -> MassBudget:
    """The slider-crank's engine mass, on the EX-link's terms.

    Every item is computed by the same function :mod:`exlink.mass_budget` uses
    for the EX-link -- shafts, bearings, crankcase, cylinder head, flywheel --
    so that the comparison is a comparison of mechanisms and not of modelling
    conventions.  Two entries differ, and only because the mechanism does:
    there is no gear pair, and the flywheel has to carry the engine through
    two revolutions on one firing rather than one.

    Args:
        result: A solved slider-crank.
        fluctuation: Allowed cyclic speed fluctuation.
        material: The material.
        safety: The design factors.
        spec: Fixed engine data.

    Returns:
        The itemised budget.
    """
    peak_reaction = result.peak_bearing_load
    peak_torque = float(np.max(np.abs(result.torque)))

    journal = shaft_diameter(
        peak_reaction, peak_torque, result.mechanism.crank, material, safety
    )
    case_depth = spec.bore + 2.0 * (0.55 * journal + 3.0) + 4.0 * CASE_CLEARANCE
    shaft_length = case_depth + 2.0 * (0.55 * journal + 3.0)
    shaft = material.density * math.pi * journal**2 / 4.0 * shaft_length
    bearings = bearing_mass(journal, 2, material.density)

    # The turning-moment diagram, over the full two-revolution cycle.  As in
    # the EX-link budget this uses the gas torque, not the total: the inertia
    # part is energy traded with the mechanism's own masses and is accounted
    # for by the inherent rotating inertia below.
    torque_gas = -result.gas_force * spectral_derivative(result.lam, 1)
    required, _swing = flywheel_requirement(
        torque_gas, result.speed, fluctuation, span=4.0 * math.pi
    )
    crank_inertia = result.member_mass["crank"] * result.mechanism.crank**2 / 3.0
    deficit = max(required - crank_inertia, 0.0)
    radius = min(max(0.45 * result.width, 30.0), 150.0)
    flywheel = FLYWHEEL_WEB_FACTOR * deficit / radius**2

    case, _wall = crankcase_mass(result.height, result.width, case_depth, peak_reaction)
    cylinder = cylinder_mass(
        result.mechanism.stroke, result.peak_pressure, spec, material, safety
    )

    items = {
        "linkage": float(sum(result.member_mass.values())),
        "piston": float(result.piston_mass),
        "gears": 0.0,
        "shafts": float(shaft),
        "bearings": float(bearings),
        "crankcase": float(case),
        "cylinder_head": float(cylinder),
        "flywheel": float(flywheel),
    }
    return MassBudget(
        items=items,
        gears=None,
        shaft_diameter=journal,
        flywheel_inertia=deficit,
        flywheel_radius=radius,
        required_inertia=required,
        inherent_inertia=crank_inertia,
    )


@dataclass(frozen=True)
class Comparison:
    """One mechanism's numbers, for the head-to-head table."""

    name: str
    speed_rpm: float
    indicated_efficiency: float
    mechanical_efficiency: float
    brake_efficiency: float
    engine_mass: float
    """[kg]"""

    brake_power: float
    """[W]"""

    km_per_litre: float
    joints: int
    members: int
    feasible: bool


def evaluate_slidercrank(
    mechanism: SliderCrank,
    speed_rpm: float,
    vehicle: Vehicle | None = None,
    samples: int = 360,
    **kwargs: object,
) -> Comparison:
    """Run the slider-crank all the way to kilometres per litre.

    Args:
        mechanism: The slider-crank.
        speed_rpm: Crankshaft speed [rev/min].
        vehicle: The car; a default Prototype-class entry if omitted.
        samples: Crank angles per revolution.
        **kwargs: Forwarded to :func:`solve`.

    Returns:
        The comparison row.
    """
    car = vehicle if vehicle is not None else Vehicle()
    result = solve(mechanism, speed_rpm, samples=samples, **kwargs)  # type: ignore[arg-type]
    loss = friction_work(result)
    brake = result.indicated_work - loss
    budget = mass_budget(result)
    efficiency = brake_efficiency(brake, result.heat_release)
    # Two revolutions per cycle, so half a cycle's work per revolution.
    power = brake / 1000.0 * (speed_rpm / 60.0) / 2.0
    outcome = best_strategy(car, budget.total_kg, power, efficiency)
    return Comparison(
        name="slider-crank (Otto)",
        speed_rpm=speed_rpm,
        indicated_efficiency=result.indicated_work / result.heat_release,
        mechanical_efficiency=brake / result.indicated_work,
        brake_efficiency=efficiency,
        engine_mass=budget.total_kg,
        brake_power=power,
        km_per_litre=outcome.km_per_litre,
        joints=3,
        members=2,
        feasible=result.converged and outcome.feasible and brake > 0.0,
    )


OBLIQUITY_BOUNDS = (0.12, 0.45)
"""Search bounds on ``r / l`` for the optimised baseline.

The upper end is where the rod is short enough that it approaches alignment
with the crank -- the slider-crank's own transmission-angle singularity, and
the same pathology :mod:`exlink.metrics` guards against with ``W``.  The lower
end is where the rod is nearly four times the stroke: :func:`solve` still sizes
it, but the engine is by then taller than the car, and the model prices rod
length only through mass, not through the packaging that would really stop it.
Practical engines occupy 0.25 to 0.33, comfortably inside.
"""

SPEED_BOUNDS = (800.0, 3200.0)
"""Search bounds on crankshaft speed [rev/min] for the optimised baseline."""


@dataclass(frozen=True)
class OptimisedSliderCrank:
    """The best conventional engine this model can build, and how it was found."""

    mechanism: SliderCrank
    speed_rpm: float
    comparison: Comparison
    evaluations: int
    starts: int


def optimise_slidercrank(
    ratio: float = 16.0,
    vehicle: Vehicle | None = None,
    samples: int = 360,
    spec: EngineSpec = DEFAULT_SPEC,
    obliquity_bounds: tuple[float, float] = OBLIQUITY_BOUNDS,
    speed_bounds: tuple[float, float] = SPEED_BOUNDS,
    starts: int = 4,
) -> OptimisedSliderCrank:
    """Maximise a conventional engine's range over the variables it has.

    Why this exists
    ---------------
    Comparing an *optimised* EX-link against a slider-crank whose proportions
    were written down by hand does not measure the topology; it measures the
    optimization.  The honest comparison is optimum against optimum, and this
    is the other optimum.

    What is free and what is not
    ----------------------------
    The compression ratio is held at the value the EX-link is *required* to
    reach, so both engines trap the same charge in the same clearance volume
    and burn the same fuel per cycle.  It is not made a design variable,
    because this model has no knock limit: left free it would rise without
    bound and the comparison would measure a missing sub-model rather than a
    mechanism.  With the ratio fixed and the clearance volume fixed, the stroke
    follows, and the slider-crank has exactly two degrees of freedom left --
    the rod length, as the obliquity ``r / l``, and the speed it is run at.

    Why a derivative-free method is admissible here, having been rejected there
    ---------------------------------------------------------------------------
    The EX-link's feasible set is a manifold of measure zero, which is what
    excludes sampling methods from that problem.  This one is a box in two
    variables: every point in it is feasible to evaluate, so Nelder-Mead is a
    reasonable choice and no gradient has to be derived for a baseline.  The
    restarts are there because the objective is not concave -- range rises with
    speed through the vehicle's operating point and falls with it through
    friction and the inertia loads that set the mass.

    Args:
        ratio: Compression ratio to hold, matching the EX-link's requirement.
        vehicle: The car; a default Prototype-class entry if omitted.
        samples: Crank angles per revolution.
        spec: Fixed engine data.
        obliquity_bounds: Search interval for ``r / l``.
        speed_bounds: Search interval for the speed [rev/min].
        starts: Nelder-Mead restarts, spread over the box.

    Returns:
        The best mechanism, its speed, and its scored row.

    Raises:
        RuntimeError: If no point in the box produced a feasible engine.
    """
    from scipy.optimize import minimize

    car = vehicle if vehicle is not None else Vehicle()
    calls = 0
    best: tuple[float, SliderCrank, float, Comparison] | None = None

    def score(vector: FloatArray) -> float:
        nonlocal calls, best
        calls += 1
        obliquity = float(np.clip(vector[0], *obliquity_bounds))
        speed = float(np.clip(vector[1], *speed_bounds))
        mechanism = SliderCrank.for_compression_ratio(ratio, obliquity, spec=spec)
        try:
            row = evaluate_slidercrank(mechanism, speed, vehicle=car, samples=samples)
        except (ValueError, FloatingPointError):
            return 0.0
        if not row.feasible:
            return 0.0
        if best is None or row.km_per_litre > best[0]:
            best = (row.km_per_litre, mechanism, speed, row)
        return row.km_per_litre

    # Nelder-Mead has no bounds of its own in every SciPy it runs on here, so
    # the objective clips and the simplex is started well inside.
    low, high = obliquity_bounds
    slow, fast = speed_bounds
    for index in range(max(int(starts), 1)):
        fraction = (index + 0.5) / max(int(starts), 1)
        guess = np.array([low + fraction * (high - low), slow + fraction * (fast - slow)])
        minimize(
            lambda vector: -score(vector),
            guess,
            method="Nelder-Mead",
            options={"xatol": 1.0e-4, "fatol": 1.0e-3, "maxiter": 200},
        )

    if best is None:
        msg = "no feasible slider-crank in the search box"
        raise RuntimeError(msg)
    _range, mechanism, speed, row = best
    return OptimisedSliderCrank(
        mechanism=mechanism,
        speed_rpm=speed,
        comparison=row,
        evaluations=calls,
        starts=int(starts),
    )


def optimise_slidercrank_constrained(
    ratio: float = 16.0,
    vehicle: Vehicle | None = None,
    samples: int = 360,
    spec: EngineSpec = DEFAULT_SPEC,
    obliquity_bounds: tuple[float, float] = OBLIQUITY_BOUNDS,
    speed_bounds: tuple[float, float] = SPEED_BOUNDS,
    max_iterations: int = 80,
) -> OptimisedSliderCrank:
    """The same baseline, optimised the way the EX-link now is.

    Why this exists alongside :func:`optimise_slidercrank`
    ------------------------------------------------------
    A comparison is only worth making if both sides get the same treatment, and
    that has to include the *method*, not only the models.  The EX-link's best
    design comes from an SQP that holds every constraint at every step (§3.10);
    scoring the baseline with a derivative-free search that simply rejects
    infeasible points would leave the comparison measuring the optimizer again,
    which is the error :func:`optimise_slidercrank` was written to remove one
    level up.

    So this poses the baseline the same way: maximise range subject to the
    constraints as constraints, by SLSQP.  The slider-crank has two of them --
    the engine must produce net work, and the vehicle must meet the average
    speed rule -- against the EX-link's fourteen, because a slider-crank has
    one dead centre, no transmission angle to lose and no gear train.

    What is *not* matched, and cannot be here
    ------------------------------------------
    The EX-link's final formulation also constrains a system reliability index
    over its eleven dimensions (§3.10).  There is no equivalent for the
    slider-crank in this package: :mod:`exlink.robustness` builds its
    covariance and its constraint Jacobians for the EX-link's design vector
    specifically.  Constructing the analogue for a two-variable mechanism is
    tractable but is not the same model, and asserting a reliability comparison
    across two different uncertainty models would be worse than declining one.
    The comparison in §6.3 is therefore between two *nominally* optimised
    engines, and says so.

    Args:
        ratio: Compression ratio to hold, matching the EX-link's requirement.
        vehicle: The car; a default Prototype-class entry if omitted.
        samples: Crank angles per revolution.
        spec: Fixed engine data.
        obliquity_bounds: Search interval for ``r / l``.
        speed_bounds: Search interval for the speed [rev/min].
        max_iterations: SLSQP iteration budget.

    Returns:
        The best mechanism, its speed, and its scored row.

    Raises:
        RuntimeError: If no point in the box produced a feasible engine.
    """
    from scipy.optimize import minimize

    car = vehicle if vehicle is not None else Vehicle()
    calls = 0
    best: tuple[float, SliderCrank, float, Comparison] | None = None

    cache: dict[bytes, tuple[SliderCrank, float, Comparison, float] | None] = {}

    def state(vector: FloatArray) -> tuple[SliderCrank, float, Comparison, float] | None:
        """Mechanism, speed, scored row and average road speed at one point."""
        key = np.ascontiguousarray(vector, dtype=float).tobytes()
        if key in cache:
            return cache[key]
        obliquity = float(np.clip(vector[0], *obliquity_bounds))
        speed = float(np.clip(vector[1], *speed_bounds))
        mechanism = SliderCrank.for_compression_ratio(ratio, obliquity, spec=spec)
        try:
            row = evaluate_slidercrank(mechanism, speed, vehicle=car, samples=samples)
            result = solve(mechanism, speed, samples=samples)
            loss = friction_work(result)
            brake = result.indicated_work - loss
            budget = mass_budget(result)
            outcome = best_strategy(
                car,
                budget.total_kg,
                brake / 1000.0 * (speed / 60.0) / 2.0,
                brake_efficiency(brake, result.heat_release),
            )
            answer = (mechanism, speed, row, outcome.average_speed)
        except (ValueError, FloatingPointError):
            answer = None
        cache.clear()
        cache[key] = answer
        return answer

    def objective(vector: FloatArray) -> float:
        nonlocal calls, best
        calls += 1
        found = state(vector)
        if found is None:
            return 0.0
        mechanism, speed, row, _average = found
        if row.feasible and (best is None or row.km_per_litre > best[0]):
            best = (row.km_per_litre, mechanism, speed, row)
        return -row.km_per_litre

    def constraint(vector: FloatArray) -> FloatArray:
        found = state(vector)
        if found is None:
            return np.full(2, -1.0e3)
        _mechanism, _speed, row, average = found
        # The two things ``Comparison.feasible`` folds into a bool: the engine
        # must make net work, and the car must hold the average-speed rule.
        return np.array([row.brake_power, average - car.minimum_average_speed])

    low, high = obliquity_bounds
    slow, fast = speed_bounds
    outcome = minimize(
        objective,
        np.array([0.5 * (low + high), 0.5 * (slow + fast)]),
        method="SLSQP",
        bounds=[(low, high), (slow, fast)],
        constraints=[{"type": "ineq", "fun": constraint}],
        options={"maxiter": int(max_iterations), "ftol": 1.0e-9},
    )
    del outcome
    if best is None:
        msg = "no feasible slider-crank in the search box"
        raise RuntimeError(msg)
    _range, mechanism, speed, row = best
    return OptimisedSliderCrank(
        mechanism=mechanism,
        speed_rpm=speed,
        comparison=row,
        evaluations=calls,
        starts=1,
    )


SLIDERCRANK_CONSTRAINTS: tuple[str, ...] = (
    "ratio_upper",
    "ratio_lower",
    "rod_angle",
    "side_load",
)
"""The requirements a slider-crank actually has, in the order reported.

Four, against the EX-link's seven, and the difference is the mechanism rather
than the treatment.  A slider-crank has one dead centre, so there is no
top-dead-centre gap; its coupler cannot approach a transmission-angle
singularity, so there is no compatibility condition; and its expansion and
compression strokes are the same length by construction, so ``STE = 74`` is not
a requirement it can be asked to meet -- that asymmetry is the whole reason the
EX-link exists.  What both must meet is the compression ratio, the rod-angle
limit and the side-load limit, and those are what is compared.
"""


@dataclass(frozen=True)
class SliderCrankReliability:
    """FORM reliability of a slider-crank, on the same terms as the EX-link's."""

    mechanism: SliderCrank
    grade: int
    value: FloatArray
    """Constraint values, negative meaning satisfied."""

    sigma: FloatArray
    """First-order standard deviation of each constraint."""

    correlation: FloatArray
    per_constraint: dict[str, float]
    system: float
    """Probability that *any* constraint fails, correlation kept."""

    independent_bound: float

    @property
    def beta(self) -> FloatArray:
        """Reliability index of each constraint."""
        safe = np.where(self.sigma > 0.0, self.sigma, np.inf)
        return -self.value / safe

    @property
    def system_beta(self) -> float:
        """System reliability index."""
        from scipy.stats import norm

        return float(norm.isf(min(max(self.system, 1e-16), 1.0 - 1e-16)))

    def binding(self) -> str:
        """The constraint contributing most of the failure probability."""
        return max(self.per_constraint, key=lambda name: self.per_constraint[name])


def _slidercrank_constraints(
    crank: float,
    rod: float,
    band: float,
    samples: int,
    targets: DesignTargets,
    spec: EngineSpec,
) -> FloatArray | None:
    """The four constraints, negative meaning satisfied.

    ``gamma`` is taken from the solved model rather than from a closed form.
    The obvious closed form -- ``tan`` of the maximum rod angle -- is wrong by
    80 %, because the peak side load and the peak gas force occur at *different*
    crank angles: the gas force peaks at top dead centre, where the rod is
    nearly upright and the side load is near zero.  A ratio of two maxima is
    not the maximum of a ratio.
    """
    mechanism = SliderCrank(crank=crank, rod=rod)
    try:
        result = solve(mechanism, 1.0, samples=samples, spec=spec)
    except (ValueError, FloatingPointError):
        return None
    if not result.converged:
        return None
    swept = spec.piston_area * mechanism.stroke
    ratio = (spec.dead_volume + swept) / spec.dead_volume
    obliquity = mechanism.obliquity
    if obliquity >= 1.0:
        return None
    rod_angle = math.degrees(math.asin(obliquity))
    gas = float(np.max(np.abs(np.asarray(result.gas_force))))
    liner = float(np.max(np.abs(np.asarray(result.liner_force))))
    side_load = liner / gas if gas > 0.0 else 1.0e3
    difference = ratio - targets.compression_ratio
    return np.array(
        [
            difference - band,
            -difference - band,
            rod_angle - targets.max_rod_angle,
            side_load - targets.max_side_load,
        ],
        dtype=float,
    )


def slidercrank_reliability(
    mechanism: SliderCrank,
    grade: int = 8,
    band: float = 0.05,
    samples: int = 360,
    targets: DesignTargets = DEFAULT_TARGETS,
    spec: EngineSpec = DEFAULT_SPEC,
    step: float = 1.0e-4,
) -> SliderCrankReliability | None:
    """Probability that a built slider-crank misses one of its requirements.

    The same method as :func:`exlink.robustness.failure_probability` -- ISO 286
    tolerances on the machined lengths, first-order propagation to each
    constraint, and the correlated multivariate-normal orthant for the system
    probability -- applied to the baseline so that §6.3 compares reliability
    with reliability rather than reliability with silence.

    What is the same, and what cannot be
    ------------------------------------
    Same: the IT grade, the ``sigma = half-width / 3`` convention, the
    definitions of the compression ratio, the rod angle and the side-load
    ratio, and the orthant integral.

    Not the same, and unavoidably so: the *dimensions* carrying the tolerance.
    The EX-link's covariance is over eleven lengths and two clocking angles;
    a slider-crank has two lengths.  A mechanism with fewer dimensions has
    fewer ways to be wrong, and a like-for-like reliability comparison between
    two mechanisms is therefore never wholly like-for-like.  What the
    comparison does establish is each engine's probability of missing *its own*
    requirements when machined to the same grade, which is the question a
    builder would ask.

    Gradients are central differences on two variables, not the analytic
    Jacobians §3.5 needs for the EX-link.  The reason §3.5 rejects differences
    is extremum switching, and it applies here too -- ``gamma`` is a ratio of
    two maxima -- so the step is checked against a tenth and ten times its
    value in the tests rather than assumed safe.

    Args:
        mechanism: The slider-crank.
        grade: ISO 286 IT grade of the machined lengths.
        band: Half-width on the compression-ratio requirement.
        samples: Crank angles per revolution.
        targets: Constraint right-hand sides.
        spec: Fixed engine data.
        step: Relative step for the central differences.

    Returns:
        The reliability, or ``None`` if the mechanism cannot be solved.
    """
    from scipy.stats import multivariate_normal, norm

    from .robustness import IT_FACTORS, SIGMA_PER_HALF_WIDTH, tolerance_unit

    lengths = np.array([mechanism.crank, mechanism.rod], dtype=float)
    nominal = _slidercrank_constraints(
        float(lengths[0]), float(lengths[1]), band, samples, targets, spec
    )
    if nominal is None:
        return None

    jacobian = np.zeros((nominal.size, lengths.size))
    for index, value in enumerate(lengths):
        delta = step * max(abs(float(value)), 1.0)
        rows = []
        for sign in (+1.0, -1.0):
            moved = lengths.copy()
            moved[index] += sign * delta
            row = _slidercrank_constraints(
                float(moved[0]), float(moved[1]), band, samples, targets, spec
            )
            if row is None:
                return None
            rows.append(row)
        jacobian[:, index] = (rows[0] - rows[1]) / (2.0 * delta)

    half = IT_FACTORS[grade] * np.array([tolerance_unit(float(value)) for value in lengths])
    sigma_matrix = np.diag((half / SIGMA_PER_HALF_WIDTH) ** 2)
    covariances = jacobian @ sigma_matrix @ jacobian.T
    sigma = np.sqrt(np.clip(np.diag(covariances), 0.0, None))
    outer = np.outer(sigma, sigma)
    correlation = np.divide(covariances, outer, out=np.eye(sigma.size), where=outer > 0.0)

    safe_sigma = np.where(sigma > 0.0, sigma, np.inf)
    beta = -nominal / safe_sigma
    per_constraint = {
        name: float(norm.sf(value))
        for name, value in zip(SLIDERCRANK_CONSTRAINTS, beta, strict=True)
    }
    independent = float(1.0 - np.prod([1.0 - value for value in per_constraint.values()]))

    finite = np.isfinite(beta)
    if not np.any(finite):
        system = 0.0
    else:
        reduced = correlation[np.ix_(finite, finite)]
        ridge = reduced + 1.0e-8 * np.eye(reduced.shape[0])
        try:
            safe = float(
                multivariate_normal(mean=np.zeros(ridge.shape[0]), cov=ridge).cdf(beta[finite])
            )
        except Exception:
            safe = float(np.prod([1.0 - value for value in per_constraint.values()]))
        system = float(min(max(1.0 - safe, 0.0), 1.0))
    return SliderCrankReliability(
        mechanism=mechanism,
        grade=grade,
        value=nominal,
        sigma=sigma,
        correlation=correlation,
        per_constraint=per_constraint,
        system=system,
        independent_bound=independent,
    )


def firing_frequency_sensitivity(
    performance: object,
    vehicle: Vehicle | None = None,
) -> dict[str, float]:
    """Test the comparison's single most load-bearing assumption.

    In this model the EX-link completes all four strokes in one crankshaft
    revolution, because that is what :func:`exlink.cycle.find_phases` demands
    of the piston motion.  A conventional four-stroke needs two.  The EX-link
    therefore accumulates, *per cycle*, half the journal rotation and half the
    piston sliding distance -- and that, not extended expansion, turns out to
    be the larger part of its advantage.

    So the comparison is re-run as if the same linkage drove a conventional
    four-stroke gas exchange: twice the friction per cycle, and half the power
    at the same speed, since an engine that fires half as often makes half the
    power.  Both are re-scored through the vehicle rather than compared as
    efficiencies, because halving the power moves the operating point the
    burn-and-coast strategy can use, and that is where the effect shows up.

    If the advantage survives, it comes from the thermodynamics; if it does
    not, it comes from the firing frequency, and the conclusion has to be
    stated as being about firing frequency.  Measured against
    :func:`optimise_slidercrank`'s baseline it does not survive -- it reverses
    -- which is why that baseline has to be optimised rather than assumed.

    Args:
        performance: An :class:`~exlink.performance.Performance` to re-score.
        vehicle: The car; a default Prototype-class entry if omitted.

    Returns:
        ``brake_efficiency`` and ``km_per_litre`` as modelled, and the same two
        under the doubled-friction assumption.

    Raises:
        ValueError: If the performance carries no friction breakdown.
    """
    car = vehicle if vehicle is not None else Vehicle()
    loss = getattr(performance, "friction", None)
    if loss is None:
        msg = "cannot re-score a design that was never sized"
        raise ValueError(msg)

    quantity = performance.heat_release  # type: ignore[attr-defined]
    speed_rpm = performance.speed_rpm  # type: ignore[attr-defined]
    mass = performance.engine_mass_kg  # type: ignore[attr-defined]

    def score(brake_work: float, revolutions: float) -> tuple[float, float]:
        efficiency = brake_efficiency(brake_work, quantity)
        power = brake_work / 1000.0 * (speed_rpm / 60.0) / revolutions
        outcome = best_strategy(car, mass, power, efficiency)
        return efficiency, outcome.km_per_litre

    as_modelled = score(loss.brake_work, 1.0)
    # Two revolutions per cycle: twice the sliding and twice the rotation for
    # the same indicated work, so twice the friction.
    doubled = score(loss.indicated_work - 2.0 * loss.total_work, 2.0)
    return {
        "brake_efficiency": as_modelled[0],
        "km_per_litre": as_modelled[1],
        "brake_efficiency_four_stroke": doubled[0],
        "km_per_litre_four_stroke": doubled[1],
    }
