# Appendix A. Theory

The derivations behind `exlink`: the kinematic inversion, the cycle model, the
load analysis, the sizing criteria, and the exact derivatives of all of it.
Symbols are defined where they first appear.

Units throughout are millimetres, radians (degrees only at the API surface),
newtons, $\mathrm{MPa} = \mathrm{N/mm^2}$ and $\mathrm{N\,mm}$.

This appendix carries the derivations only. The argument they support, and every
measured figure, belong to the paper and are not restated here: a second copy
would have to be kept in step with the first, and a derivation that quotes a
stale number is worse than one that quotes none.

| for | see |
|---|---|
| why the conventional objective fails | §2.2 |
| the range objective and the model chain | §3.2, §3.3 |
| the coupling, and MDF against IDF | §4.3 |
| the discrete gear choice | §4.4 |
| tolerance and reliability | §4.5 |
| every measured result | §5 |

## A.1 Parametrisation

The mechanism is a two-shaft linkage. Shaft $R_1$ carries the crank $q_1$ ending
at $Q$; shaft $R_2$ carries $q_2$ ending at $D$. $R_2$ sits at distance $I$ from
$R_1$ along the direction $\theta_r$. A pair of gears of primitive radii

$$r_1 = \tfrac{2}{3} I, \qquad r_2 = \tfrac{1}{3} I, \qquad r_1 / r_2 = 2$$

ties the two shafts, so $R_2$ turns twice for every turn of $R_1$.

**Which shaft is the crankshaft.** The four strokes complete in one turn of
$R_1$, hence in $720^\circ$ of $R_2$ — what a conventional four-stroke
crankshaft does. $R_2$ is therefore the crankshaft, and the shaft power is taken
from; $R_1$ is the half-speed shaft. The analysis is nonetheless parametrised on
$\theta_1$, because one turn of it is exactly one cycle and the closed-form
inversion below is written in it. Everything downstream follows: the code's
`speed_rpm` is $R_1$'s, the engine speeds quoted are $R_2$'s at twice that, and
$M_r$ below is the whole engine torque referred to $\theta_1$, so referring it to
$R_2$ halves it and leaves the power alone.

The swing rod $a$ runs $Q \to A$; the *trigonal link* is the rigid triangle
$A\text{–}D\text{–}E$; the piston rod $e$ runs $E \to P$; the piston crown $H$
sits $p = 16$ mm above $P$, on the cylinder axis $x = x_1$.

**Angle datums.** $\theta_r$ is measured from $+x$, but the two crank angles
$\theta_1$ and $\theta_2$ are measured from $+y$, so that

$$R_1 \to Q = q_1(-\sin\theta_1,\; \cos\theta_1),$$

and $\theta_1 = 0$ puts the crank straight up. That is what the $q_1\sin\theta_1$
and $-q_1\cos\theta_1$ of the loop closure (A.3) are saying, and it is easy to
read past: an arc drawn from $+x$ instead misses the member it is meant to
measure by a right angle.

Describing the triangle by its three sides $b$, $c$, $d$ would force the design
space to respect the triangle inequality, and would leave the sign of $\theta_b$
undetermined. $E$ is therefore placed in the frame carried by $c = AD$:

$$b = \sqrt{x_b^2 + y_b^2}, \qquad
  \theta_b = \operatorname{atan2}(y_b, x_b), \qquad
  d = \sqrt{(x_b - c)^2 + y_b^2}.$$

Now $x_b$ and $y_b$ range freely over $\mathbb{R}$, $\theta_b$ carries its own
sign, and the design space is a plain box. The two routes to $\theta_b$ — this
one and the Carnot expression
$\theta_b = \arccos\bigl((b^2 + c^2 - d^2)/(2bc)\bigr)$ — are checked against
each other in `tests/test_design.py`.

## A.2 Kinematics

Six degrees of freedom $(\theta_1, \theta_2, \theta_T, \theta_a, \theta_e,
\lambda)$ and five constraints, so one input $\theta_1$ fixes everything.

**Gear relation.** External gears turn opposite ways, at the inverse ratio of
their radii:

$$\theta_2 = -2\theta_1 + \theta_f. \tag{A.1}$$

**Loop closures.** Two vector chains close on themselves,

$$R_1 \to Q \to A \to D \to R_2 \to R_1 = 0, \tag{A.2}$$

$$R_1 \to Q \to A \to E \to P \to H \to R_1 = 0.$$

Projecting the first on the axes and isolating the terms in $a$ and $c$,

$$\begin{aligned}
A &= q_1 \sin\theta_1 - q_2 \sin\theta_2 + I\cos\theta_r, \\
B &= -q_1 \cos\theta_1 + q_2 \cos\theta_2 + I\sin\theta_r.
\end{aligned} \tag{A.3}$$

Squaring and adding eliminates $\theta_a$ and $\theta_T$ separately, leaving
their difference $T = \theta_a - \theta_T$:

$$A^2 + B^2 = a^2 + c^2 + 2ac\cos T,$$

$$T = \arccos\!\left(\frac{A^2 + B^2 - a^2 - c^2}{2ac}\right). \tag{A.4}$$

**Compatibility.** That arccosine argument must stay inside $(-1, 1)$. Define

$$\delta_{c1} = \max_{\theta_1}
  \left|\frac{A^2 + B^2 - a^2 - c^2}{2ac}\right|. \tag{A.5}$$

If $\delta_{c1} \ge 1$ for even one crank angle, the four-bar cannot pass that
angle: the shafts rock instead of turning. This is the Grashof condition for the
sub-mechanism $(a, c, q_1, I, q_2)$, written as something an optimizer can read.

Then, with $q = \operatorname{atan2}(a\sin T,\; a\cos T + c)$,

$$\theta_T = \operatorname{atan2}(B, A) - q, \qquad \theta_a = \theta_T + T.
\tag{A.6}$$

**Piston rod.** From the horizontal projection of the second chain,

$$\cos\theta_e = \frac{q_1\sin\theta_1 - a\cos\theta_a
  - b\cos(\theta_b + \theta_T) + x_1}{e}, \tag{A.7}$$

$$\delta_{c2} = \max_{\theta_1} \left|\cos\theta_e\right|. \tag{A.8}$$

**Piston height.** From the vertical projection,

$$\lambda = q_1\cos\theta_1 + a\sin\theta_a + b\sin(\theta_b + \theta_T)
  + e\sin\theta_e + p. \tag{A.9}$$

The constant $p$ makes $\lambda$ the height of the crown $H$ rather than of the
wrist pin $P$; being constant it cancels out of every stroke and volume, so only
the absolute datum depends on it.

**Why the analytic inversion matters.** Newton–Raphson would solve the same
equations. But $\delta_{c1}$ and $\delta_{c2}$ exist only because the inversion
is explicit, and handing those two numbers to the optimizer — instead of letting
the analysis diverge — is what turns a problem full of hard failures into one
with a usable search landscape. It is the single most consequential modelling
decision here.

**Critical configurations.** $|\cos T| = 1$ means $a$ and $c$ are parallel: the
swing rod stops working as a rod and the mechanism gains a degree of freedom.
The constraint is set at $W \le 0.985$, keeping $T$ inside
$[10^\circ, 170^\circ]$. The other pair, $\theta_e \in \{0, \pi\}$, is already
excluded by the $10^\circ$ rod-angle limit.

Verified in `tests/test_kinematics.py`: every link keeps its length to
$10^{-9}$ mm, and $P$ and $H$ stay on $x = x_1$, over the whole revolution.

## A.3 The Atkinson cycle

Over one revolution of $R_1$ — one full cycle — $\lambda(\theta_1)$ must have
four monotone phases: two maxima (top dead centre, reached twice) and two
*different* minima. The deeper minimum ends expansion, the shallower one ends
intake.

| from | to | stroke |
|---|---|---|
| TDC | deep BDC | expansion |
| deep BDC | TDC | exhaust |
| TDC | shallow BDC | intake |
| shallow BDC | TDC | compression |

$$\mathrm{STE} = \lambda_{\mathrm{TDC}} - \lambda_{\mathrm{deep}}, \qquad
  \mathrm{STC} = \lambda_{\mathrm{TDC}} - \lambda_{\mathrm{shallow}},$$

$$g = \left|\lambda_{\mathrm{TDC},1} - \lambda_{\mathrm{TDC},2}\right|,$$

with $\lambda_{\mathrm{TDC}} = \max(\lambda_{\mathrm{TDC},1},
\lambda_{\mathrm{TDC},2})$ — the higher of the two, since that is the one which
sets the clearance volume. That maximum is not differentiable where the two
coincide, which §5.3 shows is not a technicality.

`find_phases` counts sign changes of $\mathrm{d}\lambda/\mathrm{d}\theta_1$ and
refines each extremum with a parabola through its three neighbouring samples —
necessary because $g$ is constrained at $0.01$ mm on a $0.5^\circ$ grid.

Volume and pressure follow, with $A_p = \pi\Phi^2/4$ the bore area:

$$V = V_0 + A_p(\lambda_{\mathrm{TDC}} - \lambda), \qquad
  V_1 = V_0 + A_p\,\mathrm{STC}, \qquad \varepsilon = V_1/V_0,$$

$$P = \begin{cases}
P_0 & \text{intake, exhaust},\\
P_0 (V_1/V)^\gamma & \text{compression},\\
P_3 (V_0/V)^\gamma & \text{expansion},
\end{cases}
\qquad P_3 = k P_2, \quad P_2 = P_0\varepsilon^\gamma,$$

with combustion instantaneous at top dead centre and blow-down instantaneous at
the deep bottom dead centre. The gas force on the crown uses the *gauge*
pressure, $P_{\mathrm{gas}} = (P - P_0)A_p$, so intake and exhaust are unloaded.

Note the sense of the compression exponent: $P = P_0(V_1/V)^\gamma$, so pressure
rises as the charge is compressed, and at $V = V_0$ it gives
$P_2 = P_0\varepsilon^\gamma$ as it must.

Given $\Phi = 32$ mm, $V_0 = 3$ cc and $\varepsilon = 16$, the compression stroke
is pinned at $\mathrm{STC} = 15 V_0 / A_p \approx 55.95$ mm against
$\mathrm{STE} = 74$ mm — the asymmetry the linkage exists to produce.

Designs failing the phase test are penalised rather than rejected, with
$\eta = 0$ and $H = B = 1000$: a piston that goes up and down once per
revolution is a plain Otto engine. Both failure modes are exercised in
`tests/test_model.py`.

## A.4 Quasi-static loads

Inertia is neglected here — this is a first sizing iteration, and the masses are
not known until the parts have a shape.

**Piston.** With $P$ the gas force and $\theta_e$ the rod angle,

$$C = \frac{P}{\sin\theta_e} \quad\text{(rod load)}, \qquad
  D = P\cot\theta_e \quad\text{(side load, reacted by the liner)}.$$

**Trigonal link.** Force balance at $A$, $D$ and $E$ plus the moment about $D$.
Writing $\overrightarrow{DE} = b\,u(\theta_b + \theta_T) - c\,u(\theta_T)$ and
$u_e = (\cos\theta_e, \sin\theta_e)$,

$$\overrightarrow{DA} \times F_A + \overrightarrow{DE} \times F_E = 0
  \;\Longrightarrow\;
  -cA\sin(\theta_a - \theta_T) - C(\overrightarrow{DE} \times u_e)_z = 0,$$

$$A = \frac{-C\,(\overrightarrow{DE} \times u_e)_z}{c\sin(\theta_a - \theta_T)},
\tag{A.10}$$

$$Q_x = C\cos\theta_e - A\cos\theta_a, \qquad
  Q_y = C\sin\theta_e - A\sin\theta_a.$$

The leading minus sign in (A.10) matters. Drop it — an easy slip when inverting
the moment equation — and the computed torque disagrees with the principle of
virtual work by a factor of about $-4$, with the efficiency coming out negative.
With the sign above, agreement is exact; §A.6 is what catches it.

Note $\sin(\theta_a - \theta_T) = \sin T$ in the denominator: the swing-rod load
diverges at exactly the critical configurations (A.5) excludes — a second,
independent reason to enforce it.

**Shafts.** With $\alpha = 20^\circ$ the standard involute gear pressure angle,

$$T_{\mathrm{gear}} =
  \frac{-q_2 (Q_y \sin\theta_2 + Q_x \cos\theta_2)}{r_2\cos\alpha},$$

$$M_r = q_1 A\cos(\theta_a - \theta_1) + r_1 T_{\mathrm{gear}}\cos\alpha.
\tag{A.11}$$

The two gear torques are $r_1 T\cos\alpha$ and $r_2 T\cos\alpha$ with the *same*
sign, so with $\omega_2 = -2\omega_1$ and $r_1 = 2r_2$ the pair transmits no net
power — checked explicitly in `tests/test_loads.py`, since a gear pair that
generated power would silently inflate the efficiency.

The bearing reactions carry the tooth-force direction in full,

$$R_{1x} = A\cos\theta_a + T\sin(\theta_r + \alpha), \qquad
  R_{1y} = A\sin\theta_a - T\cos(\theta_r + \alpha),$$

so the line of action rotates with the shaft axis $\theta_r$, not with $\alpha$
alone.

## A.5 The conventional objectives and constraints

**Efficiency.** The average mechanical efficiency is

$$\eta = \frac{\oint M_r \,\mathrm{d}\theta_1}
  {2(\mathrm{STE} + \mathrm{STC})\langle P\rangle}
  = \frac{\langle M_r\rangle}{\langle P\rangle}
    \cdot \frac{\pi}{\mathrm{STE} + \mathrm{STC}}, \tag{A.12}$$

a ratio of two works — the torque's on $R_1$ over the gas force's on the piston.
It measures the linkage's aptitude for turning force into torque, and it grows
without bound as the mechanism grows. That unboundedness is exactly why $H$ and
$B$ must enter the problem.

A companion measure $\varphi = \langle M_r\rangle / \langle p\rangle$, the
torque per unit pressure, has the dimension of a volume and is proportional to
$\eta$ at fixed stroke, so it adds no independent information.

**Envelope.** $H$ along the stroke and $B$ across it: the bounding box of every
body over every configuration — joints, both gear primitives, and the piston
over its full travel.

**Clearance.** The trigonal link must stay $10$ mm clear of the cylinder.
`cylinder_clearance` models the liner as the half-strip
$x \in [x_1 \pm \Phi/2]$, $y \ge y_{\mathrm{bottom}}$, and minimises the distance
from the three triangle edges over the revolution. It is monotone in the right
direction and vanishes on contact, which is what the constraint needs; its
numerical value is not a detailed CAD clearance.

**The conventional formulation.** Collecting these,

$$\begin{aligned}
\min_{X} \quad & f(X) = (-\eta,\; H,\; B)^{\mathsf{T}} \\
\text{s.t.}\quad & c(X) = (\mathrm{mra} - 10,\; W - 0.985,\; g - 0.01,\;
  10 - d,\; \gamma - 0.02)^{\mathsf{T}} \le 0, \\
& c_{\mathrm{eq}}(X) = (\mathrm{STE} - 74,\; \varepsilon - 16)^{\mathsf{T}} = 0,\\
& l_b \le X \le u_b,
\end{aligned} \tag{A.13}$$

with $W = \max(\delta_{c1}, \delta_{c2})$ and $\gamma = \max(D)/\max(P)$. §3
replaces it.

## A.6 Verification

The force chain is not merely transcribed; it is pinned by an independent
identity. In a massless, frictionless, quasi-static mechanism the instantaneous
power in equals the power out:

$$M_r(\theta_1)\,\omega_1 = P(\theta_1)\,v_{\mathrm{piston}}
  \quad\Longrightarrow\quad
  M_r = -P\,\frac{\mathrm{d}\lambda}{\mathrm{d}\theta_1}. \tag{A.14}$$

Every step of the chain — piston, trigonal link, both shafts, the gear pair —
must conspire to satisfy this at *every* crank angle, and `exlink` reproduces it
to machine precision (`tests/test_loads.py`). That is what exposed the sign slip
in §A.4.

A second, independent route: the mean torque must equal the indicated $p$–$V$
loop area divided by $2\pi$. It is also checked, and it never touches the force
chain. Together with the rigid-link check of §A.2 these leave very little room
for the model to be wrong in a way the tests would not see.

## A.7 Sizing the parts, and the coupling that creates

A quasi-static study cannot size the parts, and the obstruction is circular: the
inertia loads need the part masses, the masses need the cross-sections, and the
cross-sections need the loads. `exlink.dynamics`, `exlink.sizing` and
`exlink.coupled` close that loop:

$$\text{diameters} \;\longrightarrow\; \text{member masses}
  \;\longrightarrow\; \text{inertia forces}
  \;\longrightarrow\; \text{internal loads}
  \;\longrightarrow\; \text{diameters}.$$

Neither half can go first. That is a genuine multidisciplinary coupling, and it
has to be solved rather than sequenced.

**Why sequential elimination cannot be extended.** Without inertia every rod is a
two-force member: the forces at its two ends are equal, opposite and collinear
with the rod, which is exactly what lets the loads be eliminated one body at a
time from piston to shafts. Add mass and that collapses — a rod with a
distributed d'Alembert load has end forces that are neither collinear nor equal,
so no body can be solved before its neighbours.

Counting with Grübler over 7 links (6 moving plus ground), 8 lower pairs (7
revolutes and the piston's guide) and 1 higher pair (the gear mesh),

$$M = 3(7 - 1) - 2(8) - 1 = 1,$$

so the load problem is statically determinate: 18 unknowns against 6 bodies of 3
equilibrium equations each. `exlink.dynamics.solve` assembles and solves that
$18\times18$ system at every crank angle. Its determinant is the same quantity
(A.5) protects: at a critical configuration the mechanism gains a degree of
freedom, the matrix goes singular and the internal forces diverge. The condition
number is reported so that the connection is visible rather than implied.

**Accelerations.** Every history is smooth and periodic on a uniform grid, so the
Fourier derivative is exact where a finite difference is only $O(h^2)$. That
matters because accelerations are *second* derivatives: on a $0.5^\circ$ grid a
finite-difference second derivative loses about six digits, enough to pollute the
inertia forces and, through them, the sizing loop. Angles that accumulate whole
turns, $\theta_2 = -2\theta_1 + \theta_f$, are split into a ramp plus a periodic
part first (`exlink.derivatives`).

Constant shaft speed is assumed throughout, which simplifies the bookkeeping in
two ways worth stating: both shafts turn at a constant rate, so neither has
angular acceleration nor contributes an inertia couple; and shafts and gears are
concentric with their own axes, so their centres of mass do not move and they
exert no inertia force at all. Only the offset crank throws do, which is why the
modelled bodies are the crank arms rather than whole shaft assemblies.

## A.8 Failure modes and cross-sections

Each link is a round bar whose outer diameter is solved for — the smallest
section satisfying all three modes over the whole revolution:

| mode | criterion |
|---|---|
| static | peak fibre stress against $S_y/n_y$; the state is uniaxial, so von Mises reduces to $\lvert\sigma\rvert$ |
| fatigue | Goodman, $\sigma_a/S_e + \sigma_m/S_u \le 1/n_f$, with a Marin-corrected $S_e = k_a k_b k_c k_d k_e \cdot 0.5 S_u$ |
| buckling | Euler on the peak compressive load, with $K = 1$ for links and $2$ for a crank throw |

Fatigue is evaluated per extreme fibre, so that $\sigma_a$ and $\sigma_m$ are
taken at a fixed material point rather than at whichever fibre happens to be
worst at each instant. A compressive mean stress is not credited as beneficial —
the Goodman line is truncated at $\sigma_m = 0$ — which matters here because the
connecting links swing between tension and compression every revolution. The
required diameter comes from bisection rather than a closed form, because the
fatigue size factor $k_b$ itself depends on the diameter being solved for.

**Solid or hollow.** The study is reported with solid bars, but the members that
would really be made from tube — the rods and the trigonal link, not the crank
throws — may be bored to a ratio $k = d_i/d_o$. Writing $R^2 = r_o^2 + r_i^2$,
that changes four constants and nothing else:

$$A \propto 1 - k^2, \qquad Z \propto 1 - k^4, \qquad I \propto 1 - k^4, \qquad
  \frac{m(3R^2 + L^2)}{12} \text{ with } R^2 \propto 1 + k^2.$$

Every quantity proportional to a power of $d$ keeps that power, so the analytic
derivatives of §A.10 are unchanged in form — the shape factor cancels out of
$\mathrm{d}(1/A)/\mathrm{d}d = -2/(Ad)$ and
$\mathrm{d}(1/Z)/\mathrm{d}d = -3/(Zd)$. A wall floor $2t/(1-k)$ bounds the outer
diameter below, which is what stops the model discovering free mass at small
sections.

**Internal loads.** For a member spanning two joints, the force and moment at a
section a fraction $s$ along it follow from the free body of $[0, s]$. Because a
rigid body's acceleration varies *linearly* along any straight line through it,
that load is linear in $s$ and both integrals close in form:

$$F(s) = m\left[a_1 s + (a_2 - a_1)\frac{s^2}{2}\right] - F_1,$$

$$M(s) = -m\left[(\Delta r \times a_1)_z \frac{s^2}{2}
  + (\Delta r \times \Delta a)_z \frac{s^3}{6}\right]
  + s\,(\Delta r \times F_1)_z.$$

**One idealisation.** The trigonal link is a single rigid part, so as a frame it
is three times statically indeterminate. It is treated as a pin-jointed triangle:
each side takes an axial force from joint equilibrium, which *is* determinate,
and bends only under its own distributed inertia as a simply supported beam. That
captures the axial load exactly and under-estimates bending at the corners.

Shafts are not sized. Their bearing span is not a design variable, so any section
would be arbitrary, and being concentric with their axes they do not feed the
inertia loop. The bearing reactions are computed and constrained.

## A.9 Convergence of the sizing loop

A scaling argument settles it. A bending-critical member needs
$d \sim F^{1/3}$, so its mass goes as $m \sim d^2 \sim F^{2/3}$, and the inertia
force it then creates is $F \sim ma$. Composing,

$$m \sim (Ca)\,m^{2/3},$$

so the loop gain is sub-linear, a fixed point exists at

$$m = (Ca)^3, \tag{A.15}$$

and plain Gauss–Seidel reaches it — slowly, and with a mass growing as the *cube*
of the acceleration level, hence as the **sixth power of engine speed**. That is
why `MDAGaussSeidel` is the right MDA here (no coupled Jacobians are needed, and
the bisection has no useful derivative to give a Newton method anyway), and why
engine speed, which does not appear in the quasi-static problem at all, becomes
one of the strongest drivers once inertia is present.

Two consequences of constant speed are worth separating, both checked in
`tests/test_dynamics.py`. The **mean torque does not move**: at constant speed the
mechanism returns to its starting state every revolution, so its kinetic energy
is unchanged and the inertia forces do no net work. Efficiency, being a
mean-torque quantity, is therefore untouched by speed. Only the peaks change —
and the peaks are what size the parts. With the gas load switched off, every
reaction scales as exactly $\Omega^2$.

## A.10 Derivatives

The feasible set is a sliver: at the reference design the two equalities leave a
band $0.1$ mm wide on $\mathrm{STE}$ and $0.1$ wide on $\varepsilon$ inside an
eleven-dimensional box with sides of tens of millimetres, while $W$ and $\gamma$
sit within $0.4\,\%$ and $7\,\%$ of their bounds. A derivative-free method cannot
work in that. The whole analysis chain is closed form, so its derivatives are
too.

**Forward-mode chaining.** Each intermediate carries $\mathrm{d}/\mathrm{d}X$ as
an $(n_\theta, 11)$ array, differentiated straight through the closed forms of
§A.2. One pass produces all eleven components.

**The envelope theorem.** Most metrics are extrema over the crank angle. For a
maximum attained at $\theta^\star$,

$$\frac{\mathrm{d}}{\mathrm{d}X}\max_\theta f(X, \theta)
  = \left.\frac{\partial f}{\partial X}\right|_{\theta^\star}, \tag{A.16}$$

because the term through the moving maximiser carries
$\partial f/\partial\theta = 0$, so no derivative of the *location* is needed.

That is not a convenience — it removes a failure mode. Because these are maxima,
the sample attaining them **switches** as the design moves, and a difference
quotient taken across the switch is wrong: on $\gamma$ at a $10^{-4}$ mm step, by
$25\,\%$. Both `tests/test_jacobian.py` and `tests/test_dynamics_jacobian.py`
pin that behaviour deliberately.

The extremum does not sit on a grid point, so the partial derivative is
interpolated to the parabolically refined location with a three-point Lagrange
quadratic, matching the order of the refinement that located it; linear
interpolation leaves $\mathrm{STE}$ a factor of fifty worse.

The same switching, at the level of *which of two top dead centres is higher*
rather than which grid point attains a maximum, is what §5.3 is about. It is the
one place where the envelope theorem does not rescue the derivative, because
there the two branches do not share a value of $\partial f/\partial\theta$.

**Through the coupling.** `exlink.dynamics_jacobian` differentiates the MDA
itself, so GEMSEO assembles the coupled derivative from local Jacobians instead
of differencing the fixed point:

| piece | route |
|---|---|
| accelerations | $D^2$ is linear, so $\mathrm{d}a/\mathrm{d}p = \Omega^2 D^2(\mathrm{d}r/\mathrm{d}p)$ |
| mass properties | direct, with $m = \rho(\pi d^2/4)(1 - k^2)L$ and mass-weighted centres |
| $18\times18$ solve | $\mathrm{d}x/\mathrm{d}p = A^{-1}(\mathrm{d}b/\mathrm{d}p - (\mathrm{d}A/\mathrm{d}p)x)$, reusing the factorisation |
| internal loads | closed form; the trigonal truss via $\mathrm{d}z = (M^{\mathsf{T}}M)^{-1}M^{\mathsf{T}}(\mathrm{d}r - \mathrm{d}M\,z)$ |
| sizing bisection | implicit function theorem on $U(d, N, M) = 1$ |

Verified against converged central differences: mass properties and
accelerations to round-off, the equilibrium solve to $5\times10^{-7}$, the
internal loads to $4\times10^{-7}$, the sizing by directional derivatives, and
GEMSEO's own `check_jacobian` on both disciplines.

Four quantities are left differenced deliberately: $\eta$, $H$, $B$ and the
clearance. All are smooth, none is tight, and $\eta$ would additionally need
$\mathrm{d}\theta_{\mathrm{TDC}}/\mathrm{d}X$, because the combustion pressure
jump puts moving-boundary terms in its integral.
