# Localized filament wave and an exact restoring deficit in its derived model

This candidate changes the geometric setting to a localized excitation of a
vortex background. The exact results below concern binormal curvature flow
and its Hasimoto equation. Their Euler transfer is stated separately; no
filament approximation is presented as an exact finite-core Euler particle.

## 1. Execute the embedded geometric soliton

Use geometric time units in which gamma_t=gamma_s cross gamma_ss and s is
arclength. Choose eta>0, xi real, set y=eta(s-2xi t),
theta=xi s+(eta^2-xi^2)t and A=2eta/(eta^2+xi^2). Define

    gamma=(s-A tanh y, A sech y cos theta, A sech y sin theta). (1)

Direct differentiation gives

    |gamma_s|=1, gamma_t=gamma_s cross gamma_ss,
    curvature=2eta sech y, torsion=xi.

Its Hasimoto field is

    psi=2eta sech y exp(i theta),
    i psi_t+psi_ss+|psi|^2 psi/2=0.                     (2)

This is a finite-amplitude exact geometric solution. Taking xi^2>eta^2 makes
gamma_1,s>=1-2eta^2/(eta^2+xi^2)>0, hence the curve is an embedded graph over
the straight axis. Its transverse displacement and curvature are exponentially
localized. The finite excess of arclength over axial projection is

    Delta L=int(1-gamma_1,s)ds=4eta/(eta^2+xi^2).        (3)

The planar xi=0 member is a poor finite-core candidate: its symmetric loop
self-intersects because eta s=2tanh(eta s) has a positive root as well as zero.
The embedded xi^2>eta^2 subclass avoids that concrete defect without changing
the model or fitting a particle observable.

## 2. Exact positive nonlinear deficit, derived rather than assumed

For a decaying H1 field psi on the line define

    M=int |psi|^2, P=Im int conjugate(psi) psi_s,
    H=int (|psi_s|^2-|psi|^4/4), M>0.

Let q=exp(-i P s/M)psi, F(s)=int_-infinity^s |q|^2, W=F-M/2. Then

    D=H-P^2/M+M^3/192
     =int |q_s+(W/4)q|^2 ds >=0.                      (4)

To prove it, expand the square. The cross term integrates as
int W (|q|^2)'/4=-int |q|^4/4. Also
int W^2 |q|^2 ds=int_0^M(F-M/2)^2 dF=M^3/12.
The Galilean phase removes exactly P^2/M from the gradient energy. Smooth
decaying fields justify the calculation directly; approximation gives the
usual finite-mass H1 version. There is no fitted coefficient or numerical
spectral assertion.

Equality in(4) solves q_s=-(F-M/2)q/4. Its nonzero solutions are exactly

    q=(M/4) sech[M(s-s0)/8] exp(i theta0),              (5)

so(1)--(2) realizes the minimizer at M=8eta and P=8eta xi. The actual values are

    H=8eta xi^2-(8/3)eta^3, D=0.

Under the focusing NLS equation(2), M,P,H and hence D are conserved by direct
integration by parts. This is a nonlinear restoring deficit around a movable
localized wave in this geometric model; it is not an Euler energy bound.

## 3. A quantitative shape observable controlled by that deficit

For a regular strictly positive density, use cumulative mass F as coordinate,
write a(F)=|q(s(F))|^2 and phase vartheta(F). The ground-state density is

    a_*(F)=F(M-F)/4.

Changing variables in(4) gives the exact identity

    D=(1/4)int_0^M |a_F-a_*,F|^2 dF
                         +int_0^M a^2 vartheta_F^2 dF. (6)

Since a and a_* vanish at the mass endpoints, one-dimensional Poincare and
the fundamental theorem of calculus imply

    ||a-a_*||_L2(0,M)^2 <=4M^2 D/pi^2,
    ||a-a_*||_infinity <=2sqrt(M D).                    (7)

Fix 0<delta<1/2 and let a_min=delta(1-delta)M^2/4. If
2sqrt(MD)<a_min/2, then a>=a_min/2 on
[delta M,(1-delta)M]. Align the physical median s(M/2) with that of the
soliton. Since s_F=1/a,

    |s(F)-s_*(F)| <= 4M sqrt(M D)/a_min^2               (8)

on this central mass interval. Thus the conserved deficit controls the
nonlinear central-mass shape and its associated phase-gradient cost. It
does not assert uniform control of arbitrarily remote tails or silently turn
a quantile norm into full Euler orbital stability. Equations(6)--(8) use a
positive-density chart; general nodal perturbations retain(4) but require a
separate compactness/generalized-quantile argument for this observable.

## 4. Actual source transfer to Euler

[Jerrard--Seis1603.00227v1](https://arxiv.org/abs/1603.00227v1), Theorem2,
controls a closed Euler filament relative to a smooth binormal curve if
vorticity concentration persists at every time, with fixed curve length and
regularity bounds. Its error is order sqrt(k_epsilon), with
k_epsilon=4pi/|log(epsilon/L)|. It does not prove that concentration hypothesis.
Our straight-background soliton is nonclosed, so a relative-energy/noncompact
version is additionally needed. The source's curve/tangent norm does not
control the curvature/torsion derivatives entering(4). Those are specific
remaining transfer estimates, not a refutation of the geometric soliton.

[Davila et al2007.00606](https://arxiv.org/abs/2007.00606), Theorems1--2,
provide actual smooth all-time rotating/translating helical Euler filaments
with axial periodicity. They establish a finite-core background supplier,
not a localized envelope or its neighboring-state stability. A nonzero
periodic field repeated along an infinite axis is not a finite-total-energy
isolated object; finite excess energy must be proved relative to a specified
background.

[Cao et al2411.02055](https://arxiv.org/abs/2411.02055), Theorem1.3 and its
doubly connected counterpart, give helical patch waves via bifurcation at
specified angular frequencies, with stated mode/pitch conditions. Their
helical no-swirl condition has velocity orthogonal to the screw direction
while vorticity is parallel to it, so u dot omega=0 pointwise. Geometric
helicity of a screw-shaped axis cannot substitute for physical fluid helicity.
These are exact Euler waves, but neither a localized envelope nor quantum spin.

The classic Hasimoto primary publisher body was inaccessible in this pass.
Equations(1)--(8) were independently differentiated/derived rather than
attributed to unseen text. The modern viscous simulation2602.22439 remains
an inventory lead and has supplied no Euler theorem here.

## 5. Physical energy and continuation

In the filament limit, the Euler kinetic Hamiltonian has leading line-energy
scaling; for a relative straight background the formal leading localized
contribution is proportional to(3) times the finite-core logarithm. The NLS
Hamiltonian in(4) is instead a curvature/torsion functional. Their parameter
dependences differ already on(1). An exact coordinate evolution does not by
itself identify the two actions or make D an Euler conserved quantity.

The geometric soliton, exact nonnegative deficit and quantile-shape control
are established at binormal/NLS scope. The direct declaration that this is
already a robust finite-core Euler object is blocked by an explicit positive
construction: build an actual Euler solution on a justified vortex background
with persistent thin-core localization and sufficient derivative control to
transfer(4)--(8), or construct a localized modulation directly on the exact
helical supplier and derive its actual invariant energy. Electron/neutrino,
quantum action/statistics and current identification remain separate.
