# 0038 — direct stationary Beltrami tube action

Parent: P251 / issue #198, continuation contract 0035. Original smooth
stationary incompressible-Euler to Cosserat objective is unchanged.
Owner: Codex `/root/construction_review`, now implementer of this new attempt;
this agent's earlier review does not independently review attempt 0038.
Write surface: this directory only. Base campaign release v0.171.0,
checkpoint 3626fbf; accepted-main context v0.174.0.

## Positive child deliverable and frozen candidates

Construct a directly derived physical angle sector and its translation
coupling for a genuinely stationary smooth Euler tube. Candidate A uses an
explicit axisymmetric cylindrical Beltrami tube with axial flow and exact
isovortical/Jacobi variations. Candidate B uses a nonaxisymmetric Beltrami
polarization so the angle is a physical field coordinate rather than rotation
of an axisymmetric gauge. EPS sources 1210.6271 and 1505.01605 fix the actual
stationary-tube applicability boundary. A plane wave is not tube evidence.

Selection is structural: exact Euler, physical angle observability, admissible
volume-preserving variations, full kinetic/pressure action, coupling to material
translation, and smooth tube realization. No empirical comparator or target
coefficient is used. The admitted variables and symmetries determine the angle
action; a rate norm cannot be substituted for a static angle potential.

Requires: exact stationary field, domain and boundary data, divergence-free
material displacement, and an explicit physical angle map. Pass licenses:
only the exact action and coupling obtained from that same construction.
The result advances N2/N3/N4/N7 only at those explicit licenses. Unclosed
ensemble and homogenization constructions remain parent obligations.

## Analytic specification

Use constant density rho, curl u0=lambda u0, div u0=0, and
p0=pstar-rho |u0|²/2. Distinguish a material displacement xi from Eulerian
velocity delta u=(partial_t+u0·grad)xi-(xi·grad)u0. The stationary Euler
Jacobi action and relabeling quotient are imported from the parent 0037
construction only when their exact statement is available; direct Euler
linearization provides an independent route. Differentiate full fields and
constraints before restricting a proposed invariant sector. If axisymmetry
makes rigid spin unobservable, execute nonaxisymmetric polarization as the
materially different candidate in the same attempt.

All current calculations are exact algebra/calculus. There is no frozen
production numerical remainder, mesh, eigenvalue decision, or numerical
stability claim. A symbolic oracle derives field residuals and exhibits a
wrong-sign or omitted-pressure mutation if a positive action is constructed.
Full small-ratio prescriptions bind before any later spectral numerics.

Status: active. Source inspection and analytic construction follow this
child contract. Route outcomes and first-run evidence are appended below.

## Exact source applicability

Theorem 1.1 of Enciso–Peralta-Salas, [1210.6271](https://arxiv.org/pdf/1210.6271),
constructs smooth stationary Beltrami fields with thin closed vortex tubes of
prescribed knot/link type, up to small diffeomorphism. The construction gives
1/r decay and infinite total energy on R³. Its local tube boundary problem
fixes the harmonic part; the KAM step preserves invariant vortex surfaces.
These are geometrical field-flow statements, rather than a perturbation-action
or Euler linear-stability theorem. Sections 2, 6–8 and Theorem 1.1 were used
to identify these inputs.

Theorems 1.1 and 2.1 of Enciso–Peralta-Salas–Torres de Lizaur,
[1505.01605](https://arxiv.org/pdf/1505.01605), give compact-manifold Beltrami
tubes and a local high-eigenvalue approximation theorem. The torus knot
construction assumes the prescribed set lies in a contractible region.
Our explicit periodic tube below is a direct Euler construction with a
noncontractible axial core; its proof does not invoke that existence theorem.
Local approximation by the theorem does not identify the surrounding action
or convert this periodic example into the prescribed contractible EPS tubes.

## A: smooth axisymmetric tube, full pressure cancellation

The Lundquist field is an exact analytic stationary Euler tube:

    u0 = V(r) e_theta + W(r) e_z,
    V=A J1(lambda r), W=A J0(lambda r),
    V'=lambda W−V/r, W'=−lambda V,
    curl u0=lambda u0, p0=pstar−rho(V²+W²)/2.

The axis is regular, V=O(r), and every cylinder is an invariant vortex
surface. A periodic axial coordinate closes these surfaces into tube tori.
Radial pressure satisfies p0'=rho V²/r. For the divergence-free material
angle displacement xi=r q(z,t) e_theta, the exact Jacobi terms are

    Dt xi = −V q e_r + r(q_t+W q_z) e_theta,
    xi.Hess(p0).xi = rho V² q²,
    L2 density = rho r²(q_t+W q_z)²/2.

Thus the apparent rotation potential from |Dt xi|² cancels exactly. For
uniform static q, the Eulerian velocity perturbation is zero: this is a
relabeling direction. The radial Jacobi residual before the pressure response
is −2 rho V(q_t+W q_z); variation in z can require additional pressure and
shape responses. Such responses can be named in a constrained ensemble.
They are not the reason this particular angle-stiffness derivation fails:
the projected angle coefficient itself is exactly zero.

Route verdict: refuted for the candidate assertion that free axisymmetric
director rotation supplies positive static locking. Mechanism: physical
relabeling and exact pressure cancellation. This is not a verdict on N3.

## B: a physical nonaxisymmetric Beltrami tube and its translation action

Execute the different candidate in which the core cross-section has an
observable orientation. Let physical coordinates be (ell x,ell y,ell z),
each 2 pi ell-periodic, with a,b>0 having velocity units:

    psi=a cos x+b cos y,
    u0=(−b sin y, a sin x, psi),
    curl u0=u0/ell,
    p0=pstar−rho ab cos x cos y.

Direct differentiation proves stationary Euler and u0.grad psi=0. At
(x,y)=(0,0), Hess psi=diag(−a,−b)/ell². Therefore nearby level curves are
closed ovals and their product with the axial circle gives actual smooth
invariant tube tori. This uses a vortex surface, not a plane-wave regression.
When a differs from b, the cross-section is noncircular.

The globally smooth periodic generator

    chi=ell(−sin y,sin x,0), xi=q(t) chi

is divergence-free and is a rigid rotation to first order at the core. Its
static Eulerian perturbation is

    delta u/q=(b−a)(sin x cos y, −sin y cos x, sin x sin y).

It is physical for a != b and becomes gauge for a=b. Average the COMPLETE
Jacobi action over the periodic cell before assigning coefficients. It gives

    M_q/volume=rho ell²,
    K_q/volume=−rho(a−b)²/4,
    L_q/volume=rho ell² q_dot²/2 + rho(a−b)² q²/8.

The single-coordinate gyroscopic integral is zero by periodic advection.
This is an exact conditional action on the stated constrained displacement
ensemble. No all-wave-number equality to microscopic Euler is required to
state it. Its static angle coefficient has the wrong sign for the desired
positive spring, even at that restricted boundary. The result does not claim
Euler instability: pressure, other shape coordinates and the relabeling
quotient affect full dynamics, and the present scalar ensemble is not a
stability analysis.

Uniform material translation is retained in this same action by
xi=X(t)+q(t)chi. Every translation has mass rho times volume. Uniform
translation/angle mass, gyroscopic, and spring cross integrals vanish; the
translation spring tensor also vanishes. This is the required k=0 symmetry
limit, not an inserted independent mass. The exact Euler family is
u_X(t,x)=u0(x−X(t))+X_dot with X_ddot=0 and p_X=p0(x−X(t)). Nonuniform
macro displacement and a microscopic cage remain new conditional ensemble
constructions, rather than consequences of uniform translation alone.

Route verdict: refuted for this explicit generator's positive-angle-spring
candidate, by its exact negative static coefficient. Positive constructions
established along the route are the actual smooth nonaxisymmetric tube,
physical angle map, and its same-action translational mass.

## C: physical axial polarization, exact smooth Euler continuation

Change the physical degree of freedom instead of rechecking the director.
For any smooth axisymmetric stationary V(r),W(r), including the Beltrami
background above, the field

    u=V(r)e_theta+[W(r)+F(r,theta−Omega(r)t)]e_z,
    Omega(r)=V(r)/r

is an exact time-dependent Euler solution with the same radial pressure.
The z-independent axial component is transported by the transverse flow and
does not feed its momentum equation. Choose
F=A(r) cos(m(theta−Omega t)), with A(r)=r^m f(r²) for axis regularity,
or smooth annular support. This is a physical vorticity polarization, since
its transverse gradient changes curl u; it is not an unobservable director.

Nevertheless, its phase is advected, not statically locked. For m=1 the
angle-integrated kinetic-energy increment per r dr dz is

    integral rho[(W+A cos(theta−phi))²−W²]/2 dtheta
      =rho pi A²/2,

which is independent of phi and nonzero for A!=0. Thus physical polarization
and nonzero energy alone still do not derive alpha. The exact phase advection
is positive construction evidence for the original mode question, while the
candidate phase-spring route is refuted by its exact phase independence.

## Evidence, scope and next construction

`verify_stationary_tube.py` derives the fields and full action terms and
rejects both omitted and doubled pressure-Hessian coefficients. Its first
execution used Python 3 in the existing repository venv with PYTHONPATH=src;
process status zero, 30 checks passed. Output is preserved in stdout.txt.
There is no numerical spectrum, approximation floor or stability claim.
These route verdicts have `EXACT_CONDITIONAL_ACTION_AND_EULER_FIELD` scope.
No accepted claim changes or parent completion follow.

The continuation earned by the failures is an interaction-bearing constrained
ensemble: move the noncircular core relative to a separately represented
smooth cage, retain both in the microscopic pressure/kinetic action, and
eliminate its conjugate shape coordinate or retain the resulting coupled
matrix. A stiffness cannot be supplied by the same free-rotation or passive
phase construction again. The smooth six-core interaction worker owns that
parallel construction in 0036. For a direct EPS route the corresponding
object is the harmonic-part-fixed multiple-tube boundary construction with
the core and exterior/cage variations included in the same action. Its local
pressure/coupling coefficients, rather than KAM trajectory stability, are
the quantities that can establish the next positive angle sector.

Child slice status: the three stated candidate computations are complete;
the desired positive locking action remains open on the parent campaign.

Validation note: Ruff is provided by `/home/dan/anaconda3/bin/ruff`, not
the repository venv. Its first actual lint pass requested replacing two
assigned lambdas with equivalent local functions. This syntax-only repair
changes no expression or scientific predicate; the recorded scientific run
remains the receipt. The corrected lint pass and `git diff --check` pass.
