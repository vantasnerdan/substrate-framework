# Target reconstruction and constructive nonlinear routes

The user asks for a campaign audit and public discussion restoring Federico's
nonlinear rotational-Euler objective. This is separate from claiming that the
full nonlinear theory has now been derived. Source pin: main b6fc902,
v0.183.0. The exact Vikulin paper Federico intended is still unconfirmed;
the primary candidates below are explicitly provisional identifications.

## Historical target and present meaning

Issue198's original intake asks for Euler coarse-graining with vortex-carried
rotational degrees of freedom. The actual September3 comment
https://github.com/vantasnerdan/substrate-framework/issues/198#issuecomment-5526924549
specifies a linear micropolar target. Later broad completion messages in
issue198 assert the straight-tube formula alpha=L_v*T/6; that route was
subsequently repaired/replaced, and the final registry does not accept it.
Issue200 and C-CST-018 expressly state prepared nonzero-K finite-window
linear response/second variation, not an autonomous finite-amplitude law.
No occurrence of Vikulin was found in the P251 source corpus or its effort
memory. That search result concerns this repository corpus, not every private
conversation or unrecorded author intention.

Federico's later quoted request changes the scientific target: identify and
derive finite-rotation nonlinear dynamics, including the intended stress and
angular-current representation, from the Euler fluid. The existing linear
preparation cannot supply that merely because its microscopic background is
an exact solution of nonlinear Euler.

## Primary-source distinctions

1. Enciso and Peralta-Salas, *Existence of knotted vortex tubes in steady
   Euler flows*, arXiv1210.6271 / Acta Mathematica214 (2015),61-134,
   https://arxiv.org/abs/1210.6271, Theorem1.1: thin prescribed tubes can be
   realized in a decaying Beltrami field, with rich invariant vortex geometry.
   This supplies a stationary geometric backbone; it does not derive a
   rotational constitutive law or a compact-velocity finite-density assembly.
2. Böhmer, Downes and Vassiliev, *Rotational elasticity*, arXiv1008.3833,
   https://arxiv.org/abs/1008.3833, abstract and assumptions1-2: a Cosserat
   special case already permits geometrically nonlinear rotations while
   choosing a physically quadratic constitutive energy. Therefore 'Cosserat
   theory is linear' describes the approximation P251 used, not the entire
   mathematical class. Finite kinematics and nonlinear material law are
   different questions.
3. Erofeev, Pavlov and Vikulin, *Do rotational waves really exist?* (2018),
   https://mpm.spbstu.ru/en/article/2018.59.7/, section5/equation3: Vikulin's
   alternative uses rigid blocks, their intrinsic moments, and the elastic
   host stress. It has symmetric stress and a one-dimensional sine-Gordon
   block-chain reduction. This is a source for a candidate target, not a
   theorem deriving it from incompressible Euler or a complete arbitrary
   three-dimensional finite-rotation continuum.
4. Vikulin et al., *On Wave and Rheidity Properties of the Earth's Crust*,
   Physics of the Solid State58 (2016),561-571,
   https://doi.org/10.1134/S1063783416030306: the primary author's full paper
   is also discoverable via its author-hosted ResearchGate record. It assumes
   an elastic block geomedium, computes host elastic stress and gives a
   sine-Gordon reduction in section5. Its elastic material inputs cannot be
   imported as if Euler had already derived them. Confirmation of Federico's
   intended paper/equations remains useful before a dedicated equivalence proof.

The sources support the distinctions above; their broader geophysical or
fundamental interpretations are not adopted as framework premises.

## Exact reasons the linear result does not select a nonlinear theory

For a relative angle theta, the periodic objective candidate potentials

    V1=k(1-cos theta),
    V2=k(1-cos theta)+beta(1-cos theta)^2

have identical linear stiffness V1''(0)=V2''(0)=k. Their mechanical torques,
defined as -dV/dtheta, differ by
-2beta(1-cos theta)sin theta=-beta theta^3+O(theta^5). These are admissible
finite-angle alternatives depending only on relative orientation; an
ambient-frame-relative SO(3) trace realizes the same example. Thus even a
correct linear optical gap does not select sine-Gordon or its nonlinear
coefficient. `root_nonlinear_separation.py` differentiates these expressions
and verifies the separation, rather than taking an expected force as input.
The independent audit corrected the first draft's potential-gradient/torque
sign convention; the distinct nonlinear laws and common stiffness are unchanged.

Material orientation alone also need not close under Euler. If
C=<r r^T> is the covariance of a transported material blob and
c=u-X_dot, Reynolds transport gives

    C_dot=<c r^T+r c^T>.

Rigid rotation c=Omega r would give C_dot=Omega C-C Omega. The exact steady
incompressible Euler strain u=(a x,-a y,0),p=-a^2(x^2+y^2)/2 instead gives
C(t)=diag(exp(2at),2exp(-2at),3) from C(0)=diag(1,2,3): its eigenvalues
change. The exposing script checks the full Euler residual and unit flow-map
determinant. This is a counterexample to unrestricted rigid-orientation
closure from local Euler equations alone. The affine field is unbounded,
nondecaying and nonperiodic, with infinite whole-space energy; it is not a
counterexample within a compact, periodic or finite-energy subclass, nor to
P251's specially prepared finite-window linear theorem. A nonlinear route needs evolving shape/
inertia variables or a proved slaving estimate, rather than assuming them away.

## Three constructive route families

### A. Exact nonlinear moment and memory equations first

For constant-density Euler on a fixed filtering scale, define a commuting
spatial filter and U=bar u. The exact momentum balance is

    partial_t U+div(U tensor U)+grad bar p=-div tau,
    tau=overline(u tensor u)-U tensor U.

The subfilter tensor tau is symmetric. Finite-cell intrinsic angular moment,
pressure torque, covariance and shape moments supply the rotational balance.
Their evolution is generally a hierarchy. Derive these exact nonlinear
relations with actual moving-boundary/pressure terms before choosing a local
closure. A projection with unresolved variables retained gives exact memory;
a local constitutive model then requires a separately proved memory decay,
scale limit, or invariant/slaved family. Symmetric microscopic momentum flux
and an effective asymmetric stress with an intrinsic-spin/couple split are
compatible when the full angular current is transformed consistently.

First executable result: one exact finite-cell nonlinear angular/covariance
hierarchy, compared against two Euler solutions with the same chosen coarse
state and different unresolved shape/stress. If the coarse derivatives differ,
add the missing state or explicit memory. This identifies the minimum physical
state for the requested theory, without guessing a nonlinear elastic energy.

### B. A controlled finite-rotation reduction on an actual vortex family

Use deformation chi, F=grad chi and a physical orientation R in SO(3), with
body angular velocity Omega=R^T D_t R and curvature Gamma_i=R^T partial_i R.
Their exact kinematic compatibility includes

    partial_i Gamma_j-partial_j Gamma_i+[Gamma_i,Gamma_j]=0.

For delta R=R eta, retain delta Omega=D_t eta+[Omega,eta] and
finite-rotation commutators; in material coordinates the analogous ordinary
time derivative is exact. Compute the full Euler action/kinetic metric on a
specified volume-preserving vortex family, including ambient fluid and
shape/reaction variables. Derive the nonlinear moment balance and coefficients
from that restricted action. Then estimate the residual normal to that family
under actual Euler evolution and prove its control for finite amplitudes on
a declared time/scale regime. Action restriction alone is not dynamical closure.

Geometry choices: repair the actual bordered Hanzawa inverse in0263, or use
a published compact Euler shell whose existence does not depend on that
unproved custom inverse. The second option requires a new response/observable
construction on that very shell; no automatic transfer of the Bessel-core
normalizer is claimed. An EPS Beltrami family is another geometric starting
point with its actual noncompact velocity tails retained.

First executable result: derive cubic and quartic action/response terms on
one justified background and independently test whether their normal residual
is controlled. This decides whether a local finite-rotation limit is plausible
before another large theorem-promotion campaign. A finite-window theorem does
not require global-in-time 3D Euler regularity.

### C. A targeted Euler-to-Vikulin bridge

After confirming the primary target, freeze its variables, finite-angle
potential, inertia, elastic-host assumptions, stress and couple-current
convention, dimension, and approximation scale. Derive the same block/cell
inertia and interaction from the Euler action and literal angular current.
In a one-dimensional invariant axis sector, check whether the finite-angle
restoring law is actually sinusoidal and whether its stiffness/length/time
scales follow from those Euler quantities. A matching linearization is
insufficient, as the exact two-potential example demonstrates.

First executable result: an equation-by-equation dictionary with an Euler-
derived finite-angle interaction or an explicit mismatch. A block-chain
sine-Gordon result is a useful bounded target; a full 3D nonlinear field theory
requires its own SO(3), shape and stress closure beyond that reduction.

## Route assessment

These are plausible construction programs, not a proof that the entire
nonlinear objective is achievable. RouteA gives the most reliable first
nonlinear object; routeB offers the strongest direct continuum result if
geometry and normal-residual control can be constructed; routeC defines the
specific Vikulin comparison and prevents a target substitution. They can
share exact moment/action calculations while keeping their missing licenses
separate. Neither a failed finite-R route nor a linear preparation is a
no-go for the broader objective.
