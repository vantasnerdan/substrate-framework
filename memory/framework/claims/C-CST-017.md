---
description: Accepted framework claim C-CST-017
author: framework-registry
created: '2026-09-05T20:07:10+00:00'
updated: '2026-09-05T20:07:10+00:00'
tags:
- substrate-framework
- accepted-claim
- C-CST-017
category: claims
confidence: established
status: active
---
# C-CST-017

## Statement
The accepted statement is reproduced exactly from the claim registry.

On the fixed C-CST-015/016 smooth stationary periodic Euler field
u=(cos Y+A cos Z,A sin Z,-sin Y), A=1/100, one positive whole-field
O(3)/time-reversal law and one common finite smooth linear Euler/Lin
preparation join actual hybrid displacement U, covariance angle Phi,
inherited second-variation action and full physical momentum/angular
current in a prepared incompressible isotropic micropolar continuum
through o(|K|^2), with any fixed finite time-derivative inventory on
each fixed compact time window.

The construction uses the actual acoustic displacement/velocity inputs
of C-CST-015/016, both quadratures of their wrapped-streamline passive
Euler sector to repair the complete point-to-hybrid acceleration, the
same-field two-positive-tag optical angle/full-spin controls, and the
complete mixed phase/Jacobi-energy normalizer. The whole initial field
is summed before averaging. All quadratic cross forms and full finite
source norms are retained. Signed controller amplitudes are physical
initial data; every probability and material mass fraction stays positive.

With C=i[K cross], transverse acoustic amplitude A0 and unrestricted
optical amplitude B, the actual observations satisfy
(U,Phi)=T(A0,B)+e,
T=[[I,-j C/(2rho)],[C/2,I]],
A0_tt+a|K|^2 A0=r_A,
B_tt+[nu^2 I+cT|K|^2 P_T+cL|K|^2 P_L]B=r_B.
Here rho is the ambient density, j>0 is the measured whole-law spin
normalization, a>0 is the actual acoustic coefficient and nu>0 is the
selected actual optical clock. The spatial coefficients cT,cL are
explicit prepared output targets with cT>j nu^2/(4rho), cL>0.
The errors e,e_t,e_tt,r_A,r_B are o(|K|^2) in the declared operator
norms along the common ordered preparation/long-wave sequence.

Set mu=rho a, alpha=j nu^2/4,
gammaT=j(cT-alpha/rho), gammaL=j cL and M=diag(rho I_T,j I).
If K2 is the canonical micropolar stiffness, direct substitution gives
M Y_tt+K2 Y=M T r+(K2 T-M T D)z+M e_tt+K2 e.
The exact algebraic mismatch is cubic in K. The acoustic and optical
zero-wave-number state Wronskians equal one, and the five-position
observation determinant is (1+j|K|^2/(4rho))^2. Thus the bounded state
inverse makes the coupled equation a statement for every retained
initial amplitude, not a check along one trajectory. Pulling back BOTH
inherited phase/kinetic and energy forms through the same physical T
gives the same M and K2 through second spatial order.

Literal total momentum retains collapsed material tags and continuous
ambient fluid. Its full stress and couple flux include pressure,
moving-boundary reaction and convective transport. For the actual
time-dependent current coefficient q, use
Q_ij=q epsilon_ijk U_t,k, S_int=S_full-div Q,
N_int=N_full-Q_t, J_int=J_full, F_int=F_full.
Then S_int=j Phi_t+o(|K|^2), while momentum is unchanged. Initial
integrated angular charge and the full q_t memory remain independent
data; a primitive of spin is not silently identified with a material
displacement dipole. The corresponding boundary action keeps its
acceleration and gradient-velocity momentum cancellation.

For canonical F_can=partial W/partial grad U and
N_can=partial W/partial grad Phi,
P_T div(F_int-F_can)=o(|K|^2),
div(N_int-N_can)-ax(F_int-F_can)=o(|K|^2).
These are the actual periodic bulk constitutive virtual-work class,
including the local angular torque. Longitudinal force is the pressure
multiplier; cut-domain surface terms remain explicit. This is not
pointwise equality of arbitrary free-surface tractions.

The same fixed cell also has an analytic finite invariant tube with
nondegenerate elliptic core and a Diophantine boundary with nonzero
flux twist. In coordinates psi=cos b+Omega^2 cos a, Omega=1/10,
the central flux derivative is
r_J(0)=(-1+6Omega^2-Omega^4)/(8(1+Omega^2)^3)
      =-235025/2060602.
The EPS normal torsion on a chosen nearby boundary is 4pi^2 r_J(J_b).
Local analytic divergence-free persistence uses the positive section
flux and its Moser identification. Its reference line density is
1/(2pi)^2, its fixed-radius tube volume fraction is positive, and a
per-axial-length tagged coefficient j0 gives whole-law density
mean_tag_fraction*j0/[3(2pi)^2]. This is the same periodic background;
the material/action normalization remains the separately computed one.

The tube core is noncontractible in T^3 and lifts to an unbounded line,
not a compact Euclidean ring. This theorem does not close that remaining
same-field geometry/density obligation. It also does not assert an
unrestricted Euler invariant manifold, acoustic-time uniformity,
nonlinear finite-amplitude stability, a universal unprepared modulus,
or a single fixed preparation smooth across all switching scales.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: C-CST-009, C-CST-015, C-CST-016. Assumptions: Actual smooth constant-density Euler/Lin on the fixed periodic cell; linear-response and second-variation scope with common laboratory amplitudes and whole-state rotations, inversion and time reversal., Actual positive material tags, regular streamline-band moments, nonzero observation gains and measured optical/spin normalization from the stated constructions. Arbitrary API matrices or coefficient values alone do not construct these suppliers., Full finite controller norms, material shape and second-momentum moments, two time derivatives of observation errors, complete phase/energy cross forms, pressure and current-memory costs enter the common ordered diagonal. The preparation may depend on the requested accuracy and fixed time window., Standard smooth finite-time Euler/Lin Sobolev estimates, exact coarea and polynomial approximation, full Bloch Leray differentiation, Haar averaging and the explicitly cited local EPS/Moser persistence theorem are declared mathematical imports., Exact analytic construction and source residual calculus are the oracle. Symbolic tests verify the finite identities and expose wrong coupling, phase, moment and torque conventions; they do not independently prove the analytic PDE estimates or Euler existence.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.182.0` with provenance `campaigns/P251-periodic-joint-continuum/adjudication.yaml`.

- `proposals/P251-cosserat-from-vortex-euler/attempts/0241/joint-residual.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0241/verify_joint.py`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0241/joint-third.stdout`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0234/current-boundary-target.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0242/review.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0243/hybrid-acoustic-repair.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0243/verify_controls.py`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0243/first.stdout`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0244/review.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0236/periodic-core-and-density.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0245/README.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0246/joint-current-bridge.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0246/verify_current_bridge.py`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0246/first.stdout`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0247/review.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0249/README.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0249/repaired-pytest.stdout`
- `src/substrate_framework/euler_joint.py`
- `tests/test_euler_joint.py`
