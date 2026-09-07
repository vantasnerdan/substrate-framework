# P253/0095 source and authority audit

## Primary sources

### Choi: Hill stability

Kyudong Choi, *Stability of Hill's spherical vortex*,
[arXiv:2011.06808v2](https://arxiv.org/abs/2011.06808), 24 January 2022.

- cached PDF: `/tmp/primary-source-cache/P253-0095/2011.06808.pdf`
- PDF SHA-256:
  `0b2fe63374837db46ac054e81440be94629352d49c6e25414885983c51726783`
- extracted text SHA-256:
  `f0d26b1941af7ebe405186cfb800a1e492aaf0327f5303e4240ec4755099d99d`
- exact locations: Theorems 1.1--1.2, PDF pp. 4--5; Remarks 1.3--1.6,
  PDF pp. 5--6; the paper's scope statement at PDF pp. 3--4.

Consumed scope: all-time stability modulo axial translation for the stated
nonnegative axisymmetric no-swirl relative-vorticity perturbations in the
`L1+L2+impulse` metric.  Theorems 1.1--1.2 do not include arbitrary
nonaxisymmetric vorticity, vortex stretching, a full three-dimensional
coadjoint Hessian, Maxwell fields, or a transported charge tag.  The source
itself describes most nonaxisymmetric and short-wave literature as linear,
approximate, or numerical; none is imported as a theorem here.

### Gavrilov: smooth compact swirl comparator

Andrey Gavrilov, *A steady Euler flow with compact support*,
[arXiv:1810.08020v1](https://arxiv.org/abs/1810.08020).

- cached PDF: `/tmp/primary-source-cache/P253-0010/1810.08020.pdf`
- PDF SHA-256:
  `fcaca85faa77e3876b11d16718037169fc026112dfd3e4248e03e963a0ebc3c9`
- exact locations: theorem on PDF p. 1; local field and flat pressure-cutoff
  construction in equations (6)--(7), PDF pp. 5--6.

Consumed scope: existence of a nonzero smooth compact finite-energy steady
Euler field with meridional flow and swirl.  The source provides no nonlinear
orbital-stability or coercive-Hessian theorem.

### Shvydkoy--Vishik: BAS corroboration only

Roman Shvydkoy and Misha Vishik, *On Spectrum of the Linearized 3D Euler
Equation*, Dynamics of PDE 1 (2004), 49--63.

- cached PDF: `/tmp/primary-source-cache/P253-0039/Shvydkoy-Vishik-2004.pdf`
- PDF SHA-256:
  `9307efc6096565f74c5a099af25c834af09ba93b3e56ff34dd076ba289551592`
- exact locations: bicharacteristic-amplitude system (1.5), PDF p. 3;
  Theorems 1.1--1.2, PDF pp. 4--5.

Consumed scope: corroboration of the standard full-pressure Euler BAS
equation.  Its spectral theorems are on `T3`; no torus essential-spectrum
conclusion is imported for the Cao whole-space free-boundary problem.

### Hattori--Fukumoto: invariant thin-ring Hill reduction

Y. Hattori and Y. Fukumoto, *Short-wavelength stability analysis of thin
vortex rings*, Physics of Fluids **15** (2003), 3151--3163,
DOI `10.1063/1.1606446`.

- cached PDF:
  `/tmp/primary-source-cache/P253-supervisor/hattori-fukumoto-2003.pdf`;
- PDF SHA-256:
  `c6a35c44d55d26b8bb3e8fc9e55b29b638d0983b646b17eb522a4192b6444fb9`;
- exact locations: BAS system (2)--(4), PDF p. 3152; transformed two-component
  system (7) and Hill equations (8)--(9), PDF p. 3152; general leading-profile
  first harmonic (22), resonance (23), and growth exponent (24), PDF
  pp. 3156--3157; the warning after equation (20) that a nonreturned
  higher-order wavevector can make first-order growth saturate, PDF p. 3154;
  arbitrary first-order particle/covector fields, initial covector
  `(sin chi,0,cos chi)+O(epsilon)`, and their substitution into (9), Appendix
  C, PDF pp. 3161--3162.

Consumed scope: the exact algebraic reduction of the full-pressure BAS
system for a thin steady no-swirl vortex ring and its first-order invariant
for an arbitrary smooth leading circular-column profile.  The paper's Kelvin
and Gaussian rings are source examples, not identified with Cao's carrier.
Application to Cao uses the independently derived compact-interior profile
topology and exact first-integral phase below.  The source's approximate
nonreturn is not silently discarded: equations (36bb)--(36bb.6) replace its
local wavevector by `d(F(P)+ell vartheta)`, prove exact return, and tune the
first diagonal detuning while preserving the paired Melnikov coefficient.
No source theorem about the Cao free boundary, dynamically
accessible closure, whole-space semigroup, Maxwell coupling, or nonlinear
orbital stability is imported.

## Reviewed framework inputs

- P253/0068 as reviewed by 0071 supplies the exact transported current,
  Euler--Lorentz equation, Maxwell equations, joint energy/momentum balances,
  and local smooth constraint propagation.
- P253/0077 as reviewed by 0081 supplies the exact axisymmetric no-swirl
  traveling reduction
  `W=r^-1 grad P cross e_theta`, `B=rH e_theta`,
  `F=E+u cross B=-grad Phi-H grad P`, and the material tag `chi(P)`.
- Corrected P253/0080 as reviewed by 0084 supplies the finite-window
  subluminal charged Cao branch, including `E_g,B_g=O(g)` and fluid/tag
  correction `O(g^2)`.  It supplies no dynamical stability.
- Corrected P253/0085 as reviewed by 0086 supplies the infinite positive and
  negative joint Hessian sectors and `iR` Fredholm-essential inclusion.  Its
  explicit warning that Hessian indefiniteness does not imply spectral
  instability is retained.
- Corrected P253/0032 and 0039 as independently reviewed in 0038 supply the
  hyperbolic full-pressure accessible cocycle and whole-space semigroup lower
  bound for one **Gavrilov** carrier only.  That carrier is not the Cao branch.
  The result is consumed solely to refute coercivity of this concrete smooth
  swirl comparator.

## Derived in 0095

Equations (8)--(15) derive the coupled Euler--Maxwell--tag principal hierarchy
from the reviewed equations.  Equations (16)--(25) expose the exact
Stokes-streamfunction-normal ambient return, including variable speed,
cylindrical radius, both magnetic polarizations, and its leading DA
degeneracy.  Equations (26)--(31) derive the first nondegenerate toroidal DA
cocycle with the tag eliminated on the same leaf.  Equations (32)--(36)
evaluate the full-pressure returned thin-column matrix and its positive
metric on a compact regular cell where `Omega,zeta_col` have positive lower
bounds.  Equations (36ai)--(36av) independently eliminate the printed
Hattori--Fukumoto two-component system, retain every parenthesis in its Hill
potential, and derive the smooth-center physical exponent `15/256`; this
supersedes the quarantined `255/1024` raw matrix-entry combination.  Equations
(36ax)--(36ba) use Cao's value-level source estimates, the full `r^2/2` speed
row, Appendix-C compactness, and interior Schauder regularity to transfer that
coefficient on a fixed compact positive-core tube without asserting endpoint
differentiability.  Equations (36bb)--(36bg) construct the globally
single-valued exact-return integer-toroidal DA packet, constrained transverse-
Maxwell/longitudinal-Gauss split, fixed-time global-input-to-local-output
Egorov estimate, raw exact coadjoint/tag/Gauss path, and essential-norm lower
bound for the fixed observed propagator `J_loc S_g` from the weighted global
complete-state input to the local same-`H^s` observation.  That observed-
propagator bound is the established Units-A--F conclusion.  On the smooth
compact packet core, declared linear moment/slice functionals have
axisymmetric coefficients and vanish by the nonzero toroidal character; no
continuous extension or submersion on the completed weighted space is
inferred.  Equation (36bh) is an exact coadjoint/tag/Gauss path, but not an
exact fixed-conserved-row leaf path.  The frozen leaf fixes total momentum,
whereas `rho_m integral u` is not continuous on the declared topology for a
generic `O(|x|^-3)` Hodge tail, and no relative-momentum domain or explicit
carrier-specific finite-row witness matrix is supplied.  Therefore the
combined nonlinear Lipschitz refutation is blocked while the linear observed-
propagator theorem survives.  For each fixed sufficiently thin Cao member and
`0<|g|<g_0(delta)` satisfying `C_mon(delta)g^2<c_HF delta`, continuity
preserves a real hyperbolic pair and expanding multiplier; equation (31) does
not prove exact reciprocity at nonzero `g`.  No delta-uniform charge wedge is
derived, and no conclusion is extended through the vanishing-vorticity
interface or a Maxwell-characteristic coincidence.  Equations
(36bj)--(36bs.1) retain the
complementary charged regime: they derive the exact fixed-delta response row,
the charged-column Maxwell/radial ODEs and constrained slow block, and the
first- and second-harmonic ordered-integral tests.  A delta-uniform charged
column/free-boundary construction remains the precise missing input.

Equations (37)--(38) are direct Fourier calculus, not sourced stability
claims.  Equation (39) is a reviewed framework theorem applied to the exact
Gavrilov carrier.  Equation (41) is the distributional coadjoint tangent for
the Hill patch and explicitly retains its interior bulk mode.

No production numerics, small spectral edge, or fitted comparator is used.
Reusable exact algebra is exposed by
`src/substrate_framework/euler_p2_principal.py` and its focused tests in
`tests/test_euler_p2_principal.py`.  The attempt verifier imports that module
and independently assembles the source-specific symbolic predicates; the
attempt-local `p2_principal.py` is compatibility re-export only.
The original 5/5 receipt remains the pre-curvature baseline.  The Luna-low
9/9 receipt and its `255/1024` predicate are preserved but superseded by the
recovery receipt.  The Recovery-7 Sol-High exact symbolic verifier checks the algebraic
hinges named in validation, including the two independent `15/256` routes,
the detuned paired Floquet generator, and the charged-column slow-block square; it
does not by itself prove the compact-interior topology, graph-domain packet
remainder, the missing fixed-row leaf correction, Hill Hessian, or nonlinear
LP2.
