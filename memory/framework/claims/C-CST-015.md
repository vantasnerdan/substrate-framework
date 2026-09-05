---
description: Accepted framework claim C-CST-015
author: framework-registry
created: '2026-09-05T12:50:00+00:00'
updated: '2026-09-05T12:50:00+00:00'
tags:
- substrate-framework
- accepted-claim
- C-CST-015
category: claims
confidence: established
status: active
---
# C-CST-015

## Statement
The accepted statement is reproduced exactly from the claim registry.

The actual periodic stationary Euler field
u=(cos Y+A cos Z,A sin Z,-sin Y), A=1/100, curl u=-u,
admits a finite correlated initial displacement preparation with
positive whole-law physical restoring response equal to its
complete inherited displacement energy through spatial order two.

Put psi=cos Y+A cos Z, alpha=cos Y-A cos Z and H=-Delta.
For n=3,5,7 set w_n=n²(n²+1)/(2(n²-1)²), S=sum n²/w_n,
c_n=2n/(w_n S). Keep c_5,c_7 fixed and solve the two exact
first-shell moments for c_1,c_3 in
f(psi)=c_1 psi+c_3 T_3(psi)+c_5 T_5(psi)+c_7 T_7(psi),
so <cos Y(alpha+f)>=<cos Z(alpha+f)>=0. The finite inverse
phi=(H-1)^(-1)(alpha+f) then exists without discarding a mode.
Its full velocity z_*=(-alpha/2+Hphi,-phi_Z,phi_Y) solves
Lz_*+F0=0 with the complete pressure, where
L=-P[(u.grad)·+(·.grad)u] and F0=P(u_Z cross u_Y).

For a common unit laboratory direction kappa and D perpendicular
to kappa, set T_D=-(D.grad)u, q_D=-kappa cross T_D,
d=kappa_Y D_Y-kappa_Z D_Z, a_u=kappa.u,
b=q_D-a_u D+d z_*, and c=a_u D+(u.D)kappa.
The complete first-shell projection Pi_-=(P-curl)/2 supplies
z_return=-Pi_- b+t Pi_- c. It is killed by the full Euler
generator, not only by its planar subsystem. Thus
z_D=d z_*+z_return is a stationary forced first velocity cell.

Let E=1+A² and R=||z_*,|wave|>1||²/5, computed by exact full
Fourier convolution from this field. The construction gives
0<R<13E/1280. Averaging one actual whole-field isotropic law with
common laboratory data yields independent energy/current coefficients

  h2=E t²/15-47E/240+R,
  R_D=E(8t+13)/120.

Choosing s=sqrt(100-960R/E), t=-(4+s)/8 gives
a=E(s-9)/120, E/240<a<E/120, h2=a and R_D=-a.
The amplitude is derived from the actual finite microgeometry,
not an empirical fitted coefficient or assigned fluid mass.

In the genuine finite-Bloch preparation eta_D(0)=D,
w_D(0)=curl_K(T_D/(-1))+ik P_K z_D and K=k kappa.
Its full constrained cotangent has zero microscopic mean.
Together with eta_V(0)=0,w_V(0)=V, the complete initial phase
is rho J and the leading physical V kinetic mass is rho.
Circulation data are explicitly prepared, not falsely identified
with the constant-D Kelvin leaf. Real conjugate bands and whole
translations, rotations, mirrors and time reversal retain the
physical normalization and transform the helicity sign together.

The actual physical observation X=D+integral(mean w)dt satisfies
X_D(T)=D-a k²T²D/2+O_T(k³), and its inherited initial energy is
rho a k²|D|²/2+O(k³). The Hamiltonian calculation includes the
actual subtraction ||a_u D||² and full pressure terms. X is not
silently replaced by mean material displacement.

This is a displacement-column theorem on fixed time windows.
Its initial phase does not fix the common-V history, observed
Wronskian or autonomous acoustic Hamiltonian. Physical angle/spin
closure, acoustic-time uniformity, stationary EPS joining and the
full coupled continuum remain distinct constructions.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: none. Assumptions: Smooth constant-density incompressible Euler on the complete periodic cell; full mean-preserving Bloch Leray and independently specified initial circulation data. The initial physical D,V and their real sideband normalization are fixed before averaging., Standard exact Fourier convolution, periodic integration by parts, smooth finite-time Euler/Lin dependence and Haar integration are declared imports. The finite polynomial inverse, full-pressure energy and actual mean-current match are constructed here., The whole-field probability law remains positive; signed microscopic initial coefficients transform with the common laboratory direction and displacement. Energy averages use the actual quadratic phase data, not signed probability., Analytic identities and exact finite algebra are the oracle; the independent review recomputes the load-bearing shell constants in Cartesian phases. No Fourier cutoff, sampled stability window or fitted spectrum is used.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.180.0` with provenance `campaigns/P251-positive-displacement-preparation/adjudication.yaml`.

- `proposals/P251-cosserat-from-vortex-euler/attempts/0188/positive-correlated-preparation.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0188/initial-phase-energy.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0188/verify_chebyshev.py`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0188/chebyshev-first.stdout`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0188/velocity-lift-diagnosis.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0188/receipt.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0197/review.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0197/verify_review.py`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0197/first.stdout`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0199/README.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0199/repaired-pytest.stdout`
- `src/substrate_framework/euler_displacement_preparation.py`
- `tests/test_euler_displacement_preparation.py`
