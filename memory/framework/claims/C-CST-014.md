---
description: Accepted framework claim C-CST-014
author: framework-registry
created: '2026-09-05T12:22:00+00:00'
updated: '2026-09-05T12:22:00+00:00'
tags:
- substrate-framework
- accepted-claim
- C-CST-014
category: claims
confidence: established
status: active
---
# C-CST-014

## Statement
The accepted statement is reproduced exactly from the claim registry.

There exist two smooth compact-vorticity stationary Euler column
families with actual nodal m=1 modes, complete decaying radial
pressure tails and positive physical tilt action per axial period.
Euler amplitude and geometric scaling give a common frequency
nu>0, nonzero carrier slopes v_r and v_1² != v_2². The full radial
matching construction supplies the mode, its first two carrier
derivatives and a fixed positive stationary material tag with
displacement and mechanical-spin rows through that order.

For positive probabilities w_1+w_2=1, prepare the same initial
physical scalar angle/rate (q0,r0) with signed amplitudes
A_r(p)=a_r+d_r(p-p_r). Set b_r=nu_r'', j_r=M_raw,r/3>0,
h_r=j_r A_r², J0=sum w_r h_r, J2=sum w_r h_r'', evaluated at
the reference carriers. Normalization and variance cancellation give

  w_1 a_1=-v_2²/(v_1²-v_2²),
  w_2 a_2= v_1²/(v_1²-v_2²).

For any declared finite preparation curvature B_*>0, the equations

  sum w_r(a_r b_r+2d_r v_r)=B_*,
  sum w_r[h_r(v_r²+nu b_r)+2nu h_r' v_r]=J0 nu B_*

uniquely determine finite slopes d_1,d_2. Their determinant is
8nu w_1w_2v_1v_2(j_2a_2-j_1a_1), which is nonzero. Actual
mode-mass derivatives remain in h_r'; signed amplitudes are squared
in the inherited positive phase and physical energy, not treated
as probability weights.

The actual time-reversed circular modes are paired with tilt axis
n perpendicular to carrier t: Theta_s=n theta+s(t cross n)theta_t/nu.
Equal half-weights cancel the conjugate tilt/current and retain
both complete positive phase forms and Euler energies. For Haar
frames the observed vector is 3 E[n theta], whereas physical energy
is averaged without that factor. Thus j=M_raw/3 and Delta=Delta_raw/3,
and T(K)=3 E[n n^T(t.K)²]=(2|K|² I-K K^T)/5.

The resulting actual observed optical action has, through spatial
degree two on every fixed time window,

  M(K)=J0 I+J2 T(K)/2,
  K_opt(K)=nu² M(K)+J0 nu B_* T(K),
  Omega_opt²(K)=nu² I+nu B_* T(K).

It is stationary at this order and its energy equals the inherited
physical Euler energy. Gradient mass is retained. Its transverse
and longitudinal squared-frequency curvatures are respectively
2nu B_*/5 and nu B_*/5, both positive. Positive tags with equal
Delta supply the complete physical rows G=Delta Phi and S=Delta Phi_t
through the same order, retaining the nonzero overlap Delta/J0.

For any positive bulk curvature coefficients C_T,C_L, put
d0=4C_T+3C_L, c_s=3C_T C_L/d0, c_a=4C_T²/d0,
c_tr=C_L/2-c_s. The canonical curvature energy
c_tr(tr grad Phi)²+c_s|sym grad Phi|²+c_a|skew grad Phi|²
is pointwise positive and has the prescribed bulk symbols.
Its null-Lagrangian change of representative retains boundary flux.

This is a positive mixture of whole column realizations, not their
superposition as a stationary field. Its finite-action mode scope
is per axial period with full radial tails and declared real Bloch
bands/phase densities. It does not assert a finite-total-energy
monochromatic R3 mode, same-cell EPS/acoustic embedding, acoustic-time
uniformity, all-K invariant microscopic plane or a unique intrinsic
modulus independent of the declared preparation.

## Status Axes
The four governance axes remain independent.

Verification is `symbolic_verified`; review is `accepted`; compatibility is `compatible_extension`; epistemic status is `active`.

## Dependency and Import Closure
The registry records the accepted closure and declared non-claim inputs.

Dependencies: none. Assumptions: Smooth constant-density incompressible Euler columns with the complete pressure and coadjoint phase; selected near-Rankine nodal branches are smoothed after fixing a simple noncritical carrier. No gap to arbitrary exterior vorticity is assumed., Standard Bessel/Sturm theory, smooth ODE parameter dependence, implicit-function theory, Euler/Lin variation, Fourier analysis and Haar integration are declared mathematical imports. The actual supplier, tag controls and correlated energy construction are proved in0181/0187 with0191's append-only geometric correction., The positive whole-realization law, real paired carrier bands, common initial physical data and carrier-dependent signed preparation are specified before response elimination. B_* is an explicit preparation input, not an empirical fitted force., The API is conditional on actual-source frequency and phase-mass jets; arbitrary numeric inputs do not prove Euler existence or material observability. Analytic construction is the primary oracle, with exact symbolic and independent direct differentiation tests exposing false algebraic closure.. Comparators: none.

## Provenance and Evidence
The accepted release and immutable campaign evidence are the authoritative pointers.

Accepted in `v0.179.0` with provenance `campaigns/P251-correlated-optical-response/adjudication.yaml`.

- `proposals/P251-cosserat-from-vortex-euler/attempts/0181/nodal-domain-mode.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0181/receipt.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0187/correlated-family-action.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0187/nodal-family-realization.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0187/receipt.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0190/review.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0190/correction-check.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0191/transverse-family-correction.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0191/first.stdout`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0192/implementation-receipt.md`
- `proposals/P251-cosserat-from-vortex-euler/attempts/0192/repaired-pytest.stdout`
- `src/substrate_framework/euler_optical_response.py`
- `tests/test_euler_optical_response.py`
