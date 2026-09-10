# 0124-beacon-sharp - Lipschitz-Weyl lemma build (beacon)

Charter: shepherd RULING on 0123 (credible x18 estimate-side path).
G-a2 physics half rides on this build.

## Frozen design (before compute — charter condition)

Object: replace the 0121 link-(ii) bound dlam <= 2*s_w*||dF|| by the
SHARP Weyl evaluation max|eig(dQ_lin)| with dQ the EXACT 3x3
perturbation from the gauge-projected linearized step (S linear in F;
no operator-norm chain). Link (i) dF = 0.91 post-S1-gauge is exact
computation, reused.

- Q/dQ assembly: tested ga_pipeline leray/make_grid (same math as
  0117 feed_member.py), n3=64, L3=8.0, fitted state + fitted J step.
- Gate G1 (charter L_J<=1): max|eig(dQ_lin)| <= 0.9 AND margin vs
  lam_min(Q_fitted) reported (>=x5 required to HOLD).
- Gate G2 (transfer proved inside): quadratic remainder bounded on
  the explicit tube E = 2*||du_est||_oo (measured max|du|); remainder
  via sup|d2zeta| on mesh x E^2 x pipeline factors; total must keep
  margin >=x2, else FALLBACK (not silent pass).
- Stop rule: G1^G2 pass -> lemma HOLDS, Maxwell H-side unblocks.
  G1 fail -> S2 quad-sensitivity probe (ONE assembly test, bounded).
  G2 fail -> S2 probe, then no-go certificate. No mesh rebuilds, no
  solves in 0124 (feed-grid assemblies + dense 1451 eig only).
