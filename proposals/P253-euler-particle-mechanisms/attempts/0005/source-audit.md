# Primary-source receipts for P253/0005

Full PDFs were cached only under `/tmp/primary-source-cache/P253-0005`; no downloaded paper
is stored in the campaign tree. The equations in `derivation.md` are rederived with one frozen
sign convention rather than copied across incompatible conventions.

## Euler dual, Hodge variable, and Hamiltonian leaf

- Boris Khesin, Daniel Peralta-Salas, and Cheng Yang, *The helicity uniqueness conjecture in
  3d hydrodynamics*, Transactions of the AMS 375 (2022), author offprint retrieved 2026-09-06
  from <https://www.math.toronto.edu/khesin/papers/helicityTAMS.pdf>.
- Cached PDF SHA-256:
  `a31b88353ec53c1a1d64a140c6e9053ae5138508f5d4255ac83d1e3d258aaafc`.
- Exact locations used: section 4.2, equations (7)-(9), journal pp. 920-922, for the
  `Omega^1/dOmega^0` dual, curl inverse, and coadjoint action; section 5.1, equation (10),
  journal pp. 922-923, for `partial_t[alpha]+L_v[alpha]=0`, kinetic Hamiltonian, and coadjoint
  orbits as symplectic leaves. Its theorem on helicity uniqueness is not needed for the action.
- Scope after the 0011 review: the dual/orbit equation supplies the complete momentum class.
  Decomposing `u^flat=delta G_2B+h` and harmonically projecting that equation gives
  `h_dot=-P_H(i_uB)`. The cited source is not read as saying that fixed-metric coefficients of
  `h` are individually conserved; Kelvin circulation applies to material loops and constrains
  the complete momentum representative.

## Clebsch action, local coverage, and helicity obstruction

- R. Jackiw, *A Particle Field Theorist's Lectures on Supersymmetric, Non-Abelian Fluid
  Mechanics and d-Branes*, arXiv:physics/0010042v1 (16 October 2000),
  <https://arxiv.org/abs/physics/0010042v1>.
- Cached PDF SHA-256:
  `3f725f9bff8795abbe1a2694cfe3873fee3c44e496893721ba1ef8ec761abca9`.
- Exact locations used: section 2.4, pp. 16-18, equations (48)-(56), especially
  `v=grad(theta)+alpha grad(beta)` and
  `L=-integral rho(theta_dot+alpha beta_dot)-H`; section 2.5, pp. 18-20, equations (58)-(71),
  for local Darboux coverage, canonical relabeling, and the total-derivative helicity density.
  The source itself explains that nonzero integrated helicity forces singular/boundary Clebsch
  data. After the 0011 review, 0005 imports this only as a local regular canonical realization:
  `H^1=0`/decay removes the harmonic-row obstruction, while global completeness still requires
  proved regular chart coverage or multiple/enlarged labels. It is not evidence that the
  displayed three bulk labels vary the harmonic row on an arbitrary nontrivial `H^1` manifold.

## KKS orbit realization and what prequantization actually establishes

- Francois Gay-Balmaz and Cornelia Vizman, *Vortex sheets in ideal 3D fluids, coadjoint
  orbits, and characters*, arXiv:1909.12485v2 (8 June 2020),
  <https://arxiv.org/abs/1909.12485v2>.
- Cached PDF SHA-256:
  `40032633529626eec1b04a9eb49c78c5352e8549bc8da4ffa8e4a5fc66600218`.
- Exact locations used: Proposition 2.4 and Theorem 2.8, pp. 6-7, for the pullback of the KKS
  form by an injective equivariant fluid momentum map; section 4, equation (46), Lemmas 4.1-4.2
  and Theorem 4.3, pp. 16-20, for a genuinely fluid coadjoint-orbit prequantization condition,
  bundle, connection, and curvature; section 5, pp. 20-24, for the separate character and
  polarization-group step. This vortex-sheet theorem is scope evidence for the geometric
  mechanism, not a smooth-carrier theorem imported into the compact swirl.

- William Gordon Ritter, *Geometric Quantization*, arXiv:math-ph/0208008v3
  (4 September 2002), <https://arxiv.org/abs/math-ph/0208008v3>.
- Cached PDF SHA-256:
  `327d7601ba7d8e59e954850afcbec5ddbd6d5781d07312278b33a0a90f52477f`.
- Exact location used: Definition 1 and the following discussion on PDF pp. 5-6. It states the
  integrality condition `[omega/(2 pi hbar)] in H^2(M,Z)`, the prequantum bundle curvature,
  and explicitly explains that the prequantum Hilbert space is too large and requires a
  polarization before the construction can be called quantization.

## Internal exact inputs

- `attempts/0001/material-balances.md`, equations (6)-(9), supplies the smooth compact radial
  swirl and its exact finite-energy ambient pressure. 0005 independently integrates its energy,
  angular momentum, and rotation-orbit KKS form.
- `attempts/0007/particle-facing-calculus.md`, section 2, supplies the centered no-swirl
  pointwise-helicity, finite-tag angular-momentum, stabilizer, and `SO(3)/SO(2)` topology
  calculation. It is used only for the narrow no-swirl route test.
- `attempts/0003` remains active attempt evidence for the Slobodeanu stationary map and its
  failed direct time-dependent transfer. It is not accepted quantum canon.

## P253/0011 correction receipt

The independent equation-level review is
`proposals/P253-euler-particle-mechanisms/attempts/0011/review.md`. It audited the prior
`derivation.md`, `source-audit.md`, and `result.yaml` hashes recorded there and preserved the
Euler--Poincare/orbit action, Hodge closure, local Clebsch KKS pullback, swirl KKS sphere, and
conditional prequantum findings. Its two required repairs are now applied:

- general `H^1` harmonic/circulation evolution is explicitly retained through representation A
  and the full Euler--Poincare equation, rather than inferred from an unvaried overlap cocycle;
- radial cutoff rotations are tangent/path lifts, while the global homogeneous sphere uses the
  ordinary rigid `SO(3)` action in the larger volume-preserving group.

No new external source body was opened for this correction. The six recorded symbolic checks
cover unchanged angular integrals, KKS area/chart data, conditional phase, and scaling; they were
not rerun because neither correction changes any checked expression.
