# Primary-source and dependency audit

Full papers remain in `/tmp/primary-source-cache`; no downloaded paper is
copied into this attempt.

## Primary sources

### Euler coadjoint orbit and KKS convention

Boris Khesin, Daniel Peralta-Salas, and Cheng Yang, *The helicity uniqueness
conjecture in 3d hydrodynamics*, Transactions of the AMS 375 (2022), author
offprint cached at `/tmp/primary-source-cache/P253-0005/helicityTAMS.pdf`,
SHA-256
`a31b88353ec53c1a1d64a140c6e9053ae5138508f5d4255ac83d1e3d258aaafc`.
Section 4.2, equations (7)--(9), journal pp. 920--922, supplies the
`Omega^1/dOmega^0` dual, curl inverse, and coadjoint action; section 5.1,
equation (10), pp. 922--923, supplies the kinetic Hamiltonian and coadjoint
orbits as symplectic leaves.  The paper does not construct Cao modes,
normalize their physical KKS covectors, or evaluate their interaction.

Francois Gay-Balmaz and Cornelia Vizman, *Vortex sheets in ideal 3D fluids,
coadjoint orbits, and characters*, arXiv:1909.12485v2 (8 June 2020), cached
at `/tmp/primary-source-cache/P253-0005/1909.12485.pdf`, SHA-256
`40032633529626eec1b04a9eb49c78c5352e8549bc8da4ffa8e4a5fc66600218`.
Proposition 2.4 and Theorem 2.8, PDF pp. 6--7, give a fluid momentum-map KKS
pullback.  Its vortex-sheet orbit is not imported as a smooth Cao-ring mode
or response theorem.

### Column eigen-equations

Thierry Gallay and Didier Smets, *Spectral stability of inviscid columnar
vortices*, arXiv:1805.05064v3, cached at
`/tmp/primary-source-cache/P253-0040/1805.05064.pdf`, SHA-256
`081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`.
Equations (2.7)--(2.16), PDF pp. 10--11, reconstruct velocity and vorticity
components and the scalar eigen-equation; equations (2.17)--(2.18), PDF
pp. 11--12, give the axisymmetric Kelvin energy identity.  These equations
determine the source mode functions but do not state a third orbit-Hessian,
mixed interaction coefficient, curved-ring response, or autonomous control.

### Cao carrier

Daomin Cao, Shanfa Lai, Guolin Qin, Weicheng Zhan, and Changjun Zou,
*Uniqueness and stability of steady vortex rings for 3D incompressible Euler equation*,
arXiv:2206.10165v2, cached at
`/tmp/primary-source-cache/P253-0040/2206.10165.pdf`, SHA-256
`6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`.
Equations (1.3)--(1.6), PDF pp. 1--2, fix the axisymmetric no-swirl
convention `zeta=omega_theta/r`, hence `bar_omega=r zeta e_theta`.
Proposition 1.4, PDF pp. 4--5, supplies the thin steady family.  Proposition
3.2, PDF pp. 17--21, and Lemma A.2, PDF pp. 49--50, supply the local scaled
profile and regular compact boundary used by the reviewed predecessor.
The construction preceding equation (3.29) and equation (3.29), PDF
pp. 27--28, define the odd first core correction through the linearized
Lane--Emden equation; it is used as `U_1^odd`, not replaced by an unspecified
Taylor coefficient.
The source does not contain a full-Euler KKS mode covector, a response
matrix, a curvature-diagonal coefficient, or a gate construction.

## Reviewed attempt inputs

- P253/0074 author derivation SHA-256
  `ab7695b026e3008a05869a9887cb114c06c6d736e5c6fc0ea95fd8459c50dace`
  and result SHA-256
  `b2ca7ddeacf3ea88e0f1e560ce03c756afa535e66788e4e712c6884c6430a09f`;
  independent P253/0078 review SHA-256
  `676fd487d310690cfe89246a5f0a8ba71c676c39c02aa69b4dd2a393a4b2c95e`
  and verdicts SHA-256
  `33725bbee5e149ca152146583614a8b8cc02bcad6ec49ed9e24ea77a06c2fe3f`.
  The reviewed scope supplies two simple positive-Krein massive modes and
  the common graph/Riesz transfer, while expressly leaving absolute physical
  KKS normalization and interactions open.
- P253/0079 corrected derivation SHA-256
  `ad80c2804798c9144b33bb78d1b9438af0e6f23928edc970825179eb597ad038`
  and result SHA-256
  `8f9747fcba0dfa7698de066d63c8641b501880b9c6274169b2db2e8f5a355cde`;
  independent P253/0082 review SHA-256
  `75e55672e359fdf91b91017d6098e266eb605f4ac7e289132153e22a4f7cf334`
  and verdicts SHA-256
  `69eee8688f986e7ee0463d92d9a115f6e397a74f2bbc13370d8d0ba414e175a4`.
  They supply the fixed-slice derivative, harmonic selection, dual-functional
  definitions, and two-sided control algebra, not response nonvanishing.
- Corrected P253/0080 derivation SHA-256
  `d2823746db0ef12e9a9d3df334031de3e3ac45e3841e402b13bbbb5f905ae479`
  and result SHA-256
  `dc84da15d2ebe68cd5a3f9330b5dd9351f8ba3ae1b8715cdb4530cfc0045a5ec`;
  independent P253/0084 review SHA-256
  `3fb7d2a69489668f6168d7f57eb9b4e86e35621064dd5f55172ca24270279aac`
  and verdicts SHA-256
  `7793688a23156d1b015ed82505f5e48b88223e19e8a36563c2f57ef63f0a74c7`.
  Only the reviewed `QAQ` HSE, bordered inverses and finite-window scopes are
  available; none supplies a response coefficient.
- P253/0083 derivation SHA-256
  `0dc3a837a201a97fd0bd89d9316992884e4b82f4d4fed8027d30f8473f5a8003`
  and result SHA-256
  `635216a567b2bdf92a12eba4a97cd46244f79fc800a55692a8f9204d45b07534`,
  with bounded correction receipt SHA-256
  `8894fdca60b3b39daad8c0759720aed27145923d5beb0583648a89892aefe256`,
  passed independent P253/0089 review SHA-256
  `da2d96384e1de7c018880d67a35f33a084515a4a80b58756233446cc2d2e01da`
  with verdicts SHA-256
  `36afcb332e6c5629b411fc84d17085be5fe3eeba03124629d4a6c323625a4be8`.
  The reviewed scope supplies the connected fixed-`(kappa,R)` thin path,
  noncircular actual-integer-fiber envelope, and exact same-family crossing.
  The byte-frozen 0088 README predates that review and therefore retains its
  activation-time author-stage wording; this body records the later reviewed
  authority without altering the replayed contract bytes.

## Derived here and boundary

The physical Hessian-row identity, exact clock/sign, skew-Hermitian
degenerate compression, full-Hodge dual density, compact patch exponent,
`(0,0)` stabilizer refutation, first allowed curvature character, and the
two-off-diagonal matrix criterion are 0088 derivations.  Gallay--Smets
equations (1.13)--(1.14), the KKS cancellation of the fixed-`k` orbit inverse,
and the physical velocity adjoint, rederived rather than imported, give the
complete near-axis interaction density.  The explicit KKS density is a
continuous `X_col^3` covector; extension through the exterior collar is
independent on the physical velocity image.  Reviewed algebraic simplicity
then identifies it with the rank-one dual Riesz row, while the regular-axis,
interface, and decaying-`K1` adjoint Sturm problem supplies its global
representative.  The exact Frobenius recurrence and core identity
`A_1^sharp=rho_0 A_1` establish the nonzero axis datum without inferring it
from Krein positivity.  Its nonzero
`A_2 A_1^sharp(k_1/k_2-1)` coefficient, punctured-annulus compact DA inverse,
and rational-character row exclusion prove `gamma_12^col>0`.
Reviewed 0074/0078 common-graph/Riesz/Hodge convergence, applied only to this
explicit smooth seed, plus direct convergence of both terms of the physical
orbit Hessian prove the source-specific fixed-`k` remainder
`M_(12,N)=M_12^col+o(1)`.  This transfers nonvanishing to sufficiently thin
exact crossings supplied by the independently reviewed 0083/0089 path, but
not to a normalized high-`N` lower bound.  The fixed-`k`
massive curvature rows are independently rederived from P253/0074's exact
cylindrical operators; P253/0044 equations (13)--(23) are used only as the
fixed-integer geometric cross-check `W/2,s Omega'/2`, not as a uniform
`n=Theta(delta^-1)` expansion.  No source is read as
proving the full order-zero curvature symbol, `D_curv!=0`, an actual-Cao
response remainder, a uniform high-index response, autonomous histories, or
a gate.  No pointwise compact remainder is promoted to collective
compactness across shrinking Sturm gaps.

No production numerics, comparator fitting, or small soft quantity is used.
The exact symbolic oracle only checks algebraic signs and scaling, so the
small-ratio numerical protocol is not activated.
