# Primary-source and dependency audit

Full papers remain in `/tmp/primary-source-cache`; no paper body is copied
into this attempt.

| Source | Exact location and use | Boundary retained |
| --- | --- | --- |
| Cao--Lai--Qin--Zhan--Zou, *Uniqueness and stability of steady vortex rings for 3D incompressible Euler equation*, arXiv:2206.10165v2 | Cached `/tmp/primary-source-cache/P253-0040/2206.10165.pdf`, SHA-256 `6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`. Equations (1.3)--(1.6), PDF pp. 1--2, give the translating no-swirl equation; Proposition 3.2 and equation (3.8), pp. 17--21, give scaled `C^1` profile convergence; Lemma 3.8, PDF p. 25, states that the bounded kernel of `-Delta-pU_+^(p-1)` is exactly the two translations; Appendix A, Lemma A.2 and (A.5)--(A.8), pp. 49--50, give the actual single outer free boundary and nonzero normal derivative. | The source does not give a nonaxisymmetric Euler graph, fiberwise torus-to-column Hodge convergence, essential-spectrum converse, Riesz projector, or nonlinear rotating branch. The higher fixed-domain profile regularity, the reduction of the nonzero-`m` Euler zero equation to Lemma 3.8, and all operator-transfer results in 0074 are derived, not attributed to Cao. |
| Gallay--Smets, *Spectral stability of inviscid columnar vortices*, arXiv:1805.05064v3 | Cached `/tmp/primary-source-cache/P253-0040/1805.05064.pdf`, SHA-256 `081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`. Their Fourier/enstrophy space is equation (1.23), PDF p. 6. Proposition 2.1 and equations (2.2)--(2.5), pp. 9--10, give the column transport-plus-compact split and explicit nilpotent resolvent. Section 2.1, equations (2.17)--(2.18), pp. 11--12, gives the `m=0`, `k!=0` energy pencil; Remark 2.5 records imaginary Kelvin modes accumulating at zero. | Their space contains vorticity on the full radial half-line and its passive exterior endpoint. 0074 does not import that endpoint into the compact-support Cao DA closure. Their results do not supply the curved-ring common graph or a Cao projector. |
| P253/0057 independent joined review | `attempts/0057/review.md` and `verdicts.yaml` at their recorded post-correction hashes. | Retains exact pointwise Piola/shape and metric-Leray algebra while classifying `HJ2`, `GR`, global graph bounds and Riesz transfer as conditional. 0074 derives a different leading massive fiber theorem and does not import those antecedents. |
| P253/0066 plus P253/0073 bounded proof completion | `attempts/0066/derivation.md`, `result.yaml`, `validation.md`, `source-audit.md`, and `0073-bounded-proof-completion-receipt.md`. | Supplies the fixed-carrier boundary jet, foliation/period bounds, compact physical DA/KKS seed and regular-cell fixed-`n` essential inclusion. It expressly leaves the converse spectrum, full gap and massive Cao modes conditional; 0074 addresses those dependencies. |
| P253/0048 plus independent P253/0053 review | `0048/derivation.md` equations (13)--(17), and `0053/review.md`/`verdicts.yaml`. | Supplies the compact `m=0` generalized Sturm--Liouville Liouville length, simple frequencies tending to zero and their spacing. It does not isolate them in the full column or transfer them to Cao. |
| P253/0044 author calculation | `0044/derivation.md` equations (2)--(9), rederived in 0074 equations (45)--(48). | Supplies a cross-check for the `m=1` translation ground-state transform, `-2F(a)` sign, exterior `K_1` row and nonzero small-`k` bending frequency. It is not used as authority for the common graph or the finite-`m` isolation theorem. |

## Derived in 0074

The following are attempt calculations rather than source claims:

- the exact normalized `n_delta` toroidal fiber, scalar/div/curl rows and
  helical masses `(k_delta+delta)^2`, `(k_delta-delta)^2`, and `k_delta^2`;
- the exact half-density toroidal Newton kernel (17), its local `K_0` limit,
  remote-axis and far-field compact-source estimates, and the resulting
  source-to-global-Hodge/Leray comparison;
- equality of the fixed nonzero-`n` ambient DA closure with the divergence-free
  `H_0^3` support space, obtained by interior inversion and density without
  dividing by the vanishing edge coefficient;
- the `k`-independent column domain and uniform graph equivalence;
- the characteristic converse for the compact DA transport space and compact
  invariance of the Fredholm essential spectrum; and
- the volume action--angle/flow conjugacy eliminating the otherwise unbounded
  transverse principal remainder, the leading graph convergence, and its
  Riesz consequence;
- the large-`|m|` Neumann exclusion, finite-`m` zero/nonidentity argument, and
  the isolation and positive-Krein-sign transfer of two simple `m=0` modes.

The last item establishes two simple positive-Krein spectral modes. The
reviewed 0048/0053 input supplies their sign, not the absolute conversion from
its Sturm amplitude to the physical density/Fourier/real-complex KKS
normalization. That conversion, the two noncommuting physical deformation
couplings needed by the 0051 analyzer, and a nonlinear branch remain open.
