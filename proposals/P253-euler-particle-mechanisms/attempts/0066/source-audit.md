# Primary-source and dependency audit

Full papers remain in `/tmp/primary-source-cache`; no paper body is copied
into this attempt.

| Source | Exact location and use | Boundary retained |
| --- | --- | --- |
| Cao--Lai--Qin--Zhan--Zou, *Uniqueness and stability of steady vortex rings for 3D incompressible Euler equation*, arXiv:2206.10165v2 | Cached `/tmp/primary-source-cache/P253-0040/2206.10165.pdf`, SHA-256 `6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`. Equations (1.3)--(1.6), PDF pp. 1--2, give the no-swirl potential-vorticity and translating equation. Proposition 1.4, pp. 5--6, gives the polynomial compact ring. Proposition 3.2 and (3.8), PDF pp. 17--21, give only scaled `C^1` convergence to the radial ground state. Appendix A, Lemma A.2, equations (A.5)--(A.8), PDF pp. 49--50, parameterizes the actual outer free boundary and proves its uniformly nonzero outward derivative. Together with fixed-carrier Schauder bootstrapping of the exact semilinear PDE, these equations support the `O(d^2)` normal jet and `O(d^(p+1))` vorticity expansion on each fixed compact boundary. They do not support uniform-in-epsilon fixed-domain `C^2` convergence, unique-center/global-foliation, or the high-`n` coefficient ledger. | The source does not state a nonaxisymmetric Euler graph, essential-spectrum theorem, Kelvin mode, Riesz transfer, nonlinear rotating branch, or charged continuation. Its axisymmetric rearrangement stability is not used as three-dimensional coercivity. |
| Gallay--Smets, *Spectral stability of inviscid columnar vortices*, arXiv:1805.05064v3 | Cached `/tmp/primary-source-cache/P253-0040/1805.05064.pdf`, SHA-256 `081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`. Proposition 2.1 and equations (2.1)--(2.6), PDF pp. 9--10, identify the column multiplication/compact split in their stated space. Section 2.1, equations (2.17)--(2.18), PDF pp. 11--12, gives the axisymmetric `m=0`, `k!=0` pencil and Rayleigh energy identity; Remark 2.5 records the imaginary Kelvin modes. | Their ambient column space includes passive exterior vorticity and its `Omega(r)->0` endpoint; it is not imported into the compact Cao DA closure. Their theorem supplies neither curved-ring graph convergence nor two isolated full-column modes. The fixed-positive-mass Piola/resolvent transfer remains a conditional 0066 route. |
| P253/0057 independent review | `attempts/0057/review.md` and `verdicts.yaml`, at their current recorded hashes. | Its blocked common graph, all-sector `GR`, and small-`k` Riesz claims remain blocked; 0066 uses a different `k_delta->k_0>0` one-contour construction. |
| P253/0048 plus independent P253/0053 review | `0048/derivation.md` equations (13)--(17) derives the generalized compact-core Sturm--Liouville Weyl law with `L_Phi=integral_0^a sqrt(Phi)dr`; `0053/review.md` and `verdicts.yaml` retain that Weyl coefficient and spacing at their recorded scope. | This supplies the limiting `m=0` frequency accumulation only. It does not isolate a mode from the whole column or transfer it to Cao. |
| P253/0062 plus the 0063 bounded correction | `attempts/0062/derivation.md`, `result.yaml`, and `0063-bounded-correction-receipt.md`. | Imported only for the exact harmonic block/intertwiner/KKS convention and corrected `s=4`, `p>=6` collar theorem. `q_*`, `V_*`, and all-sector limiting absorption remain undefined. |

## Derived rather than imported

The following load-bearing steps are calculations in `derivation.md`, not
source claims:

- (A.8) plus the exact semilinear equation is bootstrapped to the two-sided
  `d^p` carrier jet;
- nested-ball `L^infinity -> W^(2,q) -> C^(1,gamma)` estimates, interpolation
  with Cao's `C^1` convergence, and Schauder applied to the complete rescaled
  equation yield fixed-domain `C^(2,alpha)` convergence for the sufficiently
  thin family;
  nondegenerate-center persistence, the off-center gradient bound and the
  regular outer boundary then give the fixed-carrier one-center foliation and
  positive finite period bounds used in (7);
- the tangent quotient and smooth-path Hodge homotopy use the `H^5`
  volume-preserving group and `B curl F=P_LF`; their nonlinear quotient
  descent remains open;
- fixed-`n` Hodge compactness uses common compact curl support, the exact
  zero integral/first-moment low-frequency control, global decaying
  Biot--Savart, and local Rellich; together with explicit, exactly divergence-free,
  tangent-polarized DA Weyl sequences prove regular-cell transport-band
  inclusion; positive-core `C_0` is order zero in the high-`I` frequency,
  the normalized `H^3` residual is `O(epsilon+N^-1)+o(1)`, and finite-row
  corrections from the same fixed-`n` smooth DA graph are `o(1)`; a graph
  left-regularizer contradiction places the points in the upper-semi-Fredholm
  and standard non-Fredholm essential spectra; the converse
  constrained direct-integral and boundary/interface theorem remains open;
- `n_delta=nearest(k_0/delta)` supplies a useful uniform scalar elliptic mass,
  but equations (34)--(37) are targets pending the full vector Piola/Hodge
  `D->X` ledger and entire-column isolation; and
- the exact column wave energy gives positive Krein sign only to the limiting
  `m=0` column modes; no sign is transferred to an unconstructed Cao mode.

No production numerics, filament spectrum, boxed eigenvalue, fitted action,
or charged-force cancellation is used. The small-ratio numerical protocol
therefore remains a frozen contingency rather than an executed verifier.
