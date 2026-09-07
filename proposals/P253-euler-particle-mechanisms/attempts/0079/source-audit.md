# Primary-source and dependency audit

Full papers remain in `/tmp/primary-source-cache`; this attempt records only
versioned receipts and exact use boundaries.

| Source | Receipt and exact location used | Scope retained in 0079 |
| --- | --- | --- |
| Gallay--Smets, *Spectral stability of inviscid columnar vortices*, arXiv:1805.05064v3 | `/tmp/primary-source-cache/P253-0040/1805.05064.pdf`, SHA-256 `081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`. Proposition 2.1 and (2.2)--(2.5), PDF pp. 9--10, give the Fourier-reduced column transport/compact structure. Equations (2.17)--(2.18), PDF pp. 11--12, give the axisymmetric `m=0`, `k!=0` generalized energy pencil; Remark 2.5, PDF p. 12, records simple imaginary Kelvin eigenvalues accumulating at zero. | The paper does not provide a curved Cao-ring doublet, a distinct-sector frequency crossing, the compact Cao DA/KKS compression, or autonomous Euler controls. Its full-half-line vorticity space is not substituted for the compact-support Cao DA space. Equations (3)--(13) in the derivation are calculations from the stated Sturm pencil; the massive torus transfer is not attributed to this paper. |
| Cao--Lai--Qin--Zhan--Zou, *Uniqueness and stability of steady vortex rings for 3D incompressible Euler equation*, arXiv:2206.10165v2 | `/tmp/primary-source-cache/P253-0040/2206.10165.pdf`, SHA-256 `6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`. Equations (1.3)--(1.6), PDF pp. 1--2, specify the translating axisymmetric no-swirl ring. Proposition 1.4, PDF pp. 4--5, gives an existence family at fixed circulation; Theorem 1.6 and Theorem 3.1, PDF pp. 5 and 17--18, give pointwise small-parameter uniqueness up to axial translation. Proposition 3.2 and (3.8), PDF pp. 17--21, concern the scaled steady profile; Lemma 3.8, PDF p. 25, identifies the planar Lane--Emden translation kernel; Lemma A.2 and (A.5)--(A.8), PDF pp. 49--50, concern the outer free boundary. | These results identify the physical carrier and its axial symmetry. Pointwise existence plus uniqueness does not by itself prove continuity or local coverage in geometric `delta`; that is retained as an explicit IVT dependency. The source does not contain a nonaxisymmetric full Euler spectral branch, graph-Riesz derivatives, a two-frequency resonance, the constrained interaction functionals `G_12,G_3`, a two-ring relative equilibrium, or a Floquet three-wave solution. |
| P253/0074 with completed independent P253/0078 review | `0074/derivation.md` SHA-256 `ab7695b026e3008a05869a9887cb114c06c6d736e5c6fc0ea95fd8459c50dace`; `0074/result.yaml` SHA-256 `b2ca7ddeacf3ea88e0f1e560ce03c756afa535e66788e4e712c6884c6430a09f`; final correction receipt SHA-256 `1681c8f31b0fca5938e8d2f89d146a74a4cf1c902d20a44dafc7b1ca044915ad`. Independent `0078/review.md` SHA-256 `676fd487d310690cfe89246a5f0a8ba71c676c39c02aa69b4dd2a393a4b2c95e`; `0078/verdicts.yaml` SHA-256 `33725bbee5e149ca152146583614a8b8cc02bcad6ec49ed9e24ea77a06c2fe3f`. | All five scoped 0074 units are independently established: full fiber/Hodge, common action-flow graph, fixed-fiber essential gap, all-nonzero-streamline isolation, and transfer of two distinct-frequency simple positive-Krein modes. The review expressly leaves absolute physical Sturm-to-KKS conversion, frequency matching, and noncommuting physical couplings open. It does not establish exact-carrier continuity in geometric `delta`. |
| P253/0066 with independent P253/0073 review | `0066/derivation.md` SHA-256 `7bebcc7c51c93c6c28d5736730b9c1a07deb155c440942433a8ac55f54be1507`; `0066/result.yaml` SHA-256 `ceedfa810234491fa368f10310dacf35f956637a5a88042b8337d8983acad4d1`; independent `0073/review.md` SHA-256 `bed4e6e2bd2da7ce070764ca9812673cb71d3a61b395ccb7b902514ea523b76f`; `0073/verdicts.yaml` SHA-256 `ffc2f64bf8fea25b4516b53c2f3feacbf0599c80bb946a4ee1223fe88385bfda`. Review section 7 retains the limiting compact-column `m=0` generalized Sturm--Liouville modes and sectoral positive energy. | Supplies only the fixed-`k` asymptotic `sigma_J(k)=k L_Phi/(pi J)+O(k/J^2)`, with `L_Phi=integral_0^a sqrt(Phi)dr`. It does not supply whole-column isolation at that rung, absolute KKS normalization, or a joint high-`J`/thin-ring transfer with `J=Theta(delta^(-1))`; none is inferred in Route E. |
| P253/0051 conditional Schwinger--Hopf construction | `0051/construction.md` SHA-256 `841ca0caf6516e5955b2de89630bcf05402cf54473773c70b679308a23cd0103`; `0051/result.yaml` SHA-256 `f47398cb566bd3d3b4edc6952b6f5682ae4ee8c1f3b4c2161fedc73471948931`. | Supplies the finite-dimensional classical Schwinger--Hopf algebra only after a physical equal-frequency positive doublet and physical mixing maps exist. It does not supply those premises. |

## Derived here, not imported

The following are 0079 calculations:

- the classification of axial `SO(2)` irreducibles as one-dimensional
  characters, the absence of a symmetry-enforced equality among radial
  copies, and the exclusion of cosine/sine and `+n/-n` reality partners as a
  second complex oscillator;
- the Hellmann--Feynman identity (4), strict same-`k` frequency ordering (6),
  analytic matching ratio and rational-ray `C^0` IVT bracket (8)--(10i), and
  separation of the currently uncovered fixed-`n`, `k->0` transfer in
  (11a)--(11c); the IVT conclusion explicitly remains conditional on a
  continuous exact-carrier path covering the geometric-`delta` bracket;
- the exact first Euler derivative (14), its toroidal selection rule, and the
  continuous constrained response-functional definitions (17a), (18), with
  no expanded `ad/ad^*` representative inferred;
- the complex-linear/anti-complex-linear decomposition and Pauli-axis test;
- the conditional normalized high-harmonic response scale, its gate-time
  criterion, and the fixed-harmonic/high-radial-index alternative;
- the outbound/return Volterra estimates (26a)--(28b), which repair the
  preregistered one-way leakage expression and bound amplification rather
  than unsupported absolute signed drift; and
- the positive/negative-frequency projection requirement (31) for the
  three-wave route.

No nonvanishing of `G_12` or `G_3`, no exact Cao frequency crossing, and no
autonomous gate history is sourced or inferred. The symbolic oracle checks
only the finite-dimensional algebra and corrected Volterra-bound formula; it
does not validate any of those physical constructions.
