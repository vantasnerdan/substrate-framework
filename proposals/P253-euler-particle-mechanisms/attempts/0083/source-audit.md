# Primary-source and dependency audit

Full papers remain in `/tmp/primary-source-cache`; this attempt records only
versioned receipts and exact use boundaries.

## Primary source

Cao--Lai--Qin--Zhan--Zou, *Uniqueness and stability of steady vortex rings
for 3D incompressible Euler equation*, arXiv:2206.10165v2, cached at
`/tmp/primary-source-cache/P253-0040/2206.10165.pdf`, SHA-256
`6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`.

Exact locations used:

- equations (2.1)--(2.6), PDF pp. 8--9: axisymmetric operator, physical
  measure `dnu=r dr dz`, symmetric whole-space Green kernel, local logarithmic
  expansion, regular axis, and far decay;
- Proposition 1.4 and Remark 1.5, PDF p. 5: existence at fixed circulation and
  source scaling observation; neither states continuity of a fixed-mean-radius
  family;
- equations (3.1)--(3.4), PDF pp. 17--18: exact steady equation,
  circulation, concentration, and core-diameter hypotheses;
- Proposition 3.2 and equations (3.8)--(3.12), PDF pp. 18--20: local
  Lane--Emden scaling and `C^1` profile convergence;
- equations (3.14)--(3.18), PDF p. 22: refined center/scale parameterization;
- Proposition 3.10, equations (3.25)--(3.29), PDF pp. 25--27: the first
  `O(epsilon)` correction and its odd radial character;
- Proposition 3.11 and equations (3.30)--(3.33), PDF pp. 27--29: the
  `O(epsilon^2 |log epsilon|)` remainder after the odd correction;
- Corollary 3.12, equations (3.34)--(3.35), and exact auxiliary system (3.36),
  PDF p. 29: value-level local Pohozaev balance and the exact mass/scale pair;
- Proposition 3.13, PDF pp. 29--30:
  `|x_epsilon,1-r_epsilon^*|=O(epsilon^2)` and
  `|s_Cao,epsilon-s_epsilon^*|=O(epsilon^3|log epsilon|)`;
- Lemma A.2 and equations (A.5)--(A.8), PDF pp. 49--50: actual free-boundary
  normal graph and nonzero boundary derivative; and
- Appendix B, especially (B.13) and Lemma B.4, PDF pp. 56--57: the mass and
  local Pohozaev identities underlying Corollary 3.12.

The source states Proposition 3.13 for fixed `W`.  It does not state the
compact-`W` uniform version used here.  The 0083 derivation obtains that
version by observing that the proof is pointwise in `epsilon`, contains no
`W` derivative, and has coefficients uniformly bounded for `W` in a compact
positive interval and the core center in a compact annulus.  This is an
author derivation, not a quotation.

The source does not prove the exact vorticity-mean versus maximum-center row,
a connected fixed-`(kappa,R)` external-epsilon family, graph/Riesz spectral
continuity, a frequency crossing, a KKS-normalized response, or a charged
analyzer hierarchy.

## Attempt inputs and authority boundary

### Corrected P253/0080 with independent P253/0084 review

- `0080/derivation.md` SHA-256
  `d2823746db0ef12e9a9d3df334031de3e3ac45e3841e402b13bbbb5f905ae479`;
- `0080/result.yaml` SHA-256
  `dc84da15d2ebe68cd5a3f9330b5dd9351f8ba3ae1b8715cdb4530cfc0045a5ec`;
- `0080/source-audit.md` SHA-256
  `3920cc313cc54e70d4ceb80d66b2dcfb12d6c3f15a2a4e60be86a474569ac5df`;
- `0080/0084-bounded-correction-receipt.md` SHA-256
  `c107a5d2b7d841ef2decd925c43138436b9e8341dffc88cbacc0a8dfcd452e74`;
- `0084/review.md` SHA-256
  `3fb7d2a69489668f6168d7f57eb9b4e86e35621064dd5f55172ca24270279aac`;
- `0084/verdicts.yaml` SHA-256
  `7793688a23156d1b015ed82505f5e48b88223e19e8a36563c2f57ef63f0a74c7`.

Only the independently accepted circulation--mean-radius bordered
isomorphism and uniform rescaled inverse are used.  They supply local
continuity with `epsilon` external.  No 0080 sharp Schur sign, stability, or
charged analyzer conclusion is imported.  The 0083 continuation and coverage
argument remains author-stage until 0089 adjudicates it; the input inverse no
longer has a pending-review qualification.

### Reviewed P253/0074/0078

- `0074/derivation.md` SHA-256
  `ab7695b026e3008a05869a9887cb114c06c6d736e5c6fc0ea95fd8459c50dace`;
- `0074/result.yaml` SHA-256
  `b2ca7ddeacf3ea88e0f1e560ce03c756afa535e66788e4e712c6884c6430a09f`;
- `0078/review.md` SHA-256
  `676fd487d310690cfe89246a5f0a8ba71c676c39c02aa69b4dd2a393a4b2c95e`;
- `0078/verdicts.yaml` SHA-256
  `33725bbee5e149ca152146583614a8b8cc02bcad6ec49ed9e24ea77a06c2fe3f`.

The reviewed boundary supplies uniform massive-`k` common-graph/Riesz
convergence and two simple positive-Krein modes.  It does not supply absolute
KKS normalization, frequency matching, physical coupling, or carrier
continuity in geometric delta.

### P253/0079 with completed independent P253/0082 review

The post-0082-correction author hashes are:

- `0079/derivation.md`
  `ad80c2804798c9144b33bb78d1b9438af0e6f23928edc970825179eb597ad038`;
- `0079/result.yaml`
  `8f9747fcba0dfa7698de066d63c8641b501880b9c6274169b2db2e8f5a355cde`;
- `0079/source-audit.md`
  `5e445da63f4160dd095c11026aa7bb8cf5798f7f0a4b326b0dfc1cabe097a993`;
- bounded correction receipt
  `a877241cddcfbac47f98863a90f919e340656c48909b16e3f9c44477f938208f`.

Independent `0082/review.md` SHA-256 is
`75e55672e359fdf91b91017d6098e266eb605f4ac7e289132153e22a4f7cf334`;
`0082/verdicts.yaml` SHA-256 is
`69eee8688f986e7ee0463d92d9a115f6e397a74f2bbc13370d8d0ba414e175a4`.
The review accepts the bounded correction and establishes the rational-ray
endpoint construction and corrected continuous-dual interaction boundary at
their scoped levels, while leaving the exact Cao crossing and analyzer open.

## Derived in 0083

The following are not imported claims:

- simultaneous local-PDE and circulation normalization gives
  `C_s=(Lambda_p/kappa)^((p-1)/2)R^(-(p+1)/2)`;
- the refined radial/odd/remainder decomposition plus the exact mean row gives
  `R-x_epsilon,1=O(epsilon^2)` without identifying mean and center;
- compact-`W` uniformization of Cao Proposition 3.13 along the fixed-mean
  border;
- the auxiliary-pair calculation gives the refined Cao scale an
  `O(epsilon^3|log epsilon|)` error, while the source-supported safe
  free-boundary area estimate gives
  `s_area=C_s epsilon+O(epsilon^2)` without assuming cancellation of the
  entire first-order angular mean or differentiating a remainder;
- a fixed-mean seed plus uniform relative-epsilon 0080 charts, overlap, and
  local uniqueness construct one connected thin component on the exact
  fixed-compact-source positive-part space, conjugated by its explicit
  parameter-only `(R,epsilon)` dilation/centering isomorphism to one fixed
  `X_*` while the physical scale remains an output; the explicit Green/log
  residual modulus, uniform `C^(2,beta)` bound with `beta>alpha`, compact
  embedding, and uniform inverse prevent positive-epsilon chart termination;
  a preliminary `1/N` band and sequential compactness/common-domain
  resolvent contradiction, with every limit identified by normalized
  Lane--Emden uniqueness under the fixed rows, prove the noncircular
  path-uniform error envelope tends to zero on the actual integer fibers
  `n_1=LP,n_2=LQ`; the massive `k` compacts serve only as uniform containers;
  buffered endpoint
  choices on that component give exact geometric-delta coverage, with the condition
  `1/N=o(h_N)` after multiplication by `n=NP`;
- the finite subluminal ceiling and explicit hierarchy inequality; and
- the local integration-by-parts representation of the constrained response
  through `F_12`, with the adjoint `P_row^*` of the bounded primal finite-row
  projection; actual curl nonvanishing and absolute KKS normalization remain
  open.

No production numerics, empirical comparator, or fitted physical constant is
used.
