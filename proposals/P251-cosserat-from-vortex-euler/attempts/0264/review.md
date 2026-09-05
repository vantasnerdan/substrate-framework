# Independent source audit — 0264

Reviewer: herdr optical-review pane w3:p3 (fresh, non-author review process).  This is one bounded audit of the frozen `0262/local-normalizer.md` and `0262/exact-gram.md`; `0262/current-transfer.md`, full-current claims, and compact-geometry/density claims are outside scope.

## Pinned inputs

* `local-normalizer.md` SHA-256: `13a00561e486f7d305bd392bd69eb91a3d4d23d8d6471c07f2f9d1b2bef19729`.
* `exact-gram.md` SHA-256: `dcde8a0fe47aabd145909b0154fdebeff248a7706278531ac6132fb60152517d`.
* additive API `src/substrate_framework/euler_joint_normalizer.py` SHA-256: `17d27b25f63b8c250c35fa6663e5118381ffec008f0cfd0ecf7ae922e27583a2`.
* direct test `tests/test_euler_joint_normalizer.py` SHA-256: `2971ba90eaa4057befbe8d6408639a521fc898562d760f2c553f72ff4a34c824`; focused run: `8 passed` (exit 0).

## Findings

The local zero-wave construction is substantively correct under its stated hypotheses.  The Leray identity gives
`curl P(eta x omega)=curl(eta x omega)=lambda[u,eta]` for a divergence-free compact source supported where `omega=lambda u`; the signs in `Omega=-rho(v_xi,zeta)` and in the displayed Jacobi/Hessian form then agree.  The two helicity polarizations have opposite definite H signatures and opposite KKS signs, and the disjoint-profile kernel construction can impose finitely many homogeneous cross rows without destroying the positive P/E margins.  This establishes the local opposite-sign Kelvin blocks and the zero-wave finite normalizer, conditional on an actual patch/background satisfying the declared support and separation assumptions.

The common-null-frame algebra in `exact-gram.md` is exact: with `H=diag(I,-I,0,0)` and the p/q symplectic block, `V=e+ + e- + p`, `W=(e+-e-)/2`, `Q=q` gives `V*HW=I`, `V*Omega Q=I`; hence `Y=V+(Wh+Qo)/2` realizes arbitrary Hermitian `h` and skew-Hermitian `o`.  The new API is a faithful finite-dimensional implementation of this identity (including imaginary diagonal skew entries, shape/signature checks), and its direct eight-case test is green.  It does not construct the physical Euler frame, as its own contract correctly states.

## Load-bearing gap

The claimed carrier-uniform *finite-K* continuation is not yet a source-level construction.  `local-normalizer.md` writes the Bloch forms with `curl_K` and `P_K` and imposes zero-mean rows, while `exact-gram.md` then assumes uniformly bounded K-derivatives of the whitened auxiliary forms and baseline cross matrices as the carrier `N_c` grows.  No derivation is supplied that, on the same periodic background, (i) the Bloch Kelvin/Lin identity really remains `curl_K(v_K)=lambda[u,xi_K]` for the selected compact/Bloch generators, including all `iK` product terms and the pressure projector, (ii) the finite profile-constraint kernel can be chosen with a K- and carrier-uniform norm while preserving the required zero-mean rows, and (iii) the baseline/auxiliary cross forms and their parameter derivatives stay bounded when their carrier supports and remote pressure tails are present.  A local constant-curl patch or the finite algebra alone does not imply these uniform estimates.  Consequently the exact-in-K, carrier-uniform neighborhood is conditional rather than established as written; this is a missing construction, not a counterexample to the zero-wave algebra.

## Minimum repair and supported verdict

Add one explicit periodic supplier lemma (and an exposing check) for the selected source inventory: state the Bloch differential convention and prove the `curl_K` Kelvin identity; give a fixed finite profile basis whose constraint matrix has a carrier/K-uniform right inverse and zero-mean image; and bound the full `P_K` pressure-tail contributions and all baseline/auxiliary Gram derivatives on a fixed `|K|<=K_0`, uniformly for sufficiently large `N_c`.  With those estimates, the contraction/implicit Gram argument in `exact-gram.md` supplies the advertised exact affine continuation.  Until then the strongest supported statement is: local opposite-sign physical Kelvin blocks and exact simultaneous finite Gram algebra are established under the explicit periodic-background/invariant-gap hypotheses, while the finite-K carrier-uniform physical normalizer remains open.  No conclusion about the excluded current/actual-observation join or compact geometry follows.
