# P253/0104 Unit-G bounded correction receipt

## Frozen correction scope

This is the single bounded P253/0104 correction to P253/0095 Unit G and the
combined nonlinear verdict.  Units A--F, the exact raw coadjoint/tag/Gauss
path, the fixed-`J_loc` observed-propagator theorem, the reusable API, and all
verifier predicates are unchanged.  The frozen README remains byte-identical.

Pre-correction SHA-256 values were:

```text
ac00b366e829cfa89a3552a28b94a0801dde0e155b5dd5198c962f24115b0859  derivation.md
fadc47dc75af852297010b899c890ac7664356951d80d684ead21f6cc34509ce  source-audit.md
51c276ddbe798f0a65e33caccf34751da21b459ee7cd46156195f47c4e1152be  result.yaml
28e259ee7d29c4bc1e566b8a8879048eddb01d5427c3a4280e62c5bd201d1685  validation.md
7df7d9fe8848ac8d4e683c9df4e24a5d71a807bd17bfa1aa167d83b975b0b4ee  author-completion-receipt.md
b305d2ea6ded594d24a20384789f03c86d282d175a4cf3c749a52fbc53280211  artifact-hashes.sha256
```

## Correction and strongest retained claim

Option (b) was applied.  On the smooth compact packet core, the declared
linear moment/slice functionals are defined, have axisymmetric coefficients,
and vanish on the nonzero toroidal character.  Their finite-rank microlocal
rows do not affect Units A--F.  No continuous extension or finite-row
submersion on the completed weighted input space is inferred.

The path (36bh) remains exact for the coadjoint/Kelvin-circulation data, tag
distribution, charge, divergence constraints, difference Gauss law, and zero
monopole.  It is now named an exact coadjoint/tag/Gauss path, not an exact path
in every frozen conserved row.  The only Gauss-reconstructed retained finite
row is charge/zero monopole, which is automatic; no electromagnetic momentum
row is silently retained.

The frozen README also fixes total momentum.  Its fluid term
`rho_m integral u` is not continuous on the declared `sigma<3/2` weighted
space for a generic `O(r^-3)` Hodge tail, which need not be `L1`.  Removing
that row would enlarge the leaf and cannot prove a conclusion on the frozen
smaller leaf.  For the remaining declared impulse/slice/carrier rows the
required witness matrix is now displayed as (36bh.2), but no carrier-specific
compact witnesses or nonzero determinant are available.  Density alone was
deleted as a purported proof.  The next exact construction is either a
relative-momentum domain plus the full electromagnetic witness row, or a
proof that total momentum is automatic/redundant for this path followed by an
explicit determinant for the remaining rows.

Therefore Unit G and the combined fixed-member nonlinear uniform-Lipschitz
refutation are `blocked` at this missing construction.  The Units-A--F theorem
survives: for each fixed sufficiently thin carrier and signed small nonzero
charge satisfying the displayed gap inequality, the exact returned phase and
continuity-persistent real hyperbolic pair give the fixed observed-propagator
essential-norm lower bound.  Equation (31), `det M=1+O(g^4)`, is used only for
persistence; exact reciprocity is asserted for the uncharged pair, not for
the nonzero-charge pair.  The charged-column continuation likewise uses a
conjugate paired leading off-diagonal entry and leading exponents
`+/-delta|b|+o(delta)`; it makes no exact charged-monodromy reciprocity claim
without an invariant symplectic two-dimensional reduction.

## Post-correction hashes

```text
b193c60f5e43e993d35673f0c5129beaa09a32c6035915915cdc98ebcd696c48  derivation.md
ac30e355638c86434653291b9882e39896e6cc0a218d869c201d26e600723ecd  source-audit.md
889cbe25ad2731b64fed259cd10053353a8767594332184d7a08698e427b8ec3  result.yaml
f284bdab343f297cce0ca9d183f982ec0de8b5b8917c33d9fdf7925b70e6c004  validation.md
e3ad2e35a0258e7bc89291faff8eb129813e8fa5a1415dae045a88fade287c85  author-completion-receipt.md
```

The frozen README remains
`9cbc9ef9773d102b6ee1c52c6af9af8345277cb23c020e1a5e76af54674290e8`.
The verifier remains
`ecac29802577e5f1d1e29eaaa874ab69f3909689df4136c0e0859e62a21bd66a`.
No verifier or API predicate bytes changed, so no unchanged oracle or focused
suite was rerun.  The final content-addressed list is refreshed separately in
`artifact-hashes.sha256` after this receipt.
