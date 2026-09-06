# P253/0057 joined independent review of 0052 and 0054

## Frozen transaction and correction provenance

This transaction independently reviews worker-owned attempts `0052` and
`0054`. The reviewer authored or implemented neither target, its correction,
verifier, API, tests, nor receipts. The frozen review README has SHA-256
`99f0f47c10f6e04565d61d6eccba47c1ef05ee5a0f34e4ab9c9872babadabfb5`.
Central activation exited exactly `0`; the command, stdout, empty stderr, and
exit hashes are respectively
`5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f`,
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`,
and `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.
The accepted authority base is release `v0.183.0`.

The substantive pass found one joined scope defect and requested one bounded
correction. The final correction receipt has SHA-256
`4659cb522e3f46cd01099405bd12a7841aeb3cff2274e4a689f03bacc25cceeb`.
The correction check was limited to the changed claims and their direct source
and operator edges; no second substantive pass and no unchanged formula-oracle
rerun was performed.

The receipt's post hashes all match the worktree:

| Artifact | Frozen SHA-256 | Corrected SHA-256 |
|---|---|---|
| `0052/README.md` | `219334fcd5bbf40f69d691dc45596a82926b5f25e63836dc3a9eee5990ba6bd1` | `f534e6ee4b1dabb0ea3e0fdb9a0223139f5deb860571b1ee98229f06093d28ee` |
| `0052/derivation.md` | `d353a4ffb8ea1ce65b419f8152b3292add53671a3d5e94195685722ae95f1d5b` | `6ea25cd99b91686d0d7efada0bc756656b82b4ebe44f310848511ab63da23e76` |
| `0052/result.yaml` | `9bab0def7992099eace135783d77e4e2a84f96f55792089d0f51177912afc8a8` | `1d0781660745266a447bf8ecf2f2e67c7804bfb1572060ad05ce56ec7824cc7a` |
| `0052/validation.md` | `e7a9ee362443ceca42dbb631977ad87218fe6a3353ff9722079007f7810de509` | `054dd0f82cf17afc2b444d545677fbbdf5f6d1f4419cf35548bf4d55fd6eebbe` |
| `0052/source-audit.md` | `1dd942371cf4591ad991dfe9bff8ab67b9a64895ffc84f06442adc04dfd558c2` | `1da32efb5f36267c3e44717461fb178a3c3d3f1ea5c759f8db6f89275a6dace9` |
| `0054/README.md` | `a13c902a9798e48897aa787f3eecdaf331843a78d0ade783690d8ca476e96a8c` | `973729e59d2acbdf295eff2ddba3e95a9421bea08d886488bc44256361319e0a` |
| `0054/derivation.md` | `62e3d45cd8c45099c51b47467339a0824b3fdc1ddedcc371645258723575d760` | `6d4835ae109395c5e7c4be4d8d9fa285d634879e1be10ae383f350b464ca54ee` |
| `0054/result.yaml` | `f060d63e92cdf69d414a0a795dcbd35d4e5e17da6ae0557de6d5a7e4d28eebe6` | `ad6e7a19f519294d63a404c274eafe9fa3171a7dc4de65131012200cdeb6c947` |
| `0054/validation.md` | `e5ee164a75b4e75a13a48d2056c9a1f15d867c3a96ae7902881522a9ec88e829` | `0f648c06e2d5221248cb87d7cb399b5e93abe0e6502db374a5ec7e8025494543` |
| `0054/source-audit.md` | `4ba9323c8f2e16f3418d1809d1245b49e9d688f924f7a50d6b6b952e6c508f6f` | `9d5c7de1cf150438a3a610e214636855f0ea2d590340257af747c00c4dd89519` |

## Source applicability

The load-bearing primary source is Cao--Lai--Qin--Zhan--Zou,
arXiv:2206.10165v2, cached with SHA-256
`6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`.
Its exact meridional equation and Green kernel, limiting Lane--Emden kernel,
centered estimates, and uniqueness theorem apply at their printed scopes.
They do not provide a twice differentiable nonaxisymmetric Euler graph family.

In particular, source equation (3.36) defines auxiliary parameters. Corollary
3.12 gives the corresponding true parameter equations only up to
`O(epsilon^2 |log epsilon|)`, and Proposition 3.13 gives nonzero comparison
errors between the physical and auxiliary parameters. Theorem 1.6 says that
two already existing families satisfying its exact hypotheses agree up to
axial translation; it neither constructs an augmented Banach map nor proves
that map surjective or differentiable in `(q,q^2 log q)`.

Gallay--Smets arXiv:1805.05064v3, SHA-256
`081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`,
fixes the straight-column Fourier/Biot--Savart normalization but does not
supply the compact solid-torus DA graph inverse. The remaining existence
papers listed by 0054 are applicability inventory, not imports of HJ2.

## Unit A: corrected 0052 finite-core results

### Resonance and leading scales

The inversion

    sigma_n=A delta/(n+beta_0)
       +B delta/(n+beta_0)^2+O(delta/n^3+delta^3/n),
    tau_delta=c delta^2 L+d delta^2+o(delta^2)

gives

    n_*=A/(c delta L)-A d/(c^2 delta L^2)
       +(B/A-beta_0)+o(1/(delta L^2))+O(1).

The two distinct correction terms and the bounded integer choice are retained.
The leading coupling `O(delta^3 L)`, local spacing
`Theta(delta^3 L^2)`, and ratio `O(1/L)` follow at the stated sequence-level
scope. The exact Cao first and second cells and the elliptic-kernel coefficients
remain fixed algebraic atoms.

**Verdict: established.**

### Chart, Piola, and KKS normalization

Weighted enclosed area and normalized streamline time give the exact local
action form. Across unequal carrier volumes the map is correctly treated as
an orientation-preserving representation with the contravariant Piola factor
`J^{-1}D Phi`; it is not called a volume-preserving coadjoint motion. Only a
same-carrier displacement is dynamically accessible. The KKS sign follows
from `L(q,conjugate(q))=i nu Omega(q,conjugate(q))`: energy normalization gives
`Omega=-i/nu`, and multiplication by `sqrt(nu)` gives pairing `-i`.

The displayed identities do not establish the two-sided common closed DA
range or uniform generator graph equivalence. Those are correctly returned to
HJ2 after correction.

**Verdict: established at exact representation, finite-core energy, and KKS
normalization scope; blocked at common closed DA graph scope.**

### All-sector complement

The large-`|m|` triangular inverse, the translation-removed `m=1` scalar gap,
the positive fixed-`|m|>=2` form increment, and the formal `m=0` Schur scaling
are useful sectorwise calculations. They do not prove an `X_* -> D_*`
resolvent for the coupled nonnormal operator. Spectral separation cannot
replace resolvent bounds, adjoint-domain control, the low-mode kernel/cokernel
calculation, or uniform Hodge/KKS projection constants.

The correction introduces the independent hypothesis `GR` for exactly this
missing construction.

**Verdict: blocked by GR; the scalar and sectorwise calculations remain
established.**

## Unit B: actual 0054 HJ2 support

### Augmented carrier branch

Expansion of the exact rescaled equation gives the displayed Lane--Emden cells
and a formal `(q,q^2 log q)` candidate. The attempted augmented map never
specified an exact scale/speed row, Banach codomain, Fredholm index, parameter
Schur determinant, or the exact `epsilon-q-tau` restriction. Source (3.36) is
auxiliary, and source uniqueness cannot provide these missing objects.

The correction now states the exact missing construction and no longer claims
that the candidate is the actual differentiable Cao family.

**Verdict: blocked by the exact augmented-map and surjectivity construction;
the local cells are established.**

### Interface and common DA domain

For integer `p>=6`, the powers `d_+^p`, `d_+^(p-2)`, and the DA tangent's
`d_+^(p-1)` behavior support the limited conclusion that zero extension of the
base vorticity, its formal first two coefficient derivatives, and the DA
tangent produces no vortex sheet. They do not establish the traces of
`A_delta eta`, a frozen graph Sobolev interval, the reverse Hardy estimate,
closed common DA range, harmonic/zero-mode controls, or the composite
whole-space Hodge bounds.

**Verdict: limited no-sheet trace established; full trace/Hardy/common-domain
claim blocked by the named two-sided domain construction.**

### Whole-space shape and Leray identities

The Piola-conjugated Biot--Savart formula cancels the volume Jacobian exactly.
Its first and second shape derivatives retain the differences
`h(y)-h(y')` and therefore preserve the pointwise order-minus-one kernel. The
metric divergence, gradient, inverse-Laplacian, and Leray differentiation
formulas have the correct algebraic signs and include pressure rather than
discarding it. The exact circular-tube distance and phase
`exp(i ell vartheta)=exp(i ell delta t)` are also established.

These identities do not bound the global decaying/zero-moment composite Hodge
operator or prove graph-norm differentiability. The correction makes this
distinction explicit.

**Verdict: established at exact pointwise shape, toroidal geometry, and
metric-Leray algebra scope; operator/graph estimates blocked.**

### Matched remainder and two-index decay

The scalar matching integral has size `O(delta^3 L)`, but it does not prove
the front/side/overlap cross-sectional symbol seminorms needed for an
`H^s -> H^(s+1)` estimate. Those anisotropic estimates remain open.

Likewise,

    (lambda_a-lambda_b)^N <e_a,B e_b>
      =<e_a,ad_Lambda^N(B)e_b>

is exact, but one `Lambda` controls only its combined eigenvalue difference.
Changes in `m` and radial index can cancel, so this identity cannot yield the
formerly asserted product decay in both indices. Separate and mixed
commutators, or a directly proved two-dimensional summable bound with uniform
eigenbasis constants, are required.

**Verdict: the scalar scale and one-generator identity are established; the
symbol remainder is blocked, and the claimed one-generator-to-product-decay
inference is refuted by spectral-difference cancellation.**

### Contour and Riesz transfer

If HJ2 supplied the common graph remainder and GR supplied the full nonnormal
complement resolvent of size `C/r_delta`, then the displayed Neumann resolvent
identity and contour integration would give the advertised conditional
resolvent/Riesz estimates. Neither antecedent is constructed in 0054.

**Verdict: the conditional algebra is established; the actual Cao
norm-resolvent and Riesz transfer is blocked by HJ2 and GR.**

## Oracle scope

The frozen 0052 symbolic receipt has source SHA-256
`1caf2ac69fdd619b9ebcf01b98ea61d9e97328680f870e52f2c54f12be4d2884`
and stdout SHA-256
`ea2b14df27d533e1b4810f1e80e9a9fa7e7abf79404955461c51ffda7345c182`.
The frozen 0054 verifier has source SHA-256
`0f04bad95ee0954050490a372f16bea651c0e0f92838d9d6dfa36f26aabae83c`
and stdout SHA-256
`e1250273994c4ce3baa04a9b3308424f0fcd10e651a4d2115efef3e367160e8b`.
Both recorded executions exited `0`.

Their predicates support the exact algebra, cells, signs, tensor formulas,
and scale arithmetic described above. They do not prove a Banach IFT,
two-sided Hardy/common-range theorem, global symbol bound, nonnormal graph
resolvent, or Riesz projection. No formula or implementation changed, so the
review correctly reused these receipts rather than rerunning them.

## Correction check and final verdict

The joined correction fully repairs the load-bearing scientific
classifications: `0052` now separates HJ2 from the independent nonnormal
resolvent hypothesis GR, and `0054` preserves the exact cells, shape
identities, toroidal geometry, metric-Leray algebra, and limited no-sheet fact
while returning every unproved global bridge to blocked or conditional status.

The same-transaction clerical reconciliation removed the one stale
arXiv:1910.07493 inventory sentence from `0054/source-audit.md`. Its final hash
is `9d5c7de1cf150438a3a610e214636855f0ea2d590340257af747c00c4dd89519`;
the final receipt pins both the intermediate and final hashes. The requested
independent sentence/hash verification confirmed the removal and receipt hash.
The bounded stale-claim scan therefore passes with no known contradiction
across the claim-bearing target artifacts.

The final joined scientific verdict is therefore:

- **Unit A:** established for the corrected resonance/scaling, exact local
  cells, unequal-volume representation, toroidal distance, and physical KKS
  normalization; blocked for the coupled all-sector graph inverse by GR.
- **Unit B:** established for exact local cells, pointwise whole-space shape
  identities, toroidal geometry, algebraic metric-Leray differentiation, and
  the limited `p>=6` no-sheet result; blocked for HJ2.
- **Joined transfer:** established only as the conditional implication
  `HJ2 + GR => graph-resolvent/Riesz transfer`. No actual Cao Riesz family or
  nonlinear branch is established.

The next scientific construction is not another scaling oracle. It is an
exact augmented Cao branch on declared Banach spaces, a two-sided common DA
graph/Hodge theorem with complete interface traces, genuine two-index symbol
bounds, and the full nonnormal `X_* -> D_*` Grushin resolvent. The parent
campaign remains active. Nothing here licenses nonlinear stability, P2/P4,
particle, or quantum conclusions.
