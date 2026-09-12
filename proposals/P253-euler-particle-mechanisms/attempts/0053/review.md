# P253/0053 independent review of the 0048 finite-window carrier block

## Frozen boundary and correction provenance

This is the preregistered non-author/non-implementer review of root-owned
attempt `0048`. The frozen review README has SHA-256
`eb216dc45df9b98f4e9935bdaf96b780c43506281458d6d9efc1cb2ad5777de3`.
Central activation exited exactly `0`; the command, stdout, empty stderr, and
exit hashes are respectively
`5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f`,
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`,
and `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.
The accepted authority base is release `v0.183.0`.

The substantive review requested one bounded claim-domain correction. Its
receipt is `0053/target-correction-receipt.md`, SHA-256
`8c4437042bfe37270673bc295653b931ab700bf6be988f185eb4d3c4b2a49d55`.
The single bounded check verifies the complete correction map:

| Artifact | Frozen SHA-256 | Corrected SHA-256 |
|---|---|---|
| `0048/derivation.md` | `3b35b7cb56089a4fb4e78fff535f5f72128fc39c9d2cba9d1299f2d0a32560fd` | `381c59e5debe1886847acd296c795192605576cdf55afc721688d2aea9783b25` |
| `0048/result.yaml` | `13569f910e968012815e68be56cdeb8eff6f5286e03a525e6055ce9692687b29` | `e6b09d0f5e01d80c03aeba2987f0406695ca83f7abf546524d16e837d1cf894b` |
| `0048/source-audit.md` | `d465aab71ef242dd62597086f8401816b4a7ac40eda6e4799d73cb09d23e2577` | `4db69ddded6fd2d84e749271e2f13efb293c00645a8c8efe311cb35820971b03` |
| `0048/validation.md` | `3107ddb505a0cf16273d182176e186115b4e9c9ed21ce727d192bf46c3db2202` | `9ac3c21c514e4024a48f57b0fcc30599f64316aec037008357a87b15ab6490b5` |

The target README, verifier, supplier API, and supplier tests remain at their
frozen hashes. The verifier was deliberately not rerun: its existing fifteen
predicates remain algebraic regression evidence and do not prove the missing
operator constructions. Active `0052` was not inspected or used.

## Source applicability

The three primary sources are correctly pinned and narrowly used:

- Cao--Lai--Qin--Zhan--Zou, arXiv:2206.10165v2, cached PDF SHA-256
  `6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`,
  supplies the exact rescaled meridional equation, Lane--Emden core,
  limiting nondegeneracy, local free-boundary graph facts, and the stated
  `O(epsilon^2 |log epsilon|)` profile remainder. It does not supply a
  nonaxisymmetric Euler graph-Riesz family or the full moving-boundary
  Green/Leray jet.
- Gallay--Smets, arXiv:1805.05064v3, cached PDF SHA-256
  `081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`,
  supplies the column Fourier equations, axisymmetric energy identity, and
  cylindrical Biot--Savart estimates. Its enstrophy-space resolvent is not a
  compact dynamically-accessible graph-contour resolvent for the 0048 ring.
- Fischer--Schopohl, arXiv:cond-mat/0008215v3, cached PDF SHA-256
  `0b69b8cae5e758d860dcf8188260bfd79839555d173797d0a74ee487f30e5977`,
  is used only as a filament-spectrum comparator and supplies neither a
  finite-core coefficient nor a Riesz projection.

The corrected source audit now makes these transfer boundaries explicit.

## Check 1: resonance index and finite spectral-window scaling

For fixed azimuthal `l>=2`, the column asymptotics are

    sigma_n = l delta L_Phi/(pi n)+O(delta/n^2),
    tau_l = c_l delta^2 L_delta+O(delta^2),
    c_l=l sqrt(l^2-1),
    L_delta=log(1/delta).

Solving `sigma_n=tau_l` therefore gives

    n_*=l L_Phi/(pi c_l delta L_delta)
         +O(1/(delta L_delta^2)+1).

The first error is the relative `O(1/L_delta)` correction in the bending
frequency and the `O(1)` term includes the Weyl phase and integer rounding.
The pre-correction expression omitted the factor `l/c_l` before substituting
`c_l` and understated the absolute remainder. The corrected equation (33)
is supported.

The energy-normalized first curvature map has
`g_n=O(delta^2/n)`. At `n_*`, this is `O(delta^3 L_delta)`, while adjacent
Kelvin levels are separated by `Theta(delta^3 L_delta^2)`. Hence the
coupling-to-spacing ratio is `O(1/L_delta)`, and the finite spectral coupling
window contains at most one positive axisymmetric Kelvin mode together with
its Hamiltonian conjugate. The resulting one-crossing interval has width
`O(delta^2)`; summing over the `O(1/(delta L_delta))` crossings in a dyadic
shell gives bad-set measure `O(delta/L_delta)` and relative measure
`O(1/L_delta)`. These are asymptotic spectral-window statements, not temporal
Duhamel propagation.

**Verdict:** established at the corrected fixed-`l`, thin-ring asymptotic
scope. It does not construct a graph projection or invariant Euler block.

## Check 2: fixed-domain Hanzawa/Piola chart and graph equivalence

No orientation-preserving Hanzawa map, Piola pullback, transformed pressure
and Leray operator, axis/interface/exterior graph domain, or two-sided graph
norm equivalence is constructed in `0048`. The corrected derivation collects
exactly those requirements in the explicit hypothesis `HJ2` and makes the
operator expansion (31) conditional on it. Local Cao profile coordinates and
coefficient regularity do not fill this gap.

**Verdict:** blocked by the named `HJ2` construction. The correction removes
the prior unsupported implication without weakening the exact local cells or
spectral scaling.

## Check 3: full complement resolvent

The sequence algebra proves the axisymmetric estimates

    sum_(n outside I_delta) |g_n|^2/|sigma_n-tau_l|
      =O(delta^3 L_delta),
    sum_(n outside I_delta) |g_n/(sigma_n-tau_l)|^2
      =O(delta/L_delta).

Thus the `m=0` Schur tail is one factor `delta` below the bending scale and
the first homological transform is `O(sqrt(delta/L_delta))` in the declared
energy sequence norm. This is a useful tame sequence result.

It does not control every azimuthal sector, the axis, collar and transmission
conditions, the exterior-pressure tail, neutral modes, adjoint graph domain,
or a uniform contour inverse. The corrected text now defines `P_delta` only
formally and states that (44) is a conditional analytic-Fredholm
reconstruction; no Riesz projection is claimed as established.

**Verdict:** the `m=0` sequence Schur estimate is established; the full Euler
complement resolvent is blocked by an all-sector graph-domain construction.

## Check 4: KKS, energy, and generator normalization

The axisymmetric eigen-equation and normalization give

    ||v_n||_2^2=M_n/sigma_n^2,    1<=M_n<=2,

and therefore positive constrained-energy signature for every nonzero
axisymmetric pair. With the corrected generator convention
`A e_n^+=-i sigma_n e_n^+` and an energy-unit complex vector,
`Omega(Ae,conjugate(e))=H(e,conjugate(e))=1` gives

    Omega(e_n^+,conjugate(e_n^+))=i/sigma_n.

The sign is therefore consistent. The magnitude scaling
`sqrt(|sigma_n|) e_n^E` and same-Krein character of the conditional finite
matrix are supported. However, `0048` does not carry the density,
azimuthal/axial Fourier volume, Jacobian, and real-to-complex conventions
through both the bending and Kelvin coordinates. The corrected artifacts now
leave that exact physical normalization open.

**Verdict:** positive energy/Krein sign and frequency-magnitude scaling are
established; exact physical KKS normalization and the associated physical
Feshbach matrix remain conditional.

## Check 5: exact Cao cells versus the global operator jet

Expanding the exact source equation

    -Delta w+[q/(1+qx)] partial_x w=(1+qx)^2(w_+)^p

with `w=U+qV+q^2Z+o(q^2)` gives

    L_U V=2xU^p-partial_x U,

    L_U Z={p(p-1)/2}U^(p-2)V^2+2pxU^(p-1)V+x^2U^p
          -partial_x V+x partial_x U.

The parity assignments—`V` odd and `Z` even—and the centered limiting
nondegeneracy are consistent. A local centered fixed-circulation
`q^2 log(q)` cell solves the homogeneous even equation and hence vanishes.
These are exact local profile-cell conclusions. They do not improve Cao's
stated remainder into a twice differentiated whole-space Euler generator.

The corrected derivation now separates these cells from the conditional
`HJ2` expansion containing the Hanzawa/Piola, full Green/Leray, pressure,
interface, collar, and exterior contributions.

**Verdict:** the first and second Cao cell recurrences and parity are
established; the moving-boundary global graph jet is blocked by `HJ2`.

## Check 6: fixed-order residual and critical-harmonic boundary

The product estimate, transport divisor

    D_(m,j)(r)=m Omega(r)-j l Omega_pat,

and fixed-order condition
`N l |Omega_pat|<Omega(a)/2` support a formal induction for each preassigned
`N<=p-2` after assuming the missing chart, complement inverse, and physical
Riesz projection. Equations (52)--(53) therefore specify a conditional formal
residual of order `O(|a|^(N+1))`; they do not construct coefficients on an
actual Euler branch. At fixed nonzero `delta`, the first internal critical
harmonic has scale

    j_crit=Theta(1/(delta^2 L_delta)).

This explains why every fixed order is noncritical as `delta` tends to zero
but does not establish convergence or transparency at all orders.

**Verdict:** the divisor scale and conditional formal fixed-order residual are
established as algebra; an approximate Euler evolution, nonlinear branch,
and all-time result are not established.

## Verifier and supplier-code scope

`0048/verify_finite_window.py` remains at SHA-256
`3762c1d3c238f3e7d5e02fe7a2d2700547cdb0b41da5d9d95863e6a1514221cf`.
Its final repository receipt records fifteen exact symbolic predicates with
stdout SHA-256
`0edb7b7a3e962cb34a3a72dbec5e96b9bb4531555bd9aab80701125bd28f24f4`
and exit `0`. Those checks expose sign, recurrence, scaling, and discriminant
errors; they are regression evidence for the finite algebra only.

The unchanged supplier API and tests have SHA-256 values
`f2a2c58235a0f2c8ee358123b383b8e919d22ad095bc0f70b67142c7fdbd5214`
and `f7ae11d6a2bfc292675740e43357d9140257f5a206a582fafb903d03042956cc`.
They do not materialize the missing Hanzawa/Piola or graph-resolvent theorem.
No review rerun was needed because the correction changed no equation, API,
test, or verifier predicate.

## Final verdict

The corrected `0048` transaction is **established at scope
`EXACT_CAO_CELLS_AND_FINITE_SPECTRAL_WINDOW_SEQUENCE_SCALING`**. Its strongest
result is the exact first/second Cao profile-cell calculus together with the
corrected resonance location, `O(delta^2/n)` coupling, `O(1/log(1/delta))`
coupling-to-spacing and bad-measure ratios, and the
`O(delta^3 log(1/delta))` axisymmetric Schur tail. The same-Krein
avoided-crossing calculation is a correct conditional finite-matrix theorem.

The stronger fixed-domain Euler graph-Riesz and branch objective remains
**blocked by named constructions**, not refuted: construct `HJ2`, prove the
all-sector/collar/interface/exterior-pressure complement and exact physical
KKS normalization, then address critical-layer transparency or an enlarged
streamline-vorticity range before convergence can be claimed. Nothing here
licenses a temporal finite-window theorem, nonlinear stability, P2/P4,
particle, or quantum conclusion.

The single correction check found no stale claim in the corrected
claim-bearing derivation, result, source audit, or validation that promotes a
graph-domain Riesz family, exact physical KKS normalization, global graph jet,
or constructed Euler branch. The unchanged preregistration README remains an
aspirational target, not a result artifact. YAML parsing, terminal-newline
checks, and the scoped final `git diff --check` pass. No further correction is
needed.
