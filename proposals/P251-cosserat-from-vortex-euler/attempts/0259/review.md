# Independent review of the exact 0250 isovortical channel

Reviewer: `herdr geometry-review pane w3:p4`, fresh non-author and
non-implementer of the 0250 response construction.  Review boundary:
commit `47b1d25b5e01dc58762310d65499ff824940a5d3`, specifically the frozen
`0250/README.md`, `0250/receipt.md`, and `0250/isovortical-channel.md`.
Only the last artifact's claimed preparation, pressure, observation, and KKS
rows are adjudicated here.  Later `0250` working-tree additions, joint gain,
normalization, geometry/density, and the parent objective are outside this
review.

`route_verdict: established as stated at exact isovortical preparation and
finite-window asymptotic-response scope`.

`evidence_scope: independent analytic review; no new numerical oracle and no
claim promotion`.

## Claim and positive role

The reviewed result fixes one sufficiently large but finite 0211 stationary
Euler ring and a compact regular action annulus on which
`curl u=lambda u`, `lambda ne 0`, and the flux rotation has strict twist.  It
claims:

1. a smooth, compactly action-supported, divergence-free Lin generator for
   every selected nonzero toroidal harmonic, with exact initial Euler velocity
   `P(xi cross omega)` and polynomial carrier cost;
2. a nonempty open transport-frequency subband and an arbitrary-order
   full-pressure finite-time Euler/Lin expansion on that same fixed ring;
3. exact initial tagged centroid, covariance-angle, angular-displacement `G`,
   and mechanical-spin/current `S` formulas, including specifically exhibited
   transverse nonzero coefficients but no joint-rank claim; and
4. an exact zero KKS pairing between the real and imaginary quadratures of one
   transport source.

Its positive role is to replace the singular `k dot omega=0` principal Kelvin
inversion by an actual finite-`N` isovortical initial preparation.  It does not
claim an exact Euler spectral band, an invariant nonlinear manifold, joint
acoustic/optical range, a nondegenerate action normalizer, the longitudinal
covariance component, transfer to 0145/0147, or any geometry/density result.

## Frozen transaction and provenance

The 0250 schema receipt at the reviewed commit reports 270 accepted claims and
exit zero.  Independent hashes reproduce the two content pins relevant to the
claim:

| Frozen artifact | SHA256 |
| --- | --- |
| `0250/README.md` | `28f206cf88723733422e85c4295a8f220a2075671c7ad80c746fb67e45a83c6c` |
| `0250/isovortical-channel.md` | `492ecd9b5d4e27bed71757e15a7ff2c049e39134c408bbd337d205ea5f3c54ff` |
| `0250/receipt.md` | `01df693448d3afe4c39473b24e0f24490558bfc01fe0708da6d75fdc5269a34e` |

The isovortical hash agrees with the author receipt.  The reviewed delta is
proposal evidence only; no accepted claim, source module, registry, release,
or generated consumer changes at this boundary.

## Strongest supported positive statement

On the exact fixed 0211 constant-curl annulus, choose a compactly
action-supported complex amplitude

\[
 f=A(I)e^{iNI}g(I,\theta)
\]

with one fixed nonzero toroidal harmonic and angular mean zero.  There is a
smooth compactly supported divergence-free generator

\[
 \xi=-\frac{\nabla I}{\lambda|\nabla I|^2}f+\beta u,
 \qquad
 u\mathbin\cdot\nabla\beta=-\operatorname{div}
 \left(-\frac{\nabla I}{\lambda|\nabla I|^2}f\right),
\]

on a nonempty open action subannulus, and it satisfies exactly

\[
 \xi\times\omega=-\frac{\nabla I\times u}{|\nabla I|^2}f.
\]

The associated initial velocity is the actual whole-space Leray field
`v_N(0)=P(xi cross omega)`.  For every fixed finite time window and finite
derivative inventory, taking the pressure hierarchy to sufficiently high
order and then choosing a finite carrier gives the stated asymptotic Euler/Lin
response, with the generator's `C^s` cost bounded by a fixed polynomial
`N^(s+1)`.  Varying the localized action center gives a nonconstant open
transport-frequency interval without changing the ring or tag.

The material formulas (16)--(25) and the zero in (26) are exact at their
stated scopes.  The strongest honest reading of “response band” is therefore
an exact family of finite-`N` Euler/Lin preparations whose finite-window
responses have an arbitrary-order transport-frequency asymptotic.  It is not
an assertion that the linearized Euler generator has an exact absolutely
continuous spectral band; the artifact already excludes that stronger
reading.

## Evidence map

| Evidence | Proposition established | Role | Bridge to reviewed claim | Limit |
| --- | --- | --- | --- | --- |
| 0211 construction and 0215 review | One fixed global steady ring has an open literal constant-curl inner annulus, regular flux action, and strict rotation twist | Reviewed exact dependency | Supplies `omega=lambda u`, action-angle coordinates, compact regular annulus, and varying frequency ratio | Retains its translating-frame far velocity; supplies no Lin response or observation gain |
| 0250 equations (3)--(13) | Exact cross-product sign, divergence correction, open-annulus cohomological inverse, and actual Leray initial datum | Exact proof | Replaces the singular principal Kelvin lift on the same fixed field | Does not compute a physical joint gain or a canonical pair |
| 0112 recurrence (5)--(8), independently rechecked here | Arbitrary finite pressure hierarchy and exact-solenoidal completion for a transported noncritical phase | Exact proof imported at its recurrence scope | The recurrence divides by `|grad I|^2`; its later steps do not require `grad I dot omega ne 0`.  The new generator supplies the initial lift formerly obtained from that condition | Supplies an asymptotic finite-time family, not spectral stability or an invariant finite-dimensional subspace |
| Friedlander--Strauss--Vishik, Theorem 4.1 | The standard Euler bicharacteristic amplitude system used as a convention check | Primary-source corroboration | Agrees with the leading full-pressure transport equation rederived in 0112/0250 | Its instability and essential-spectrum theorems are not imported |
| 0250 equations (16)--(25), with 0226 material definitions | Exact scalar moment, centered-reference shape, `G`, `S`, and `G_t-S` rows | Exact proof | Converts the actual divergence-free generator into literal material observations | Nonzero individual transverse entries do not imply common matrix rank or the absent longitudinal angle |
| 0250 equation (26), with the repository KKS convention | The two quadratures of one source have zero mutual KKS pairing | Exact proof | Both generators lie in `span{grad I,u}`, so their cross product is orthogonal to `omega=lambda u` pointwise | Does not decide cross-pairings with separate normalizer columns or the Jacobi energy |
| 0254 fixed-band bridge | Frequencies in an open interval densely synthesize finite-window `C^r` histories once normalized physical gain/residual hypotheses hold | Exact conditional theorem | Supplies the later finite synthesis and retained coefficient cost | Supplies no Euler gain, detector, normalizer, or geometry |

The primary source checked for the bicharacteristic convention was
[Friedlander--Strauss--Vishik, *Ann. IHP Analyse non linéaire* 14
(1997)](https://content.ems.press/assets/public/full-texts/serials/aihpc/14/2/4398291/online/10.1016-s0294-1449-97-80144-8.pdf).
The EPS source underlying 0211's already reviewed persistence license was
[Enciso--Peralta-Salas, arXiv:1210.6271](https://arxiv.org/pdf/1210.6271);
no new EPS conclusion is inferred in this review.

## Oracle audit

The strongest practical oracle is independent exact rederivation.  No
floating-point or long-script remainder is present.

### Preparation and cohomology

With `k=grad I`, `a=-k cross u/|k|^2`, and
`xi_perp=-k f/(lambda|k|^2)`, direct use of
`omega=lambda u` gives `xi_perp cross omega=a f`; adding `beta u` changes
neither side.  Since `div u=0`, `div xi=0` is exactly the cohomological
equation displayed above.  In volume action-angle coordinates,

\[
 \int_{\mathbb T^2}J(-\operatorname{div}\xi_\perp)d\theta
 =\lambda^{-1}\partial_I\left(JA e^{iNI}
   \int_{\mathbb T^2}g\,d\theta\right)=0.
\]

For a fixed toroidal number `n ne 0`, the denominators are
`Omega_theta(m+n rho(I))`.  On a compact regular annulus,
`Omega_theta` stays away from zero and `rho` is bounded; only finitely many
integers `m` can resonate.  Strict twist makes each such zero isolated.
Deleting those finitely many points and then taking a compact subinterval
gives one uniform inverse for all `m`, with large `|m|` denominators growing
linearly.  This establishes an open smooth subannulus, not merely a Cantor
set.  The action derivative of `e^{iNI}` accounts for the stated one-power
source loss.

### Full pressure and quantifiers

The leading divergence of `a f` has no order-`N` term because `k dot a=0`.
Whole-space nonstationary phase therefore starts the Leray expansion at
relative order `N^-1`.  The 0112 recurrence thereafter uses only the
nonvanishing of `|k|`; `k dot omega ne 0` entered its old construction of a
Kelvin displacement, not its pressure recursion.  The exact 0250 generator
provides the missing initial isovortical lift, and the neutral-polarization
identity supplies the leading transported amplitude.  Differentiated energy
estimates give any fixed Sobolev inventory after taking correspondingly more
recurrence orders; Sobolev embedding supplies a requested finite observation
derivative inventory.

The order is the frozen one: fix the 0211 ring and open subannulus; for a
requested `T,r,epsilon`, first choose finite synthesis/localization data and a
sufficient recurrence order, then choose finite carriers.  The background
radius never follows `T`, `r`, or `epsilon`.  The potentially large synthesis,
tag, localization, and `N^(s+1)` preparation costs remain finite and are not
renormalized away.

### Material and action rows

For a stationary tag `chi(I)`, compact support and `div xi=0` give

\[
 \rho\int\chi\,\xi\mathbin\cdot\nabla F
 =-\rho\int F\chi'(I)dI(\xi)
 =\frac\rho\lambda\int J\chi' fF.
\]

This proves the centroid and second-moment rows and shows exactly why the
large tangential correction cancels from scalar configuration moments.  The
displayed `G` and `S` rows agree with the repository material definitions;
substituting the Lin reconstruction gives
`G_t-S=-2rho integral chi xi cross u`.  For a pure harmonic,
`beta=Nf/(lambda nu)+O(1)`, yielding the claimed leading `G/S` relation.
The artifact correctly stops before promoting a nonzero individual
coefficient to joint range.

Finally, real and imaginary parts of the same complex generator remain
pointwise combinations of `k` and `u`.  Their cross product is parallel to
`k cross u`, whose dot product with `omega=lambda u` vanishes.  The mutual KKS
entry is therefore identically zero with the exact normalization, not merely
small at large carrier.

## Findings

| Finding | Direct evidence | Boundary classification and minimum correction | What would resolve or overturn it |
| --- | --- | --- | --- |
| No load-bearing contradiction or missing step occurs inside the stated isovortical preparation claim | Exact rederivations above; matching frozen hashes and schema receipt | Established as stated; no correction requested | A counterexample to the open-annulus inverse, the pressure recurrence at `k dot omega=0`, or one of the exact material identities would reopen this decision |
| “Response band” could be overread as an exact spectral band, but the proof explicitly defines a finite-window asymptotic preparation family and disclaims spectral stability | 0250 sections 2 and 4; 0112 recurrence scope | No text correction is necessary.  Downstream consumers must use only the finite-window preparation meaning | An actual spectral/adjoint theorem would support the stronger spectral reading |
| Joint acoustic/optical gain, the longitudinal covariance row, and the full KKS/Jacobi/current normalizer are absent | 0250 equations (18), (21)--(26) and its explicit remaining-construction paragraph | These are declared frontier, not defects in the reviewed statement; do not narrow the positive preparation result | Evaluate the common frequency-dependent gain matrix and all cross-form columns, then apply 0254 at their normalized gains |
| Transfer from the exact axisymmetric 0211 annulus to 0145/0147 is not supplied | 0250 section 2 explicitly excludes automatic transfer | Outside this reviewed claim and not a correction | A finite-order invariant-action/cohomology normal form with the carrier, pressure, and gain diagonal |

No prior independent review has narrowed this exact 0250 isovortical claim, so
there is no burden-shift conflict or recursive review shrinkage.

## Compatibility and consumers

The construction preserves the exact 0211 translating-frame field, its
constant nonzero curl factor, flux-action convention, and fixed finite radius.
It does not cancel the background's far velocity or relabel compact vorticity
as compact velocity.  The Lin preparations are finite-energy even though the
whole-space pressure projection has its physical exterior tail.  The positive
stationary tag is one fixed smooth function; its sparse carrier coefficients
and all preparation norms remain in the later diagonal.

The direct consumer is the 0254 conditional finite-band synthesis after a
physical nonzero gain is supplied.  Because this review changes no claim or
implementation, no source or downstream replay is triggered.  C-CST-017,
geometry/density, the acoustic hybrid map, and the joint normalizer remain
unchanged inputs or active frontier rather than review debt.

## Four-axis decision

- **Verification:** exact analytic identities and source-scope reconciliation
  pass; frozen schema and hashes agree.  No numerical verifier was needed.
- **Review:** accepted at the stated attempt-evidence boundary.
- **Compatibility:** compatible with the exact 0211 ring/action conventions
  and 0254's conditional synthesis theorem; no unclaimed transfer inferred.
- **Epistemic:** established proposal evidence, not an accepted-registry claim
  or parent theorem.
- **Relationship:** additive response construction on the fixed 0211
  background.
- **Strongest proposed statement:** the statement in “Strongest supported
  positive statement” above, with “band” understood in its explicit
  finite-window asymptotic-response sense.

## Promotion transaction and correction check

No promotion transaction is authorized.  The reviewed artifact remains under
the active P251 proposal, and the accepted registry, source code, release,
generated documentation, and central manifest are unchanged.

`correction_check: not needed`; no scientific correction was requested.

## Result and frontier

The exact divergence-free isovortical initial lift, open transport-frequency
subannulus, arbitrary-order full-pressure finite-time response, literal
material rows, polynomial source cost, and null mutual KKS entry are supported
as stated.  The decisive remaining 0250 construction is not more preparation
theory: it is the common frequency-dependent physical gain calculation,
including the ambient point/pressure acoustic row, the longitudinal optical
observation, sufficient independent `G/S` columns, and separate normalizer
columns with their complete KKS/Jacobi/current cross forms.  Geometry,
normalization, density, and parent completion remain outside this review.
