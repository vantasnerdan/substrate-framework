# P253/0098 independent final review of P253/0094

Reviewer: `particle-balance-review`  
Target owner: `root`  
Activated review README SHA-256: `2cab0310da57712bc684811dcd30aadff15b699ebc06ec232131ad0ddb2ff17e`  
Frozen target manifest SHA-256: `7d60b60e930e65ecad8d9ded97834095228ce57f897f0485f19e89f72415e0c3`

## Boundary and provenance

This was a content-blind, independent, non-author review of final P253/0094.
All frozen target hashes matched. P253/0095, P253/0097, and P253/0100 were not
opened or used.

The manifest intentionally mixes attempt-relative paths with repository-root
`src/` and `tests/` paths. A one-directory `sha256sum -c` invocation is
therefore not its convention. Resolving bare entries against `attempts/0094`
and `src/`/`tests/` entries against the repository root verifies every entry.
The recorded exact oracle has twelve PASS lines and exit zero; the focused
suite has eight passing tests. Neither was rerun because the analytic claims,
not a changed implementation, were under review.

## Source and dependency applicability

The reviewed suppliers are used within their adjudicated boundaries. In
particular, the positive-Krein Cao modes do not provide an absolute nonlinear
reduction; the prescribed-current Maxwell shell theorem does not provide its
principal-value real self-energy; the strong-saddle result prevents a
two-mode energy from becoming a full persistence theorem; and the compact
phase/BF atoms from P253/0092 are consumed only through the final P253/0093
review.

The new positive formulas are derived directly rather than attributed to an
external theorem: the scalar implicit-function argument, phase moment ledger,
Alexander-duality ranks, closed two-label form, Euler similarity, and the two
regularized circular double integrals are all exposed. No source is asked to
supply the missing same-carrier invariant manifold, compact two-form field,
or finite-core charged branch.

## A0 — local reduced-energy sign theorem

For

    E_red=E_0+nu*J+a_2*J^2+R(J),  J>=0,

with `a_2>0`, `R` of class `C^2`, `|R'(J)|<=C J^2`, and
`|R''(J)|<=C_1 J`, the scalar equation
`nu+2*a_2*J+R'(J)=0` has the unique local branch
`J(nu)=-nu/(2*a_2)+O(nu^2)`. For sufficiently small negative `nu` this is a
positive strict local minimum; for positive `nu`, zero is the sufficiently
local boundary minimum. The theorem is conditional on the displayed reduced
energy and its radius, exactly as stated.

**A0 verdict: established.**

## A1 — actual Cao reduction and coefficients

The target correctly identifies, but does not construct, the load-bearing
coefficient calculation. The quartic coefficient requires the direct fourth
variation together with complement slaving on a common constrained domain.
The joint strong saddle and Maxwell continuum make a naive full-space inverse
or two-mode truncation invalid. The real Maxwell shift separately requires a
renormalized principal-value reciprocal Feshbach construction; outgoing shell
power supplies only the dissipative boundary value. No existing positive
oscillator or Coulomb term earns the negative linear coefficient needed by
A0.

**A1 verdict: blocked.** The missing construction is an actual same-carrier
invariant or constrained reduction, its third/fourth variations and complement
solution, and the real principal-value Maxwell term with analytic remainder
bounds.

## B — compact phase character

The typed first-order phase action gives

    J_phase=S_0*N,
    Q=g_0*m*N=(g_0*m/S_0)J_phase.

The integer `m` is chosen representation data, whereas the classical `N` and
the action coefficient `S_0` remain continuous. The family
`S_0=C_0*g_0^2/(4*pi*epsilon_EM*c_EM)` is invariant under the stated field
rescaling for every dimensionless `C_0`; dimensional analysis therefore does
not select it. At fixed charge, changing `S_0` changes `J_phase` continuously.

**B verdict: refuted as a classical action selector.** The mechanism is the
free continuous action coefficient and classical number, not the absence of a
valid compact character.

## C1 — ordinary Euler vorticity as compact BF connection

On smooth whole space the Euler vorticity two-form is a real exact advected
form. Its circulation/flux data are continuous and it has neither a compact
two-form gauge identification nor an integral BF level. Treating the real
form as a gerbe connection would add structure not supplied by Euler.

**C1 verdict: refuted.**

## C2 — smooth two-circle-label pullback

Alexander duality gives `H^2(R^3-{point};Z)=Z`, whereas a ring complement has
`H^1=Z` and `H^2=0`. The ring row is linking circulation, not a point-charge
Gauss surface. For smooth circle labels,
`B=d alpha wedge d beta` satisfies `dB=0` on every defect-free chart, so
variation of `B wedge dA` supplies no smooth bulk electric BF source.

**C2 verdict: refuted as a smooth electric BF supplier.** Distributional
label defects are outside this route and activate C3 rather than invalidating
the calculation.

## C3 — genuine compact carrier defect

A singular compact two-form connection could evade the smooth `dB=0`
mechanism. The target does not construct its compact periods, distributional
source, finite core energy, gauge boundary terms, BF level, or same-carrier
action coefficient.

**C3 verdict: blocked at the named compact-defect construction.**

## D — coupled defect-energy realization

An integer character may label a branch but does not by itself add a negative
detuning or binding term. With only positive oscillator and Coulomb blocks,
the internal action remains minimized at zero. A chemical potential, binding,
anomaly-inflow term, or discrete admissibility row would have to be derived
from the same carrier and inserted consistently into A0.

**D verdict: blocked.** The strongest result is the exact identification of
the additional signed energy row that a successful continuation must supply.

## E0 — fixed-aspect leading ring balance

The circular filament double integral gives the leading slender-ring fluid
energy

    E_fluid=(rho_m*kappa^2/2) R L,

and the regularized circular line-charge Newton integral gives

    E_C=Q^2 L/(8*pi^2*epsilon_EM*R).

At fixed aspect/profile the common leading logarithm `L` is constant under
the model radius variation. Therefore `L*(A R+B/R)` has the unique positive
critical radius

    R_*=|Q|/(2*pi*|kappa|*sqrt(rho_m*epsilon_EM)),

positive second derivative, and thin-ring axial impulse

    I_z,*=Q^2/(4*pi*epsilon_EM*kappa).

The signs, density, `pi` factors, and signed circulation dependence check.
This is a leading fluid-plus-electrostatic, fixed-aspect asymptotic. It selects
radius and translational impulse only. It is not a minimum of the full charged
Cao energy, an internal KKS action, an action quantum, or P2.

**E0 verdict: established at exactly that leading-asymptotic scope.**

## E1 — exact finite-core charged Cao scale selection

The full theorem would need a fixed-label admissible radius path or a direct
finite-core expansion including the actual tag profile, comoving
magnetic/current energy, finite logarithmic parts, varying core/aspect rows,
and the subluminal admissibility window. It would then need a separate map
from selected translational impulse to an internal mode action. None is
constructed here, and E0 is correctly presented as a comparison target.

**E1 verdict: blocked at the fixed-label finite-core energy and constraint
calculation.**

## Evidence and combined verdict

The importable helpers and tests accurately encode the finite algebraic
ledgers. The API deliberately restricts the leading-ring helper to the
positive branch; the analytic derivation carries the absolute values and the
signed impulse. The scripts do not pretend to construct A1, C3, D, or E1.

P253/0094 is **established at its route-scoped author boundary**. A0 and E0
are genuine positive theorems in their stated conditional/asymptotic classes;
B, C1, and C2 are exact mechanism-specific refutations; A1, C3, D, and E1 are
blocked at explicit constructive dependencies. No bounded correction is
required. The strongest next achievement is the exact fixed-label charged
Cao energy derivative, in parallel with either a same-carrier nonlinear
reduction supplying the A0 coefficients or a genuine compact defect with a
derived BF/action coefficient.
