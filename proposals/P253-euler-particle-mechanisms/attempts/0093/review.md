# P253/0093 independent review of corrected 0092

## Integrity and bounded correction

This review executed against the centrally activated contract at README
SHA-256
`da86576b5c1e80c25aaf2915a20a8b286af5571bdb62e7833527d43c70e26f20`.
`particle-balance-review` authored or implemented none of 0092. P253/0088
and P253/0094 remained unopened and excluded. No production numerics ran.

The single bounded correction passed. Its final receipt
`0092/0093-bounded-correction-receipt.md` has SHA-256
`15c9d1cab6b221f4bda166d6f7330905de8753d161b4e2e70e32dbf25d6d7fe8`.
Final corrected hashes are:

- README: `2e8138c4871f6548fc9f8a215fd747252a33bf2d4f6099afc1d7f04c7ae85dd0`;
- derivation: `3ce9062d93a43ea3d3f65eb817ff2fc3d1c9e98da47f83ded583ea4e184f7cbc`;
- source audit: `dd170288f4f7d3b82b214532d9bd808d5b556c2afc15dd0cd9afa5aee48f3daa`;
- result: `a28ea21fcbc42b383489089ff475176a54d63d98313e2ba513075436f7b5b678`;
- validation: `29797a58b37fa5f587769c57c9250f4cdb135f45f623f6dfd2fe759b95edd4e6`;
- verifier: `86ca3ae3859abb42ce2db0a6c572b96e62d498ebeefc6f97caa79a48257b2509`;
- public API: `c52dc54d5e51f85cf65ffca95ccac7eed518c6ab3c6717f3c5bb63c173854172`;
- public test: `807d345b9d311df88a98dd893b4d32f3e147029e7b8d5817ea6b4654a673089d`;
- author completion: `3d72173f592ade4096ebd230a1d6f76d83f683901f3c81ee690a753258c739c4`;
- manifest: `f5421864a03c018e2fa7e0136c2c47e2bdca39bad96c4f79bd7752dbb96c4ff0`.

Because the API and verifier wording changed, their focused and exact checks
were rerun by the author: eight tests and twelve exact rows passed with exit
zero. The final `Q/g_0` repair was wording-only and caused no second rerun.

## Source and evidence applicability

The imported campaign results are used only at their independently reviewed
scopes: bare-Euler similarity, Schwinger--Hopf action normalization, the
classical transported tag current, and the parameterized charged Cao branch.
None supplies a compact matter representation, a quantum momentum spectrum,
a compact two-form gauge field, a BF level, or a same-carrier action-selection
law.

The classifications of principal `U(1)` bundles by `H^2`, the cellular
cohomology of spheres, and `T^*S^1=S^1 x R` are standard exact mathematics and
are applied correctly. The API's topology and BF booleans are typed ledgers,
not executable proofs of those facts. The symbolic verifier independently
checks finite algebra and catches normalization drift; it cannot quantize a
classical phase space or establish a physical topological supplier.

## Unit A — dimensions, gauge rescaling, and real-tag continuum

With charge-explicit dimensions,

    S_EM(Q)=Q^2/(4 pi epsilon_EM c_EM)

has dimensions of action. Under `A'=aA`, `phi'=a phi`, the coefficients
transform as

    epsilon'=epsilon/a^2,
    mu'=mu a^2,
    g_tag'=g_tag/a.

Thus both `c_EM^2=1/(epsilon mu)` and `g_tag^2/epsilon` are invariant. A
positive `a` is only the normalization magnitude; the signed physical tag
coupling remains separate.

The passive material equation is linear in the real tag. Hence every real
`lambda` gives another admissible tag leaf with

    Q_lambda=g_tag lambda C_chi.

Fixing one leaf conserves its charge but neither selects that leaf nor creates
a smallest nonzero value. The charged Cao branch likewise takes its tag
integral as input.

**Unit A verdict: exact dimensions, rescaling invariants, and signed current
are established; charge quantization and nonzero selection by the existing
real-tag route are refuted by the continuous family.**

## Unit B — compact phase and continuous classical momentum

The corrected convention is dimensionally complete:

    S_phase=integral [S_0 n D_t theta
      -g_0 m n(phi-u dot A)-H] dx dt,
    kappa_m=g_0 m/S_0.

Here `n` is a dimensionless number density, `S_0` is an imported action,
`g_0` is an imported charge unit, and `m` is the chosen integer character.
The transformation `theta -> theta-kappa_m Lambda` cancels the Maxwell gauge
change exactly. The physical symplectic form and moment map are

    Omega_phase=S_0 integral delta n wedge delta theta,
    J_phase=S_0 N_phase,
    Q=g_0 m N_phase=kappa_m J_phase.

The circle coordinate and integer character are genuine, but the cotangent
fiber and `N_phase` remain real and continuous. The character is chosen
representation data; it does not select `S_0`, `g_0`, or a classical momentum
lattice. The correction also correctly keeps `g_tag` unrelated to
`kappa_m` absent a separate same-carrier map.

**Unit B verdict: the compact character and gauge cancellation are
established, while classical charge quantization from phase periodicity is
refuted.** Prequantization or a quantum charge spectrum needs an externally
normalized action unit and an actual quantization rule.

## Unit C — whole-space topology and electric boundary charge

For a smooth bundle extending over compactified whole space,

    H^2(S^3;Z)=0,  pi_3(U(1))=0.

The spatial bundle is topologically trivial. By contrast,
`H^2(S^2;Z)=Z` classifies a magnetic Chern flux on a linking or asymptotic
sphere. A nonzero value cannot extend smoothly through all of `R^3`; it needs
a removed/singular core or a nonextendible bundle.

Electric charge

    Q=integral_(S^2_infinity) epsilon_EM *F

is the Gauss-law boundary moment map. It is not that magnetic Chern class.
Relating the two needs additional duality, monopole, higher-form, and physical
normalization data.

**Unit C verdict: established as stated.** Smooth compactified whole-space
`U(1)` topology supplies no electric Chern quantization, while the magnetic
`S^2` route remains a distinct candidate with an altered core/domain.

## Unit D — Euler vorticity and the BF supplier boundary

For a smooth whole-space Euler velocity, the vorticity two-form is the real
exact form `d u^flat`. Material advection preserves its real continuum data;
it does not create large compact two-form gauge transformations, integral
periods, or an integer quantum BF level. In a compact `3+1` BF action, the
two-form connection, compact gauge equivalence, exponentiated-action
normalization, and level are extra structures.

**Unit D verdict: the direct Euler-vorticity BF route is refuted at its stated
supplier claim.** This does not refute a newly constructed compact defect or
higher-form sector; it identifies exactly what that continuation must add and
derive.

## Unit E — fixed-charge action rank and conditional root

On the product classical directions `(lambda,J)`, the map

    (lambda,J) -> (Q,J,S_EM(Q))

has rank two when `g_tag C_chi` is nonzero. Its second column is exactly
`(0,1,0)^T`, so the fixed-charge tangent `(0,delta J)` leaves `Q` and
`S_EM(Q)` unchanged while varying the carrier-mode action. The existing
charged branch contains no row that removes this independent direction.

If a different common-phase construction explicitly imposes

    Q=kappa_m J,  J=S_EM(Q),

then the nonzero algebraic root is

    J_*=4 pi epsilon_EM c_EM/kappa_m^2.

It is positive for either sign of `kappa_m` and has the correct action units.
But `kappa_m` is a charge-per-action coefficient, not `g_tag`, and the second
equation is precisely the missing selection law rather than an Euler--Maxwell
equation. Positive uncoupled classical energy instead has its minimum at
`J=0`.

**Unit E verdict: fixed-charge selection is refuted for the existing
classical product extension; the common-phase nonzero root is established
only as conditional algebra.** No existence, accessibility, attraction, or
physical selection theorem follows.

## Unit F — parent boundary

The exact route results neither adopt the `U(1)` extension as bare Euler nor
close the active compact-defect/BF or response routes. They supply no
nonlinear P2 persistence, Born probability, autonomous reset, exchange or
fermionic character, particle, P4/P5 completion, electron, or neutrino.

**Unit F verdict: established boundary.** These are route-scoped classical
no-selection statements, not a universal action/charge no-go or parent
campaign verdict.

## Strongest supported result and next construction

Corrected 0092 proves that `S_EM(Q)` is an exact normalization-invariant
action scale, while each tested existing classical mechanism fails to select
it as a nonzero carrier action: the real tag has a signed continuum of leaves,
a compact phase retains continuous cotangent momentum, smooth whole-space
topology has no electric Chern class, direct Euler vorticity supplies no
compact BF datum, and fixed charge leaves an independent action direction.

The next positive construction is a genuine same-carrier compact defect or
higher-form sector whose compact connection, integer level, charge unit, and
symplectic coefficient are derived rather than inserted, followed by a
nonzero constraint or energetic mechanism that removes the fixed-charge
action direction. The physical response/Feshbach, P2 persistence, and
Born/reset/exchange/neutrino obligations remain separate active routes.
