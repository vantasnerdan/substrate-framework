# P253/0069 joined independent review of 0064 and 0065

## Frozen transaction

This is one independent, non-author and non-implementer review of the
root-owned attempts `0064` and `0065`. The preregistered README has SHA-256
`7879bcdc902abc88304b4fd5397e530201e67f43b10529130a2bf0f4a38a282c`.
Central activation exited exactly `0`; the command, stdout, empty stderr, and
exit hashes are respectively
`0b8bcc78ad3e326535d62e436c7e5f623fef1e2c3c8d7c12ef1a792f07e8e1c8`,
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`,
and `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.
All preactivation target hashes matched the README freeze.

The review uses accepted release `v0.183.0`. It opened the frozen source
audits and the supervisor strategy/scope records before the claim bodies.
Attempts `0067` and `0068` and their gauge APIs/tests were not opened or used.
No production numerics or unchanged oracle was rerun.

One bounded target correction was requested and has passed its correction-only
check. It is limited to the scope of the broad 0064 headline and
zero-coefficient wording; every displayed compact-vorticity, pressure, flux,
transverse-projector, puncture, circulation, and helicity equation is retained.
The correction receipt has SHA-256
`68be057ffcc5fac71f5bc86f6459b3a46bf6f6b1bf49594802987a77e163fc88`.

## Source applicability and evidence roles

The compact-vorticity and puncture calculations are direct consequences of
the three-dimensional whole-space Biot--Savart formula, the Euler pressure
Poisson equation, de Rham cohomology of point and knot complements, the
divergence theorem, and the Fourier representation of isotropic multipliers.
No external electromagnetic equation or charge assignment is imported.

Accepted `C-CST-018` supplies a fixed prepared, finite-window linear
Euler/Lin continuum whose acoustic branch satisfies
`A_tt+a|K|^2 A=r_A` and whose physical coefficient is `mu=rho a`. It thereby
licenses the static transverse symbol `P_T/(mu |K|^2)` at the accepted
second-variation/preparation scope. It does not supply an autonomous defect,
an unrestricted invariant manifold, a scalar source, a constitutive carrier
force, recoil, a Coulomb law, or all-time propagation. The 0064 source audit
uses it only at that boundary.

The earlier impulse audit supplies corroborating convention provenance, but
the relevant moment contractions are rederived here. The 0060 compact-pressure
observation is likewise not load-bearing: exterior pressure constancy and the
virial density factor are derived directly. The 0064 supervisor strategy
record correctly removes autonomous-recoil language from the prescribed-force
quadratic comparator. The 0065 pre-review record correctly restricts the
multiplier theorem to linear translation-invariant maps at nonzero wavevector
and repairs the polar/axial parity statement.

The `euler_charge_multipole.py` API evaluates exact leading coefficients and
domain checks. Its six focused tests and the nine/five assertion verifiers are
algebraic regression evidence. They do not prove integration-by-parts
hypotheses, asymptotic remainders, de Rham classifications, local Euler
existence, autonomous recoil, or a particle interpretation. Those bridges are
supplied only where the analytic derivations below close.

## Unit A: 0064 compact multipoles and transverse-source boundary

### Compact-vorticity identities

For smooth compactly supported divergence-free vorticity,

    integral omega_i = 0,
    M_ij := integral x_i omega_j = -M_ji,
    I_a := (1/2) epsilon_aij M_ij,
    M_ij = epsilon_aij I_a.

The signs follow by integrating `partial_j(x_i omega_j)` and
`partial_k(x_i x_j omega_k)`. Expanding the Biot--Savart potential and taking
its curl gives the first permitted velocity coefficient

    u(x) = {3 n (I dot n)-I}/(4 pi r^3) + O(r^-4).

The coefficient is leading when `I` is nonzero; if it vanishes, decay begins
at a higher multipole. For two fixed compact carriers at separation `d n`,
the total kinetic energy's cross term is

    E_12 = rho I_1i I_2j partial_i partial_j(1/(4 pi d)) + O(d^-4)
         = rho {3(I_1 dot n)(I_2 dot n)-I_1 dot I_2}/(4 pi d^3)
           + O(d^-4).

The sign agrees with direct contraction of the two antisymmetric first
moments. Its translation derivative is `O(d^-4)`. The displayed `d^-3`
coefficient is the leading term only when its orientation contraction is
nonzero. This proves an anisotropic compact-vorticity interaction, not a
mechanical point force or scalar charge.

### Pressure and flux

With physical pressure and
`p=rho (-Delta)^(-1) partial_i partial_j(u_i u_j)`, two integrations by parts
remove the pressure-source monopole and dipole under the stated weighted
moment and boundary hypotheses. Moving both derivatives onto
`G=1/(4 pi r)` yields

    p(x) = rho {3 n_i n_j-delta_ij}
                 integral u_i u_j /(4 pi r^3) + o(r^-3).

This is the first permitted pressure coefficient, not an unconditional
nonzero leading term: it vanishes for an isotropic kinetic tensor, as the API
test itself checks. For a stationary compact-velocity field, `u=0` in the
exterior and Euler gives `grad p=0`; after subtracting the exterior constant,
the exact virial row is

    integral u_i u_k = -(delta_ik/rho) integral p.

All signs and density factors are correct.

The radial field `u_q=q x/r^3` has flux `4 pi q`, distributional divergence
`4 pi q delta_0`, and divergent core energy. Smooth source-free whole-space
incompressibility therefore excludes this isotropic flux monopole. This does
not exclude tangential, zero-flux `r^-2` tails.

### Local isotropic scalar-to-transverse obstruction and collective action

At each nonzero wavevector, a linear translation- and rotation-covariant map
from one scalar field to a polar vector has multiplier
`i k F(|k|^2)`. It is longitudinal and is annihilated by `P_T`. At `k=0`, a
constant scalar-to-vector tensor invariant under all proper rotations is zero.
For a pseudoscalar, `q k` is axial under `O(3)`; it is parity-incompatible with
a polar target and remains longitudinal if an axial target is allowed. This is
an exact class-scoped obstruction. It does not cover internal orientations,
anisotropic backgrounds, contact distributions, nonlinear/configuration-
dependent maps, or nontranslation-invariant constructions.

The transverse static kernel

    G^T_ij(x)=(delta_ij+n_i n_j)/(8 pi mu r)

is the correct inverse of the full-pressure symbol `mu |k|^2` on divergence-
free displacement. Translation invariance makes a compact internal source a
stress divergence `f_i=partial_j Sigma_ij`, hence zero net force and at most a
dipolar displacement. The quadratic action with prescribed compensated vector
forces has a reciprocal oriented `1/d` cross term, but no carrier kinetic/KKS
row, constitutive `F_a`, smooth core self-energy, or total carrier-medium
Noether recoil. The supervisor repair therefore states the exact boundary.

### Positive finite-energy escape route outside the compact-moment theorem

The substantive audit found a useful route that the compact multipole theorem
does not cover. For a fixed nonzero vector `a`, define

    u_a(x) = (a cross x)/(1+|x|^2)^(3/2).

It is smooth, divergence free, tangent to every centered sphere, and has
pointwise zero spherical flux. Direct differentiation gives

    curl u_a = { (2-|x|^2)a + 3x(a dot x) }
               /(1+|x|^2)^(5/2).

Thus `u_a=O(r^-2)` and `omega_a=O(r^-3)`; the vorticity is generically not
`L^1`, so the compact/finite-vorticity-moment expansion is inapplicable.
Nevertheless `u_a` is in `H^s(R^3)` for every finite `s>=0`, and

    ||u_a||_2^2 = (pi^2/2)|a|^2.

It is therefore admissible initial data for the standard local Sobolev Euler
class (`s>5/2`), with no all-time assertion. Writing
`phi=(1+|x|^2)^(-1/2)` gives `u_a=-a cross grad phi` and

    phihat(k)=4 pi K_1(|k|)/|k|,
    uhat_a(k) ~ 4 pi i(k cross a)/|k|^2  as k -> 0.

For translated `a` and `b` copies, Parseval and
`F^(-1)(P_T/|k|^2)=(I+n tensor n)/(8 pi d)` give

    rho integral u_a(x) dot u_b(x-d n) dx
      = (2 pi rho/d){a dot b+(a dot n)(b dot n)} + o(d^-1).

This is a genuine smooth finite-energy, local-Euler `d^-1` kinetic escape
route. It is orientation dependent and supplies neither a scalar sign, an
isotropic `q_1 q_2` law, a Gauss source/current, nor P5. Its role is to prevent
extension of the compact-vorticity theorem to all smooth finite-energy Euler
data and to motivate an explicit noncompact-tail route.

**Unit A verdict:** the corrected compact-vorticity multipole, weighted-pressure,
radial-flux, translation-Ward, and local scalar-to-transverse results are
established at their stated domains. They refute the corresponding compact or
local-isotropic scalar Coulomb routes. The target no longer makes a blanket
conclusion for all smooth finite-energy Euler fields and now preserves the
oriented noncompact-tail construction above as a positive escape route.

## Unit B: 0065 topology, helicity, and nonlocal multipliers

### Point puncture

`R^3\{0}` retracts to `S^2`, so its first and second de Rham groups are
`0` and `R`. The closed flux two-form `i_u vol` for `u_q=q x/r^3` has sphere
period `4 pi q`. On every sphere carrying that flux,

    integral_(S_r)|u|^2 dS >= (4 pi q)^2/(4 pi r^2)
                            = 4 pi q^2/r^2.

Coarea therefore gives

    (rho/2) integral_(epsilon<r<R)|u|^2
       >= 2 pi rho q^2(1/epsilon-1/R),

with equality for the radial representative. This is a class-wide lower
bound, not merely one ansatz. Filling the puncture with a smooth source-free
core makes every bounding flux zero. The signed cohomology class and its exact
energy obstruction are established; a finite-energy whole-space Euler charge
is refuted for this fixed-flux route.

### Loop circulation

For a smooth embedded closed loop or finite-radius tube core, the complement
has a meridional `H^1` period. Locally the ideal filament velocity is
`Gamma e_phi/(2 pi s)`, giving logarithmic energy per unit length; a smooth
finite-radius vortex tube removes the core divergence. A compact closed tube
still has zero total vorticity, so its first permitted far field is the impulse
`r^-3` coefficient reviewed in Unit A. When that impulse contraction is
nonzero, pair cross energy is `d^-3`; otherwise it decays faster. The
circulation period is a genuine topological Euler label, but it is not the
point-complement `H^2` flux and supplies no scalar Coulomb interaction.

The target does not confuse this fixed complement period with a new universal
particle charge. Kelvin conservation additionally requires the stated smooth
material-loop Euler evolution and pressure domain; no all-time carrier result
is inferred.

### Helicity and multiplier scope

On the declared decaying divergence-free Hodge domain,

    B = curl(-Delta)^(-1),
    H(omega)=integral B[omega] dot omega,
    delta H=2 integral B[omega] dot delta omega.

The Fourier symbol `i k cross/|k|^2` is Hermitian on the transverse paired
domain, so the factor two and sign are correct. Helicity is conserved for
smooth Euler solutions with the required decay/boundary rows and can encode
linkage, but is already nonlocal and varies continuously under Euler
similarity. Nothing here constructs an integer Hopf sector or makes helicity
the source of an independent Gauss equation.

The nonzero-wavevector classification of linear translation-invariant
isotropic scalar-to-vector multipliers is the same exact representation result
used above and remains valid for declared singular radial multipliers on their
tempered/finite-energy domains. It leaves contact distributions at `k=0`,
nonlinear/configuration-dependent topological maps, internal frames,
anisotropic backgrounds, and separately added scalar/longitudinal sectors
open. The proposed scalar action is therefore an assumption ledger for a new
degree of freedom, not an Euler derivation.

**Unit B verdict:** established as stated at its route-scoped boundaries. The
puncture flux class and universal fixed-flux core-energy divergence, smooth
loop circulation with compact-tube multipole behavior, Hodge/helicity
variation, and nonzero-wavevector linear multiplier obstruction all survive.
They do not establish electric charge, an automatic Hopf/Gauss source, P5, or
a global no-go.

## Oracle audit and findings

The 0064 verifier derives the Newton Hessian, antisymmetric-moment contraction,
displayed dipole coefficients, density-bearing pressure quadrupole, radial
flux, internal-stress zero mode, Oseen scaling, and transverse annihilation.
The 0065 verifier derives the puncture energy/flux powers, transverse scalar
annihilation, absence of a one-wavevector transverse axial source, and
logarithmic ideal-line energy. The final recorded executions exited `0`.

These predicates are sensitive to the principal signs, density, powers, and
zero-mode restrictions, but do not prove the analytic domains or nonvanishing
of a leading coefficient. The target API's isotropic pressure test positively
exposes why `r^-3` is only the first permitted pressure order. Reusing the
unchanged receipts is proportionate; no new numerical or regression run would
answer the discovered scope issue.

| Finding | Direct evidence | Minimum correction / retained result | Upgrade path |
|---|---|---|---|
| The 0064 `strongest_exact_result` overextends compact-vorticity/finite-moment power counting to “direct smooth localized” Euler velocity fields. | The explicit smooth `H^s` zero-flux field `u_a` above has noncompact `O(r^-3)` vorticity and an oriented `d^-1` kinetic cross term. | Restrict the old headline to compact-vorticity, weighted-pressure, radial-flux, and local-isotropic-source classes; retain every exact compact formula and record the noncompact tail as an escape route. | Classify admissible noncompact tails and derive an actual scalar sign/current/action or prove a covering obstruction for that broader class. |
| “Begins/leading” wording does not always state nonzero coefficient hypotheses. | `I=0`, a vanishing orientation contraction, or an isotropic kinetic tensor kills the displayed `r^-3`/`d^-3` coefficient. | Say first permitted coefficient/order, leading when nonzero, otherwise faster. | Prove nonvanishing for a specified carrier family and orientation. |

## Four-axis decision and frontier

- **Verification:** exact analytic verification for both units at the domains
  stated above; recorded symbolic/API evidence is corroborating regression.
- **Review:** Unit A is established after one bounded
  compact/weighted/local-isotropic scope correction; Unit B is established as
  stated.
- **Compatibility:** compatible with constant-density local Sobolev Euler and
  accepted `C-CST-018` only at its prepared transverse scope. No source,
  pressure, parity, density, or sign mismatch remains in the equations.
- **Epistemic:** active route results, not promoted framework claims and not
  P5 completion.

The strongest retained result is twofold. First, compact vorticity gives an
exact impulse-controlled anisotropic multipole interaction, puncture flux has
a sharp class-wide energy obstruction, loop circulation and helicity remain
real internal labels, and one scalar cannot linearly source the homogeneous
transverse channel. Second, smooth finite-energy local Euler data admit an
oriented noncompact `r^-2` tail with a `d^-1` kinetic cross term, so that wider
route remains scientifically live.

The next construction is to determine whether a persistent carrier can support
that noncompact oriented tail with a finite same-field action/current and then
whether an internal doublet or other actual symmetry converts the tensor law
into an isotropic signed observable. Contact zero modes and nonlinear
topological maps remain separate live routes. None licenses a particle,
electric charge, a global charge no-go, or campaign completion.

## Correction check

The single bounded correction is complete and passed. Its receipt has SHA-256
`68be057ffcc5fac71f5bc86f6459b3a46bf6f6b1bf49594802987a77e163fc88`
and pins the following final target hashes:

| Artifact | Corrected SHA-256 |
|---|---|
| `0064/README.md` | `3ea0f3b8f0b1697032f240b36f48e94b3f278a0160780e616e1c841199e34e7d` |
| `0064/derivation.md` | `4e2babcb76e06a0b9d53e67697dd3c861ddeadf36fabbf88ab23fb82a3386add` |
| `0064/result.yaml` | `a118834a38c970581531d54b500866aa4bcf1fc3a7294f107b91068d69cb9232` |
| `0065/derivation.md` | `71d137d5aa14737609d3484781455623e937f59c4ef4a0c61dd0b1ea88fbeca0` |
| `0065/result.yaml` | `9740bb5b86c17bcbaaeddf94a144e084ee2b8e8181e329b99138abd15194ff56` |

The check inspected only the revised scope/nonvanishing sentences and their
direct result edges. The old universal smooth-localized headline is absent;
the first-permitted/nonzero-coefficient qualifications agree across both
targets; and the exact smooth `H^s` oriented `r^-2`/`d^-1` escape route is
retained with its no-charge, no-current, no-persistence, and no-P5 boundaries.
Both corrected result YAML files parse and the affected stale-claim scan and
scoped `git diff --check` pass. The unchanged API/oracles were not rerun and no
second substantive pass was performed.
