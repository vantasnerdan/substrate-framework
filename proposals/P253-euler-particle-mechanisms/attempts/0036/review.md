# P253/0036 independent review of the 0034 canonical continuation

## Frozen transaction and independence

This review adjudicates root-owned `0034` as one coherent fixed claim: the
exact regular-strip axisymmetric Clebsch representation, finite-domain
translation-impulse boundary identity, actual `0027` solitary end-jump
asymptotic, full Hodge-coupled solitary linearization and neutral identities,
material Routh reduction, and the Varholm--Wahlén--Walsh source mismatch.
Stability and particle interpretation are excluded.

Central registration is present and `0036/activation-schema.exit` is zero.
The frozen `0036` README hash is
`a22244731b46ef25ef5e6ef41eee0858bdca0c16b65eb3c64109edce35a8201e`;
the activation exit, stdout, and empty stderr hashes are respectively
`9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`,
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
and `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`.
The observed source head is
`1d378f9aebd10a8ff0dd91f26bf15d371f3e243e`; the reviewer authored or
implemented none of the target equations, API, tests, or receipts.

The substantive boundary is:

| Artifact | SHA-256 |
|---|---|
| `0034/README.md` | `552a165559dea55526b3fb49f67cf3bd4a752d3eda24eaa218ae343d5983b500` |
| `0034/canonical-construction.md` before correction | `5a72e4e5b7bf7a83d8d5add8b0f2e87f2d6bf961f4c9a57764f40137b259ce12` |
| `0034/canonical-construction.md` after correction | `6f233ee844fb3d61c63efb4b7281ee04eaf2536f9753d3958ceba544fe8b3512` |
| `0034/review-correction-0036.md` | `16a65d783d01a84121b0303cac27b75d45f3d7f12d97adc7978485196a1fc99d` |
| `0034/source-transfer.md` | `dc20c11e0aaceed08f1939a1c210ee348a624944cba6b654315ed8416d54b74d` |
| `0034/result.yaml` | `9f1b1df7fca625a7a5698d1f971fda477a97719aa4050670be7dd1d1df5235e6` |
| `0034/validation.md` | `ab74b169384c5f59f6843dfe7e4039e9936f3827c4c560629e0ad541897ac09b` |
| `0034/coordination-model-repair.md` | `18062a99050d919897398b690e43e2893441a037e51a46fe4512027c5e2d9511` |
| `src/substrate_framework/euler_axisymmetric_action.py` | `e062c38657bfbc911b1174eea2bce77f27d6801e560aa6814cb3aac0b6db181b` |
| `tests/test_euler_axisymmetric_action.py` | `93ca6f21eea7554b84414a51b533ef3470c272879934a964194488e4cb97fc11` |
| `0034/canonical-api.stdout` | `86068682b5c90937637c78468882882450ba5bdbc7a81e22a8837393bd21c9c3` |
| `0034/full-operator-api.stdout` | `007654ebadaf0b84d5f292ad1acfa8bd7eed5fcd10cd8eafa2061eed7a5704cf` |
| `0034/jump-api.stdout` | `b3d09481c8a1a1f5f642953336e9b973edf6c09ba57f8a583edcde90501117b4` |

All three captured test exits are zero with hash
`9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`;
their stderr files are empty. The append-only receipts record four, five, and
finally six focused tests, with the final run passing in 3.44 seconds. They
were inspected and not rerun. Active `0037` and every later edit were excluded.

The primary source is Varholm, Wahlén, and Walsh,
[*On the stability of solitary water waves with a point
vortex*](https://arxiv.org/abs/1811.08024v2),
`arXiv:1811.08024v2`, PDF SHA-256
`16da491f9e6ed10e2decea12dc37971e04acb8494a6ddb6fd841cd44c80336bd`.
Accepted `0027/0028` was used only for the exact solitary family and its
weighted convergence; reviewed `0030/0035` supplied only their physical action
and accessible-tangent conventions. Neither supplier's proof was reopened.

## Verdict and strongest supported statement

The coherent `0034` claim is **established after one bounded material-domain
correction**.

For the accepted small whole-space axisymmetric Euler solitary family, on any
compact regular label strip where `xi_r>0`:

1. `zeta={xi,beta}` and the displayed physical first-order action reproduce
   the complete axisymmetric Euler swirl and poloidal-vorticity equations with
   the common factor `2*pi*rho_m`, whole-space Hodge recovery, and the correct
   Clebsch gauge.
2. The physical axial impulse on a finite meridional domain equals the
   canonical bulk translation momentum plus the stated oriented boundary
   one-form. The boundary term cannot be discarded without a fixed-leaf or
   vanishing-boundary hypothesis.
3. The mixed-Casimir distribution is the two-end Clebsch jump. For the actual
   `0027` branch it has the controlled positive leading asymptotic in equation
   (13), uniformly on compact regular label annuli, and is therefore exactly
   nonzero there for every sufficiently small nonzero branch parameter.
   Consequently no single equal-zero-end Clebsch chart represents that family.
4. Equation (14) is the full axisymmetric translating-frame Euler
   linearization, including the nonlocal Hodge perturbation. Translation and
   branch differentiation give the exact zero-mode and inhomogeneous
   speed-companion identities (15)--(16). The speed derivative changes the
   mixed-Casimir distribution and is not automatically a fixed-leaf tangent.
5. Fixing the material angular-momentum distribution gives the exact Routh
   reduction with the negative centrifugal potential and positive centrifugal
   acceleration. On the whole-space column the finite object is the corrected
   relative Routhian (19a), with the relative volume constraint and a declared
   convergence domain; the unrenormalized action exists only on finite domains
   or for compact variations.
6. The direct route importing Varholm--Wahlén--Walsh Theorem 2.4 is refuted by
   a precise hypothesis mismatch: the actual accessible traveling Hessian has
   arbitrarily many high-frequency negative directions (and is not a bounded
   one-negative-direction energy Hessian), whereas their Assumption 6 requires
   one simple negative eigenvalue, one simple symmetry zero, and a positive
   remainder separated from zero. This does not imply instability.

These results are exact canonical and operator identities plus one controlled
small-branch asymptotic. They are not an evolution estimate or a stability
theorem for the nonzero solitary wave.

## Canonical equations, signs, and normalization

Use

    dnu=r dr dz,
    {f,g}=(f_r g_z-f_z g_r)/r,
    xi=r u_theta,
    zeta=omega_theta/r,
    psi=K zeta.

The bracket is the Poisson bracket for the meridional area form. With compact
variations, its cyclic integration identity gives

    H_beta={psi,xi},
    H_xi=xi/r^2+{beta,psi}

for

    H=(1/2) integral (zeta K zeta+xi^2/r^2)dnu,
    zeta={xi,beta}.

Varying

    2*pi*rho_m integral dt [integral xi beta_t dnu-H]

therefore yields `xi_t=-H_beta` and `beta_t=H_xi`. Jacobi then gives

    zeta_t+{psi,zeta}={xi,xi/r^2}=2xi xi_z/r^4.

The sign follows consistently from
`Omega=2*pi*rho_m integral d beta wedge d xi` and `i_X Omega=dH`.
The source, kinetic energy, cylindrical measure, and density factors are all
correct. Reversing the centrifugal term fails the first focused test.

The gauge `beta -> beta+f(xi)` leaves `zeta` invariant and transforms the
canonical equation covariantly. This is a regular-strip chart: the flat core
or exterior does not become a nondegenerate canonical coordinate system, and
no global single chart or physical internal angle is inferred.

## Finite-domain impulse and two-end distribution

Because

    {xi,beta}dnu=d xi wedge d beta,

the differential-form identity

    (r^2/2)d xi wedge d beta
      =beta d(r^2/2) wedge d xi-d[(r^2/2)beta d xi]

gives, with counterclockwise boundary orientation,

    I_D=integral_D beta xi_z dnu
        -integral_boundary(D) (r^2/2)beta d xi.

This independently confirms equation (10), including its sign. At `r=0` the
factor `r^2` removes the axis contribution; other radial and axial faces need
their stated support, fixed-leaf, or convergence conditions. When the surface
term is fixed, the bulk momentum generates negative-coordinate translation,
so `H-cP` has the correct moving-frame sign.

In coordinates `(xi,z)` on a regular strip,

    {xi,beta}=(xi_r/r) partial_z beta|xi,
    dnu=(r/xi_r)d xi dz.

It follows exactly that

    J(xi)=beta_+(xi)-beta_-(xi)
         =integral zeta r/xi_r dz,
    integral zeta D(xi)dnu=integral D(xi)J(xi)d xi.

Thus all smooth test functions `D` recover the mixed-Casimir distribution,
not merely one scalar. A common gauge shift changes both ends equally and
cannot remove nonzero `J`.

For the exact solitary labels `xi=L(a)` and streamsurface `r=R(a,z)`, the
accepted frame relation gives

    f(R,z)=c(R^2-a^2)/2,
    zeta R/xi_r=L/(ca)[R R_a/a^2-R_a/R].

Writing `R=a+h`, the first variation is
`2Lh/(ca^3)`; all `h_a` terms cancel. The implicit relation gives
`h=f(a,z)/(ca)+O(f^2+|f f_r|)`. The `0027` convergence
`F_mu -> f_0 A_*` in a weighted `Z^s` space can be taken at arbitrarily high
fixed `s`; on a compact radial annulus its elliptic equation upgrades this to
the radial `C1` convergence needed here. The spatial weight supplies the
`L1` axial control. Since `integral A_* dX=6/beta`, one obtains

    J(L(a))=
      12 L(a)f_0(a)/(beta c_0^2 a^4) mu L_mu
      +o(mu L_mu).

All displayed coefficients are positive. The remainder is uniform only on a
fixed compact regular annulus, exactly as stated. Equation (13) is an
asymptotic with controlled remainder; its consequence `J>0` for sufficiently
small `mu>0` is an exact statement about the actual branch.

## Full linearization and neutral identities

The nonlinear axisymmetric equations in the frame moving at speed `c` are

    zeta_t=-{psi-c r^2/2,zeta}+2xi xi_z/r^4,
    xi_t=-{psi-c r^2/2,xi}.

Their complete first variation at `(psi_c,zeta_c,xi_c)`, with
`delta psi=K eta`, is precisely

    eta_t=-{Psi_c,eta}-{K eta,zeta_c}
            +(2/r^4)partial_z(xi_c chi),
    chi_t=-{Psi_c,chi}-{K eta,xi_c}.

No pressure/Hodge row, background-gradient coupling, or centrifugal factor is
missing. The API returns this full expression but correctly does not claim a
closed generator, exterior estimate, or well-posed evolution by itself.

Differentiating the exact steady residual in `z` gives
`L_c partial_z q_c=0`. Differentiating it in `c` gives
`L_c partial_c q_c=-partial_z q_c`; the sign is fixed by the `+c partial_z`
term in the moving-frame residual. The test checks the translation derivative
against the full residual rather than against a copied zero. The branch
identity is analytic and is not separately exercised by the test.

Because equation (13) varies along the speed branch, `partial_c q_c` generally
changes the mixed-Casimir distribution. It is an inhomogeneous companion on
the full field family, not an unqualified generalized eigenvector of the
fixed-leaf operator. Neither identity establishes semisimplicity, a resolvent,
a spectral gap, or propagation through the `k=0` threshold.

## Material Routh reduction and bounded correction

For an axisymmetric volume-preserving material map

    X=(r,theta_0+vartheta,z),
    r det[d(r,z)/d(a,b)]=a,

the physical kinetic density is
`2*pi*rho_m a(r_t^2+z_t^2+r^2 vartheta_t^2)/2`. The cyclic momentum per fixed
mass label is exactly `xi=r^2 vartheta_t`. Subtracting momentum times angular
rate at fixed `xi` yields

    R=2*pi*rho_m integral a da db
      [r_t^2+z_t^2-xi^2/r^2]/2.

The sign is correct: varying `-xi^2/(2r^2)` produces the outward acceleration
`+xi^2/r^3`, while the volume constraint produces pressure. The Hamiltonian
restores positive swirl energy. This is a genuine material Routh reduction,
not identification of the Clebsch gauge with the material angle.

The initial body overstated the domain by writing this as an absolute
whole-space action. The accepted column has `u_theta=L_infinity/r` for every
axial position, so its absolute kinetic/Routh integrals diverge. The one
authorized correction changes `canonical-construction.md` from SHA-256
`5a72e4e5b7bf7a83d8d5add8b0f2e87f2d6bf961f4c9a57764f40137b259ce12`
to `6f233ee844fb3d61c63efb4b7281ee04eaf2536f9753d3958ceba544fe8b3512`.
It scopes equations (18)--(19) to finite material domains or compact
variations and adds

    R_rel=2*pi*rho_m integral dt integral a da db
      [(r_t^2+z_t^2)/2-xi^2/(2r^2)+xi^2/(2a^2)]

against the reference map `r=a,z=b`, together with the relative volume
constraint and a declared convergence domain. The added term is fixed label
data and changes no equation. It also scopes the first-order action to
fixed-end relative or compact variations. The correction receipt hash is
`16a65d783d01a84121b0303cac27b75d45f3d7f12d97adc7978485196a1fc99d`.
The bounded diff contains exactly this repair and fully closes the finding;
the unchanged API/tests were correctly not rerun.

## Varholm source applicability

The physical application of `arXiv:1811.08024v2` is a two-dimensional
capillary--gravity free-surface Euler system below a graph, with vacuum,
surface tension, gravity, and a submerged point vortex. Those ingredients do
not describe the present three-dimensional same-fluid axisymmetric column.

More decisively, its abstract Theorem 2.4 requires Assumptions 1--6: compatible
energy and well-posedness spaces, a conserved `C3` Hamiltonian/momentum pair, a
closed injective skew-adjoint Poisson map with dense range, a unitary affine
symmetry, a regular bound-state family, and an augmented energy-space Hessian
with exactly one simple negative eigenvalue, one simple symmetry zero, and a
positive remainder bounded away from zero. Conditional stability also uses
`d''(c)>0` and a solution remaining bounded in the stronger space.

The reviewed `0030` traveling Hessian already has negative accessible
directions at arbitrarily high axial frequencies; disjoint packets/frequencies
give more than one negative direction, and the mixed term is not bounded in
the kinetic energy topology. Therefore Assumption 6 is not available. A branch
slope, the neutral identities, or the source's dense-range Poisson treatment
cannot remove that hypothesis. `0034` correctly classifies this one transfer
route as refuted while preserving the Euler wave and all canonical results.
It makes no instability or global no-go inference.

## API predicates, evidence roles, and next dependency

The six focused symbolic tests have appropriate limited roles:

- one derives the Euler vorticity equation from both canonical equations and
  rejects omission/reversal of the swirl source;
- one checks the exact Clebsch gauge covariance;
- one proves the finite-domain impulse one-form identity and evaluates a
  nonzero oriented rectangle that rejects the opposite boundary sign;
- one checks the mixed-Casimir end-jump change of variables and physical
  multiplier on an explicit field;
- one differentiates the full nonlinear moving-frame Euler residual to obtain
  the Hodge-coupled operator and translation identity;
- one derives the exact solitary jump integrand and verifies cancellation of
  `R_a-1` at first order.

They do not prove global chart coverage, the integrated small-branch error,
operator generation, propagation, Routh convergence, Varholm hypotheses, or
stability. Those conclusions were assigned only to the direct analysis above.
No production numerics, small force, soft eigenvalue, or energy splitting is
used.

Verification is exact analytic with symbolic corroboration; independent review
is established after one bounded correction; compatibility is proposal
evidence on the accepted `0027/0028` supplier; epistemic status is exact for
the stated regular-strip, relative-action, asymptotic, and formal full-operator
identities.

The next dependency is a fixed-mixed-Casimir zero-frequency
resolvent/modulation and continuum-complement estimate for the actual operator
(14), retaining the nonlocal Bessel exterior and translation mode, followed by
nonaxisymmetric sectors if required. This review does not inspect active
`0037` and licenses no stability, restoring-force, particle, spin, statistics,
electron/neutrino, or parent-campaign claim.
