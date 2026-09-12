# P253/0023 independent source, action, and Jacobi review

## Frozen transaction and independence

This review adjudicates attempt `0015` as one coherent source/action unit. The
unit consists of the García--Hassainia--Hmidi relative-periodic patch pair, its
labeled leaf and translating Hamiltonian, the corrected local KKS loop action,
the autonomous-maximizer and raw-action-Hessian obstructions, the full formal
Jacobi equation with phase/translation modes, and the source reduced-filament
period--twist/Jordan relation. A source theorem or algebraic component does not
inherit a verdict independently of the bridges joining this unit.

Central activation is recorded with exit zero. The frozen `0023` README hash is
`217f5e82527168e528443f7c696b1f2ecfd750008dc95adb77e3a7e013497178`.
The accepted base is `v0.183.0` at
`b6fc902a0942d07996f12a81028fbd3f7c909a43`. The reviewer authored neither
`0015` nor its exact verifier and inspected no `0015` scientific body before
activation. Active attempt `0021`, including its contour propagator,
monodromy, and any later result, is excluded from this transaction.

At the substantive-pass boundary, the target hashes were:

| Artifact | SHA-256 |
| --- | --- |
| `0015/README.md` | `2489b734c9fa5e067c2ecc9c6b793d3f034bfe3e47c2e3a05a66b17bb990d045` |
| `0015/action-sign-correction.md` | `a64811d54db5d2c8348e1fa6ab2098a8f886f189290c4067014772c3c1546d80` |
| `0015/construction.md` | `1d146f6e1d68175a3844e0d27ea3acda6bbe41bb85e88af93905d65fb2a38eb6` |
| `0015/source-receipt.yaml` | `c6db13f1ea859ebe2d2cc47d55f034fd7810a10fcb8a1b0b827e3ae4b4ee5474` |
| `0015/verdicts.yaml` | `2facf6023aeff46c0530cffcb8f25202df44a2ed81ed1b418b07c0f6f6099a6c` |
| `0015/validation.yaml` | `89ead1e3610b43e5ea300d521e8d382f36cba82996514053f17419ad550b5c79` |
| `0015/verify_exact_structure.py` | `0e9350ef0b7b0b7c1eca869096395c8b92cd9e429522466eac98c7fe406cf981` |
| `0015/exact-check.stdout` | `04c84f52c49e3d10b666bb3c09cf19a4652bda1a69e40bff8dbd3ddb50086287` |

The bounded correction requested below changes the final reviewed boundary;
its replacement hashes are recorded in the correction-check section.

## Primary-source applicability

The controlling primary body is García--Hassainia--Hmidi,
*Time-periodic leapfrogging vortex rings in the 3D Euler equations*,
`arXiv:2603.21644v1`. The cached PDF hash
`7ba22ef57ba9453dd2631c33b1ad094a35ce2be171426cecce80fba8a5cdfac6`
and extracted-text hash
`f6461cfc68e2aa3b08cadc3409a36f4f0267c5a7c73fdb03e46f9636c121c871`
match the source receipt.

Theorem 1.1 supports the exact global axisymmetric no-swirl solution for fixed
`0<a<b`, `kappa>0`, sufficiently small `epsilon`, and
`lambda` in the asymptotically full-measure Borel set `C_epsilon`. It consists
of two equal-strength patches of potential vorticity `epsilon^-2`, with
sufficiently high Sobolev boundary regularity. Sections 2--4 support the
canonical half-plane variables, half-period label relation, common axial drift,
and the source center system. They do not supply stability of an open set.

The source's slow time is `tau=L t`, `L=abs(log epsilon)`. Corollary 3.1 gives
the common center as `U_epsilon tau+W_epsilon(tau)`, where `W_epsilon` has slow
period `T`. Consequently the uniform physical speed and translating-frame
physical period are exactly

    c_epsilon=L U_epsilon,
    mathcal_T_epsilon=T/L.                              (R1)

Proposition 3.4 exchanges the two radial center coordinates after `T/2`; their
nonzero initial difference proves that the labeled translating-frame solution
is nonconstant. The source periodic contour spaces in Sections 5 and 7--9 have
zero area mode, reversibility/evenness, half-period label reduction, speed
modulation of the sine mode, and separate cosine-mode treatment. Their
Nash--Moser approximate right inverse solves a periodic invariance equation on
that construction space. It is not the full non-reversible initial-value
Jacobi propagator or a Floquet-stability theorem. The `0015` source receipt
states these boundaries correctly.

## Physical normalization, domain, and labeled leaf

Let `rho_m>0` be constant material density, `rho=r^2/2`, and
`q=omega_theta/r`. Since `d^3x=2 pi r dr dz=2 pi d rho dz`, direct cylindrical
integration gives

    K_phys=(rho_m/2) integral_R3 |u|^2 d^3x
          =pi rho_m integral_Pi Psi q d rho dz
          =2 pi rho_m E,

    I_z,phys=(rho_m/2) integral_R3 (x cross omega)_z d^3x
            =2 pi rho_m P,                              (R2)

where `E=(1/2) integral q Psi` and `P=integral rho q` are exactly the reduced
quantities in `0015`. For axisymmetric no-swirl displacement generators
`chi,psi`, the three-dimensional KKS form similarly reduces to

    Omega_phys=-2 pi rho_m integral_Pi q{chi,psi} d rho dz.
                                                               (R3)

Thus the displayed `0015` Hamiltonian, KKS form, and loop action are all
normalized by the same positive constant `2 pi rho_m`. Their equations and
sign verdicts are unchanged, but calling them physical action quantities
without declaring this normalization loses the density and azimuthal factor.
This is the first part of the requested bounded correction.

The source Green kernel fixes the decaying laboratory-frame velocity on
`R^3`; there is no additional finite-energy harmonic state. The constant axial
term in the translating generator is a momentum augmentation, not an ambient
uniform-flow kinetic-energy contribution. Pressure remains the exact one-form
quotient already covered by the full Euler orbit action.

Each patch has exact `(rho,z)` area `pi epsilon^2`: the ellipse map has unit
determinant and the source imposes zero angular mean on `f`. Each label
therefore has circulation `pi`. While supports are disjoint, independent
Hamiltonian boundary variations can be combined into one physical
area-preserving variation by using disjoint neighborhoods. The local product
leaf is therefore a genuine labeled chart of the total Euler leaf, not two
independently evolved fluids.

There is one small extension point. `0015` defines the leaf using compactly
supported area-preserving diffeomorphisms, whereas the literal translating
flow contains a uniform axial translation. The periodic translating-frame
supports have compact union away from the axis, so a compactly supported
Hamiltonian extension agreeing with `Psi-c_epsilon rho` near both supports
produces exactly the same patch path/tangent. Stating that extension closes the
leaf-membership implication without changing the object. This is the second
part of the same bounded correction.

## Corrected KKS sign and critical-loop calculus

With `{a,b}=a_rho b_z-a_z b_rho` and
`delta_chi q_j={q_j,chi_j}`, the reduced KKS convention is

    Omega(delta_chi q,delta_psi q)
      =-sum_j integral q_j{chi_j,psi_j},
    i_X Omega=dH.                                      (R4)

For `H=E-c_epsilon P`, symmetry of the Green kernel yields
`delta H/delta q_j=h=Psi-c_epsilon rho`, and cyclic bracket integration gives

    D H[delta_chi q]
      =-sum_j <chi_j,{q_j,h}>.

Hence `X_H,j={q_j,h}` and this is exactly the source equation after the frame
change. The sign in the retained `action-sign-correction.md` agrees with the
independent reduction and with the corrected `0005/0011` convention.

For a local potential satisfying `-dTheta=Omega`, the relative action is

    A[q;Q]=-integral_Q Omega-integral_0^mathcal_T H(q(t))dt,

and its first variation is

    delta A(eta)=integral[-Omega(eta,q_dot)-dH(eta)]dt.  (R5)

Since `Omega(X_H,eta)=dH(eta)`, equation (R5) vanishes at `q_dot=X_H`.
The published labeled pair is closed after (R1), so it is an exact critical
loop of this local/relative KKS action. A different filling may change the
action by a symplectic period but cannot change the local first variation.
Multiplying (R4)--(R5) and `H` by `2 pi rho_m` restores the physical action and
leaves the same Hamiltonian vector field.

The stabilizer consists of generators whose induced patch tangent vanishes;
quotienting it is necessary for nondegeneracy. No global primitive, global
Clebsch chart, or quantum phase is inferred. At the source's Sobolev patch
regularity the statement is an exact distributional/local contour variational
identity; it is not a claim that a strong global Darboux chart exists on every
vortex-patch orbit.

## Autonomous functional obstruction

At every source phase,

    D H(q_*)[delta_chi q_*]
      =-sum_j <chi_j,partial_t q_j,*>.                   (R6)

The half-period radial exchange proves the orbit is nonconstant. If its state
velocity vanished at one phase, uniqueness for the autonomous source evolution
would make that state an equilibrium and the entire orbit constant. Thus at
least one boundary distribution `partial_t q_j,*` is nonzero at every phase,
and a compactly supported smooth test generator has nonzero pairing in (R6).

This establishes noncriticality on the actual labeled leaf. A differentiable
local or global maximum must be critical, so neither one phase nor the
phase/translation orbit can maximize `E-cP` on that leaf. Conservation only
gives `D H[q_dot]=0`; it does not give `D H=0`. This is a decisive exact
obstruction to transferring Cao's steady-ring autonomous maximizer to this
non-equilibrium periodic pair. It is not an obstruction to the KKS loop being
stationary, to a time-slice extremum on a smaller class, or to another
modulated variational principle.

## Full-loop high-frequency indefiniteness

In a smooth periodic symplectic trivialization along the loop, the second
variation has first-order principal part

    Q(eta)=integral[-Omega_0(eta,eta_dot)-<S(t)eta,eta>]dt.
                                                               (R7)

Choose a genuine non-stabilizer spatial two-plane with `Omega_0(e,f)=1` and
periodic real variations

    eta_plus=e cos(nt)+f sin(nt),
    eta_minus=e cos(nt)-f sin(nt),
    n=2 pi k/mathcal_T.

Then `integral Omega_0(eta_plus,eta_plus_dot)=n mathcal_T` and the minus family
has the opposite value. Frame-connection and spatial-Hessian contributions on
this fixed smooth plane are `O(1)` in `k`; the symplectic terms are
`-n mathcal_T` and `+n mathcal_T`. The infinite-dimensional patch tangent
space permits the plane to avoid the finite phase/translation and conserved-
quantity conditions without leaving the labelwise zero-area leaf. Equivalently,
the same conclusion follows from a periodic smooth symplectic two-frame along
the loop, with its frame derivatives absorbed into the `O(1)` terms.

The raw action Hessian therefore has both signs at arbitrarily high temporal
frequency. This proves strong indefiniteness on the declared full local loop
space and rules out a sign-definite raw-action maximum or positive Lyapunov
deficit. It does not compute the Euler monodromy spectrum, refute spectral or
orbital stability, or rule out a reduced positive energy--momentum norm.

## Jacobi fields and KKS conservation

Linearizing the full two-label translating equation gives

    partial_t eta_j={eta_j,Psi_*-c_epsilon rho}
                    +{q_j,*,Phi_eta},

    Phi_eta=(1/(sqrt(2) pi)) integral G(x,y)(eta_1+eta_2)(y)dy.
                                                               (R8)

The shared `Phi_eta` retains all self and cross interaction. Time
differentiation and axial-translation equivariance give the exact periodic
Jacobi fields

    eta_ph=partial_t q_*,
    eta_z=-partial_z q_*.                                (R9)

They are independent because the phase field changes the relative radial
configuration whereas common axial translation cannot. They are symmetry
eigenvectors at multiplier one, not the positive/negative temporal-frequency
test families and not generalized eigenvectors.

For every smooth realization of (R8), Hamiltonian-flow invariance gives exact
KKS pairing conservation. If a bounded initial-value evolution and period map
exist on a declared full contour space, (R9) produces two unit eigenvectors and
the map is symplectic. `0015` correctly keeps that bounded realization
conditional; neither the formal equation nor the source periodic inverse
establishes it.

## Source period--twist/Jordan relation

For the source's exact planar reduced center system, let `x(t;lambda)` start at
`(lambda,0)` and have period `T(lambda)`. Differentiating
`x(T(lambda);lambda)=x(0;lambda)` gives

    (M_fil-I)v_lambda=-T'(lambda)v_ph,
    (M_fil-I)v_ph=0,                                    (R10)

where `v_lambda=(1,0)` and `v_ph=x_dot(0)` is nonzero and independent of it.
Proposition 3.2 proves uniform positivity of `T_0'` on compact lambda
intervals, while Proposition 3.3 gives
`T'=T_0'+O(abs(log epsilon)^(-1/2))`. Thus `T'>0` for sufficiently small
`epsilon`, and (R10) is a genuine size-two Jordan chain at multiplier one.

The relation is invariant under conversion from slow to physical time: the
factor in `mathcal_T'=T'/L` cancels the factor in the physical phase velocity.
It is exact for the source finite-dimensional filament/center subsystem. The
full Euler patches exist only on a Cantor set of `lambda`, and the source does
not provide a differentiable open-parameter family of full Euler flows or a
full contour return map. `0015` correctly does not transfer this Jordan block
to finite-core Euler monodromy.

## Oracle and sensitivity audit

The primary PDF and direct calculus are the strongest evidence. The recorded
exact verifier is a useful regression: it checks the translating-generator
sign, cyclic bracket identity, corrected local-action sign, opposite
high-frequency signs, and the two-dimensional symplectic Jordan normal form.
Its source and inputs are unchanged at this review boundary, so it was
inspected rather than rerun, as requested.

The verifier is sensitive to a lone KKS/action sign flip and to loss of the
nilpotent Jordan part. It does not prove source hypotheses, construction-space
coverage, patch functional analysis, or `T'>0`; those are supplied by the
primary-source applicability pass and the direct derivations above. The common
physical factor in (R2)--(R3) is absent from the verifier because it cancels
from every predicate. That explains why the run remained green while the
normalization wording required correction.

No numerical spectrum, soft Hessian eigenvalue, small force, or energy
splitting is used. The small-ratio numerical prescriptions do not bind, and a
pass count is not treated as scientific evidence beyond the exact assertions
actually executed.

## Findings and minimum correction

| Finding | Direct evidence | Disposition | Upgrade path |
| --- | --- | --- | --- |
| The displayed KKS/action quantities omit the physical common factor `2 pi rho_m` while being described as actual/physical. | Independent cylindrical reduction (R2)--(R3), consistent with corrected `0005/0011` density conventions. | Minimum bounded normalization correction; no route verdict changes. | Declare the displayed quantities reduced per `2 pi rho_m`, or multiply `E,P,Omega,A` consistently by that factor. |
| The literal translating flow is not compactly supported although the written leaf uses `Diff_c`. | The `-c_epsilon rho` Hamiltonian generates uniform axial translation. | Minimum bounded leaf-extension sentence; no object or verdict changes. | Use a compact Hamiltonian cutoff agreeing near the compact periodic supports, or state the equivalent compactly supported area-preserving extension. |
| The source periodic inverse is not a full IVP/Floquet result. | Source spaces impose reversibility, label delay, zero mode, and first-mode projections; the inverse is for the periodic invariance equation. | Already correctly scoped; no correction. | Construct the full two-label non-reversible initial-value propagator and reduced transverse power bound. |

No contrary evidence weakens the finite-core existence, exact noncriticality,
local critical-loop identity, high-frequency indefiniteness, Jacobi identities,
or filament-only Jordan relation.

## Coherent verdict and remaining dependency

After the bounded normalization and compact-extension repair, the one coherent
`0015` source/action claim is **established as scoped**. The
actual relative-periodic pair is a critical loop of the correctly signed local
KKS action; it is noncritical and hence nonmaximizing for the autonomous
`E-cP`; the raw full-loop action is strongly indefinite; the exact full formal
Jacobi equation has independent phase and axial-translation modes; and the
source reduced-filament return map has the period--twist Jordan chain (R10).

- Verification: `exact_source_applicability_and_calculus_verified`, with the
  unchanged exact script classified as regression evidence.
- Review: `established_after_bounded_correction`.
- Compatibility: `compatible_proposal_evidence`; the physical density/action
  convention is restored by a common positive factor.
- Epistemic: exact at the stated source, local variational, and filament-only
  boundaries; no stability inference.
- Relationship: Route A autonomous maximum and Route B2 raw-action coercivity
  are refuted by named mechanisms inside the coherent result; the restoring
  obligation remains open rather than inheriting those route verdicts.

The remaining dependency at the frozen `0015` boundary is a bounded
full-source-space initial-value contour propagator followed by a reduced
transverse estimate controlling powers of the finite-core patch monodromy and,
for nonlinear stability, a matching nonlinear modulation bootstrap. No result
from excluded attempt `0021` is used to alter that statement.

Nothing here licenses full finite-core Floquet or nonlinear/orbital stability,
an all-time open neighborhood, three-dimensional/swirl stability, a mechanical
force, a particle, physical spin, quantum phase/statistics, an electron or
neutrino, parent completion, or a global no-go.

## Bounded correction check

The author changed only `0015/construction.md`, from
`1d146f6e1d68175a3844e0d27ea3acda6bbe41bb85e88af93905d65fb2a38eb6`
to `90eec03cd419ca92e6796e6e4d3839b6c6d2152ffa2b04ef49048eacd3a5bff9`,
and added the append-only correction receipt
`0015/normalization-compact-extension-correction.md` with SHA-256
`ffb8d67400089e55046643d714c9c375d50460d38257daba5a849a85f6f993f5`.

The bounded check confirms that equations (8), (9), (16), and (18) now state
the reduced/physical relation with the same factor `2 pi rho_m`; consequently
both sides of `i_X Omega=dH` scale together and every sign and route verdict is
unchanged. Equation (7a) uses a compact cutoff equal to one on a neighborhood
of the compact periodic support tube, so its Hamiltonian and gradient agree
with `Psi-c_epsilon rho` on both patches and realize the same source-regular
patch transport in the compactly supported area-preserving class.

The correction receipt accurately records the scientific effect and the
unchanged hashes of `action-sign-correction.md`, the exact verifier, and its
captured output. In accordance with the frozen plan, none of those artifacts
or the primary source was reopened in this correction pass, and the unchanged
oracle was not rerun. Scoped whitespace checking of the corrected source and
receipt passed. The correction is accepted and closes both review findings.
