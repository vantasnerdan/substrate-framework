# P253/0107 — exact fixed-leaf completion versus whole-core and Hill P2 routes

Owner: `particle-foundations`

Status: **centrally registered and initially activated on README SHA-256
`dc8aa8efb88f214b897dbf95bf42a48551fb54e0e3d521be6884e68e07630e0c`,
then paused for the postactivation momentum-domain correction recorded beside
this README.**  No derivation conclusion, source body, API, test, verifier, or
numerical calculation may proceed until the coordinator schema-replays the
corrected README and `attempts/0107/activation-schema.exit` again contains
exactly `0`.  The original activation receipts remain immutable provenance.
Central proposal files, commits, integration, and review remain coordinator
owned.  This attempt will own only files under `attempts/0107` unless a later
activated contract explicitly licenses an importable helper.

## Frozen authority and reviewed boundary

The authority boundary is accepted release state plus the following reviewed
P253 artifacts at their exact scopes:

- byte-current final independent P253/0104 review, SHA-256
  `f4cfcf3284037641e6696d6fbcc86ce53bde613b59eaf268f0aa8ae1c879dbe6`;
- byte-current final P253/0104 verdicts, SHA-256
  `d6224fc58cfe6bf9f044253b6b1a8bdd2b1b02dc4b7ef819803863ed59eb8f6e`;
- byte-current P253/0095 manifest, SHA-256
  `22f0aeffa4c0b3a74b64a266eb3185beb066fc1178cdebd171c2883bcdc9f5c6`;
- P253/0095 Unit-G correction receipt, SHA-256
  `1e6afd6f20c9d41e850f1a3d8c42bab06226d1de74f9c71ba7c9817b37c82a5b`;
- P253/0095 checkpoint whitespace receipt, SHA-256
  `c67e2653b3867837adf04fc8bea72ae6605d75d44dad93b029ef4142a63bc338`.

The checkpoint cleanup changes no scientific content or oracle result.  Its
receipt and the byte-current manifest preserve provenance without altering
the reviewed scientific boundary at checkpoint `68012e81`.

P253/0104 establishes Units A--F of 0095: for each fixed sufficiently thin
Cao member and sufficiently small signed nonzero charge satisfying its
fixed-member gap inequality, an exact returned dynamically accessible phase
carries a real expanding multiplier into the fixed observed-propagator bound

    ||J_loc S_g(jT_*)||_ess,
      X_in^s -> X_DA,loc^s/sym >= c_obs lambda_+^j.       (1)

It does not establish exact nonzero-charge reciprocity, an exact curve in the
corrected fixed-`J_ren` leaf, the nonlinear uniform-Lipschitz
contradiction, a delta-uniform charge wedge, or P2.  P253/0104 separately
retains the charged-column classifier and the full-three-dimensional Hill
bulk/interface problem as active routes.

## Parent objective, attempt objective, and invariants

The unchanged parent objective P253/P2 is one actual localized finite-energy
Euler excitation with a derived restoring mechanism and an open nonlinear
same-leaf neighborhood that persists for all time modulo only its physical
symmetries.  This attempt does not replace P2 with a prepared history, a
linear packet, a finite-time bound, or a finite-codimension manifold.

The exact attempt objective is a competition among three routes:

- **Route A — exact fixed-leaf completion:** close 0095 Unit G by deriving the
  renormalized translation-momentum domain, constructing the carrier-specific
  finite-row witnesses and their bounded lower-triangular inverse, producing
  the exact nonlinear fixed-`J_ren` curve with its auxiliary initial
  carrier/slice normalization, and only then applying the reviewed same-`H^s`
  differentiation bridge to (1);
- **Route B — whole-core charged column:** construct the charged-column/free-
  boundary family and classify every regular-core point by `K(s,tau)` and
  `b(s,tau)`, including the accessible `x->0` graph/resolvent sector; and
- **Route C — materially different Hill carrier:** derive the complete
  bulk/interface/exterior constrained dynamically accessible Hessian for
  Hill's vortex after its actual translation quotient.

Each route receives exactly one verdict: `established as stated`, `refuted`
with its mechanism, or `blocked` with its missing construction.  A route
verdict does not decide a sibling, P2, a particle, or the parent campaign.

All routes preserve incompressibility, the whole-space Hodge/Leray field,
finite energy, the exact coadjoint/tag leaf, charge and Gauss constraints,
physical density and electromagnetic constants, Euclidean stabilizers, and
the reviewed fixed observed operator.  No P4/P5/P6, quantum, electron,
neutrino, particle, or parent-completion claim is in scope.

No production numerical remainder is frozen.  The planned work is exact
Fourier/multipole calculus, functional analysis, finite-dimensional row
algebra, and elliptic/free-boundary analysis; `small-ratio-numerics` is not
activated by this README.  If a later soft eigenvalue or stability-edge
verifier becomes necessary, its numerical design must first consume that
skill in a new activated attempt or append-only expansion.

## Common conventions and fixed member

Use the Fourier convention

    qhat(k)=integral_R3 exp(-i k dot x)q(x)dx,
    uhat(k)=i k cross qhat(k)/|k|^2,                     (2)

and hydrodynamic impulse

    I_h(q)=(rho_m/2)integral_R3 x cross q(x)dx.          (3)

Route A is performed at one fixed sufficiently thin **nonzero-charge** Cao
member `Z_g=(omega_g,chi_g,E_g,B_g)`, with fixed `delta` and fixed
`0<|g|<g_0(delta)`.  All inverse constants, correction radii, and nonlinear
path radii may depend on this member, on `g`, and on `delta`.  Since
`B_g=O(g)`, no `g`-uniform electromagnetic witness inverse is claimed.  At
`g=0`, electromagnetic momentum correction is a separate automatic-zero
case, not the limit of a uniform inverse.

For differences `v=u-u_g`, `q=omega-omega_g`, `e=E-E_g`, and `b=B-B_g`, first
separate the velocity-integral diagnostic

    P_vel,rel(Z;Z_g)
      =integral_R3 rho_m v dx
       +epsilon_EM integral_R3(E cross B-E_g cross B_g)dx
      =integral_R3 rho_m v dx
       +epsilon_EM integral_R3
          (e cross B_g+E_g cross b+e cross b)dx.         (4a)

On the divergence-free `L1` velocity-difference domain, (9) below makes
`integral v=0`; this domain diagnostic is not the robust joint translation
momentum.  Generic Lorentz exchange can change hydrodynamic impulse and create
an `O(|x|^-3)` velocity tail, so the `L1` slice need not be invariant.

Instead define

    Delta I_h=(rho_m/2)integral_R3 x cross q dx,
    Delta P_EM=epsilon_EM integral_R3
          (e cross B_g+E_g cross b+e cross b)dx,
    J_ren(Z;Z_g)=Delta I_h+Delta P_EM.                    (4b)

Equation (4b) is the finite relative translation moment map to be frozen.  If
`f=rho_q E+j cross B` is the physical Lorentz-force density, the route must
derive, with all signs and density factors,

    d I_h/dt=integral_R3 f dx,
    d P_EM/dt=-integral_R3 f dx,                          (4c)

and hence `dJ_ren/dt=0`.  The proof is over every finite time interval from the
**actual coupled Euler--Maxwell equations**: subtract the carrier balance on
expanding balls, display the exact fluid Lorentz-force/current exchange
cancellation against the Maxwell source, and then take the expanding-ball
limit.  The relative stress flux at infinity must vanish in the declared
decay class.  The proof must include the fluid impulse boundary term, Maxwell
stress, the affine Coulomb base, zero difference-charge monopole, and all
cross terms.  Initial `L1/L2` typing alone does not establish a conserved row.
Hydrodynamic impulse may still be imposed as an initial branch/slice row, but
only `J_ren`, not `I_h` separately, is dynamically conserved unless a zero-net-
Lorentz-force theorem is proved.  Thus the invariant nonlinear leaf fixes the
coadjoint/tag Casimirs, charge/Gauss rows, and `J_ren`; `Delta I_h=0` and the
center rows below select the initial carrier representative and modulation
slice and are not advertised as additional coupled conserved quantities.

## Route A — exact fixed-leaf completion

### A1. Fourier moments, impulse, and the `L1` velocity domain

Take compact smooth divergence-free vorticity difference `q` with
`qhat(0)=integral q=0` and finite second moment.  Define

    M_(j i)=integral x_j q_i(x)dx.                       (5)

The exact Taylor expansion and Fourier divergence row must give

    qhat_i(k)=-i k_j M_(j i)+O(|k|^2),
    k_i qhat_i(k)=0
      implies M_(j i)+M_(i j)=0,                        (6)

and the antisymmetric tensor must be identified with (3), including sign and
density:

    M_(m n)=rho_m^(-1) epsilon_(l m n) I_h,l(q),
    qhat(k)=i rho_m^(-1) k cross I_h(q)+O(|k|^2).        (7)

Thus `I_h(q)=0` implies `qhat=O(|k|^2)` and, by (2),
`uhat=O(|k|)`.  The physical-space multipole proof must retain the exact
Biot--Savart kernel and show

    u(x)=O(|x|^-4),       u in L1(R3).                   (8)

Conversely, prove the useful exact lemma: if `u in L1(R3)` and `div u=0`, then
`uhat` is continuous and `k dot uhat(k)=0`; taking `k=t e` to zero for every
unit vector `e` yields

    integral_R3 u dx=uhat(0)=0.                         (9)

Equations (7)--(9) make the velocity integral in (4a) identically zero
**while** exact vorticity mean and impulse cancellation keep the difference
velocity in `L1`.  They do not make that slice invariant and do not replace
the conserved row (4b).  The nonlinear construction is nevertheless staged:
first solve the initial impulse/center rows, then evaluate `J_ren` on that
finite-impulse domain, and only then solve its electromagnetic residual.  A
single open-map IFT that evaluates (4a) on non-`L1` intermediate directions is
not admissible.

### A2. Exact relative electromagnetic momentum and fixed-member witnesses

After the first stage has imposed `Delta I_h=0`, the residual in the conserved
row (4b) is exactly `Delta P_EM`.  Prove its field terms are `L1` by the actual
affine-base and difference estimates, for example by `L2 x L2 -> L1`,
including the Coulomb base and zero-monopole difference.  The exact nonlinear
field correction should use, if the carrier permits it,

    b(a)=0,
    e(a)=sum_(j=1)^3 a_j e_T,j,
    div e_T,j=0,
    e_T,j in C_c^infinity(R3;R3),                       (10)

so Gauss and `div B` are unchanged and the correction to (4b) is exactly
linear:

    L_EM(e_T)=epsilon_EM integral e_T cross B_g dx.     (11)

No `e cross b` term is then generated.  If (10) cannot span the target, the
attempt may add compact divergence-free `b_j`, but must then retain
`E_g cross b`, the nonlinear `e cross b` term, and its IFT remainder rather
than silently linearize it away.

The preferred onto proof must be constructive.  If a covector `y` annihilates
(11) for every compact divergence-free `e_T`, then

    P_L(B_g cross y)=0,
    curl(B_g cross y)=(y dot grad)B_g=0.                 (12)

Use the actual fixed Cao magnetic field and its support/decay to decide the
kernel of (12), not genericity.  If the kernel is zero, exhibit three compact
transverse fields `e_T,j` and the physical matrix

    G_(a j)=epsilon_EM integral
                (e_T,j cross B_g)_a dx,                 (13)

prove `det G!=0`, and record the fixed-member bound

    ||G^(-1)||<=C_EM(g,delta,Z_g).                       (14)

No bound uniform as `g->0` is permitted.  If (12) has a nonzero kernel, name
it exactly and execute the `b_j` alternative or a relative-momentum-domain
reformulation in the same run.  At `g=0`, (11) is identically zero and the
uncharged electromagnetic momentum row is automatic.

Finally, prove conservation of (4b), not only solvability at `t=0`, and show
that the compact corrections remain in the same Gauss/divergence and
fixed-charge domain.

### A3. Initial hydrodynamic-impulse slice, center, and actual parameter rows

For a compact divergence-free displacement generator `xi`, derive with the
adopted bracket sign

    D I_h[ C_0 xi ]=rho_m integral xi cross omega_g dx. (15)

If `y` annihilates (15) for every such `xi`, the solenoidal projection of
`omega_g cross y` vanishes and hence

    curl(omega_g cross y)=0.                            (16)

Relate (16) to invariance under the constant translation `y`.  Compact
support and nontriviality of the Cao vorticity should exclude a nonzero
translation stabilizer.  If so, construct three compact divergence-free
`xi_I,j` and the matrix

    H_(a j)=rho_m integral
                 (xi_I,j cross omega_g)_a dx,           (17)

with `det H!=0` and a fixed-member inverse bound.  If a kernel survives, name
the actual translation stabilizer and reduce the target row accordingly.

Construct compact divergence-free core translations `xi_C,j` that equal
constant coordinate vectors on a neighborhood of `supp omega_g`, using a
cutoff plus a controlled Bogovskii correction.  They must shift the chosen
center rows with an explicit nonsingular matrix `C` while leaving impulse
unchanged because

    rho_m integral xi_C,j cross omega_g dx
       =rho_m e_j cross integral omega_g dx=0.           (18)

Orthogonalize `xi_I,j` against the center/slice rows by subtracting the
`xi_C,j`; equation (18) preserves `H`.  List every actual residual
carrier/slice row.  For the reviewed charged Cao branch, circulation and tag
distribution are coadjoint invariants, charge/Gauss are exact, axial rotation
is a stabilizer, and the fixed circulation--impulse bordered uniqueness must
be used to decide whether any independent carrier-parameter row remains.
No anonymous “finite rows” may be added or removed.

### A4. Staged block matrix, exact IFT curve, and nonlinear test

For each fixed circuit count `j`, first choose the compact orbit tube used by
the reviewed observed-propagator construction and then choose its
high-frequency packet index `N`.  Let `q_(j,N)` denote the resulting raw linear
DA/tag/Gauss tangent.  Do not assume that this tangent annihilates impulse,
center/slice, or relative electromagnetic momentum.  Instead, first display,
after the preceding orthogonalization, the carrier-specific impulse/center
matrix

    M_IC=[[H,0],
          [0,C]],                                       (19)

including all actual surviving slice/parameter rows and a bound on
`M_IC^(-1)`.  Apply its explicit finite-rank projection to the **linear
tangent**:

    a_IC(j,N)=-M_IC^(-1) Row_IC(q_(j,N)),
    q_tilde_(j,N)=q_(j,N)+W_IC a_IC(j,N).               (19a)

For every retained finite row, use exact nonzero-character orthogonality where
available and otherwise direct oscillatory integration plus the exact
compact-support moment formulas--not continuity of total momentum on the
completed `X_in^s`--to prove

    Row_IC(q_(j,N))->0,     a_IC(j,N)->0                (19b)

as `N->infinity` for fixed `j` and tube.  Only after (19a) imposes zero initial
impulse and the center/slice rows may (8)--(9) type the velocity diagnostic
(4a) and may the conserved joint-momentum row (4b) be evaluated.

The second linear stage uses (13).  A displacement correction generally also
changes the advected tag and its Gauss-reconstructed longitudinal electric
field, so its derivative of `J_ren` need not vanish.  Write those exact lower-
left blocks as `K_I` and `K_C`.  Since the transverse Maxwell correction does
not alter vorticity, tag, initial impulse, or center, evaluate the residual
only after the first stage and set

    a_EM(j,N)=-G^(-1) Row_J(q_tilde_(j,N)),              (19c)

and add the corresponding fields (10).  Prove directly, by exact character
orthogonality where available and otherwise by the packet's
oscillatory/moment structure on the now typed relative domain, that
`Row_J(q_tilde_(j,N))->0` and hence `a_EM(j,N)->0`.  With rows
`(Delta I_h,center,J_ren)` and columns `(xi_I,xi_C,e_T)`, the full staged
derivative is lower triangular, not block diagonal:

    M_leaf=[[H,0,0],
            [0,C,0],
            [K_I,K_C,G]],
    det M_leaf=det H det C det G.                        (20)

Equation (20) is a fixed-member statement.  Its inverse and the IFT radius may
depend on `g` and `delta`.  The finite-rank correction in (19a)--(19c) changes
the raw packet by a vanishing fixed smooth component.  Prove in the same
`X_in^s -> X_DA,loc^s/sym` topology that its fixed-time propagated observation
is `o_N(1)` and that the corrected tangents remain weak null.  Therefore this
linear projection preserves the fixed-`J_loc` essential lower bound (1); the
tangent feeding (1) is explicitly `q_tilde_(j,N)` together with its field/tag
components, not the uncorrected packet.

Through each corrected tangent construct the exact coadjoint/tag/Gauss path.
The derivative of every frozen row is now zero, so the **remaining nonlinear**
carrier-specific IFT coefficients are `O(tau^2)`.  The construction must show
that exact joint momentum (4b), including any quadratic field term in an
activated alternative to (10), is zero and conserved along this nonlinear
curve by the coupled source cancellation and expanding-ball argument stated
above.

The quantifier order is frozen:

    fix j; choose its tube; choose N sufficiently large;
    construct the corrected linear tangent;
    send tau->0 in the same-H^s pathwise derivative;
    then send N->infinity.                              (20a)

Only after (5)--(20), the actual parameter rows, support/interface
admissibility, and conservation are closed may Route A consume the reviewed
same-`H^s` pathwise differentiation bridge.  It must then test whether the
hypothetical all-time Lipschitz bound differentiates to a bound on the same
fixed observed propagator in (1).  Success establishes only the fixed-member
nonlinear uniform-Lipschitz **refutation** at the reviewed observable; it is
not a one-datum escape theorem, nonlinear instability, or P2.

If any map has a kernel or any relative flux does not close, record that exact
mechanism, execute the appropriate alternative above, and give Route A the
single verdict supported by the completed ladder.

## Route B — charged-column whole-core classification

Retain the exact fixed-delta charged-column slow block from P253/0095,

    C_col(tau)=x[[0,-Z_tau],[R_tau,0]],
    C_col(tau)^2=-x^2 Z_tau R_tau I,
    K(s,tau)=4Z_tau(s)R_tau(s)/Omega_tau(s)^2.           (21)

Construct the missing charged vorticity/streamfunction constitutive row,
compact free boundary, and first-curvature whole-space Hodge/Leray graph jet.
Keep fixed-`delta` and delta-uniform conclusions separate.

Evaluate (21) on the **entire regular core**:

- `Z_tau R_tau<0` is direct column hyperbolicity;
- `Z_tau R_tau>0` and `K>=1` admits first curvature resonance at
  `|x|=1/sqrt(K)`; evaluate the paired `b(s,tau)` at every resonant point;
- a point with `K>=1` and `b!=0` supplies a refutation mechanism, while an
  isolated zero of `b` activates the ordered second harmonic and higher
  harmonics rather than a stability claim; and
- a restoring candidate requires the common positive-Krein signs and uniform
  constants `kappa_0,eta>0` with
  `0<kappa_0<=K(s,tau)<=1-eta` on its declared whole core.

Even the last inequality classifies only curvature-parametric resonance:
`nu=|x|sqrt(ZR)->0` as the accessible wave angle `x->0`.  Route B must prove
that the action-flow graph/accessibility weight controls this sector or
construct the fixed-Casimir zero-frequency complement with a
resolvent-modulation estimate.  Other Maxwell, interface, and nonaxisymmetric
channels remain part of a P2 proof.

## Route C — full-three-dimensional Hill coercivity

Use Hill's exact patch phase space rather than treating a boundary harmonic as
the whole DA leaf.  For a divergence-free displacement `xi`, retain

    delta omega=-[xi,omega_H]1_B
                 +(xi dot n)omega_H delta_(partial B).  (22)

Derive the relative energy--impulse second variation on (22), including:

- the interior bulk vorticity-stretching term;
- the moving interface and its jump/Lagrange rows;
- the whole-space exterior velocity and pressure/Hodge field;
- the exact translation kernel and quotient;
- all nonaxisymmetric sectors, not only boundary VSH modes; and
- a coercive lower bound in a named full-three-dimensional DA norm, or an
  explicit negative/null direction with its physical mechanism.

Choi's axisymmetric no-swirl rearrangement theorem remains a restricted
positive comparator.  It cannot supply (22) or full-three-dimensional
coercivity.  If Hill fails, continue to a smooth Hill-type/generalized-swirl
carrier whose complete constrained Hessian is defined; do not infer a
carrier-class no-go.

## Verification and evidence plan

Analytic work is load bearing.  A small exact replay helper may later check
the Fourier sign in (7), the block determinant in (20), or a relative-momentum
algebraic expansion, but it cannot prove decay, boundary flux cancellation,
the annihilator arguments, witness existence, an IFT, coercivity, or P2.
Reusable definitions belong in an importable framework module only after body
activation and only if they have consumers beyond this attempt; an
attempt-local verifier must import rather than duplicate them.

For Route A the strongest oracle is the explicit carrier-specific witness
package: compact formulas for every `xi_I`, `xi_C`, and electromagnetic field,
the evaluated matrices `H,C,G`, their determinants and inverse bounds, plus a
direct check that the corrected nonlinear curve satisfies every frozen row
and the relative boundary flux.  For Route B it is whole-core analytic
coverage of `K,b` and the `x->0` graph estimate.  For Route C it is the exact
bulk/interface quadratic form with a coercivity proof or an exposing
admissible direction.

No unchanged P253/0095 oracle will be rerun merely for a tally.  Commands,
stdout, stderr, exit status, source hashes, and maximum verdicts will be
captured on the first execution of any newly activated check.

## Frozen verdict and promotion boundaries

- **Route A established as stated** means (4b)--(20), renormalized joint-
  momentum conservation, an exact curve in the invariant fixed-`J_ren` leaf
  with the declared auxiliary initial impulse/center slice, and the same-target
  nonlinear contradiction all close for one fixed nonzero charged Cao member.
- **Route A refuted** requires an explicit kernel, flux, or domain mechanism
  showing that this declared fixed-leaf correction route cannot work; absence
  of a guessed witness is only a blocked route and activates its listed
  alternative.
- **Route B established as stated** means an actual charged whole-core family,
  complete `K,b` coverage, and the accessible zero-frequency complement are
  controlled in the claimed fixed- or uniform-delta scope.
- **Route C established as stated** means one full-three-dimensional Hill or
  named replacement carrier has a coercive constrained DA Hessian modulo its
  actual symmetries.

Even successful Route A establishes a route-scoped nonlinear Lipschitz
obstruction, not P2.  A successful positive Route B or C still needs the
complete nonlinear open-neighborhood/all-time persistence machinery frozen by
P253/0095.  No route licenses P4/P5/P6, charge or action quantization, quantum
measurement, electron/neutrino identification, a particle, or parent-campaign
completion.  The P253 campaign remains active after every local verdict.
