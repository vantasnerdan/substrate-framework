# P253/0084 independent review of final corrected 0080

## Integrity, independence, and correction check

This review executed against README SHA-256
`1389ed9e68cbe63fd5d6d82959f68f07c946542f8991b44d667ee6bc64ce57a0`.
Central schema activation preceded every target-body inspection;
`activation-schema.exit` contained exactly `0`, stdout reported
`WORKFLOW VALID`, stderr was empty, and every frozen target/API/test hash
matched. `particle-balance-review` authored or implemented none of `0080`.
Active `0083`, reserved `0085/0086`, and `0087` remained unopened and
excluded. No production numerical run or unchanged oracle rerun was made.

The single bounded correction passed. Its receipt
`0080/0084-bounded-correction-receipt.md` has SHA-256
`c107a5d2b7d841ef2decd925c43138436b9e8341dffc88cbacc0a8dfcd452e74`.
The corrected `derivation.md`, `source-audit.md`, `result.yaml`,
`validation.md`, and `author-completion-receipt.md` hashes are respectively
`d2823746db0ef12e9a9d3df334031de3e3ac45e3841e402b13bbbb5f905ae479`,
`3920cc313cc54e70d4ceb80d66b2dcfb12d6c3f15a2a4e60be86a474569ac5df`,
`dc84da15d2ebe68cd5a3f9330b5dd9351f8ba3ae1b8715cdb4530cfc0045a5ec`,
`a5460ad8b66d937a727d41631d2caa2eff78cd07ad069f82ab47f66e5721670d`,
and
`4a93be6526c7920454eed58af33d6b18fd50b30dd64098c4aeecff4e75c3957d`.
The refreshed manifest has SHA-256
`6ffb80da10a09c3212c8297911f0ff4c53fc45617e4dffea24f060086115ec0b`
and verifies every listed artifact. It now consistently identifies equation
(17) as the projected-complement HSE inherited as the open theorem from
`0077`, while leaving the different full fixed-`(mu,c)` `A_epsilon` inverse,
the LS full-kernel-to-sharp-Schur bridge, and the differentiated finite-core
Schur sign open. No equation, route, API, test, or oracle predicate changed.

## Primary-source applicability

Cao--Lai--Qin--Zhan--Zou, arXiv:2206.10165v2, is pinned at SHA-256
`6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`.
Direct source inspection supports the target's ledger:

- equations (2.1)--(2.6) define the axisymmetric operator, symmetric Green
  kernel, measure `dnu=r dr dz`, local logarithm, far decay, and compact-source
  estimates;
- equations (3.1)--(3.4) give the exact uncharged steady equation, axis and
  infinity conditions, circulation, concentration, and symmetry class;
- Lemma 3.8 identifies only the two translations in the bounded entire
  Lane--Emden linearized kernel;
- Lemma 3.9 is an a-priori estimate in its stated `W^(-1,q)`, local `L^infinity`
  and potential norms under a concentrated-support hypothesis. It does not
  itself give the target's `C^(2,alpha)` source-space Schauder map or
  surjectivity;
- Corollary 3.12 and Proposition 3.13 give value-level parameter asymptotics
  with nonzero remainders. Auxiliary system (3.36) is not an exact steady
  constraint or a differentiated Schur theorem;
- Appendix C supplies a symmetry-restricted weighted Fredholm construction,
  not the target's full finite-core bordered operator.

Theorem 1.6 is uniqueness within the source's concentrated uncharged class.
It supplies neither derivative injectivity nor a charged branch. The target
does not import either conclusion. Independently reviewed `0077/0081`
supplies the exact prescribed-current Maxwell reduction, regular-band tag
construction, signed-coupling parity, and affine potential domain only at its
reviewed scope; `0080` itself supplies the new nonlinear mapping and bordered
inverse needed for existence.

## Unit A — global compact-source map and regularity

The fixed annuli `K_0 compactly contained in K compactly contained in {r>0}`
remove the axial degeneracy from the local Schauder step. The global
zero-extension space

    X={zeta in C_even^(2,alpha)(R^2_+): supp zeta subset K_0}

is a closed Banach source space, and the strict negative threshold in an
outer collar persists in a small parameter neighborhood. Thus the positive
part does not create support outside `K_0`.

The source regularity is used correctly. Starting from the Cao compact
bounded density, local uniformly elliptic `W^(2,q)` estimates give
`C^(1,alpha)` for the threshold. Its positive-part right-hand side is then
Hölder, and Schauder bootstrapping gives the regularity needed for the base.
For perturbations already in `X`, the Green solve gains two derivatives on
`K`; compact embedding after multiplication by the supported coefficient
gives `a_epsilon G_1:X->X` compact. This is a separate argument from Cao
Lemma 3.9 and does not miscite that lemma as a source-space Schauder theorem.

For `p>=6`, the scalar positive-part map has enough derivatives for

    F_0(zeta,mu,c)=zeta-epsilon_core^(-2)
                       (G_1 zeta-c r^2/2-mu)_+^p

to be `C^2` from `X times R^2` to `X`, including across the moving zero set.
On the tagged regular band, the nonvanishing gradient, fixed flow box, coarea
formula, and two spare Green derivatives support two shape variations of the
enclosed-volume action and hence the stated `P -> I_P -> chi_P` map.

The Maxwell correction does not cross the free boundary: both `chi_P` and
its derivative are supported strictly inside that regular band. In the
uniformly subluminal range, the anisotropic elliptic Green inverse depends
`C^1` on its compact current and on `c`; restricting its affine potentials to
`K_0`, forming the Lorentz products, and multiplying by the interior tag
therefore yields `M:X times R->X`. The exterior `1/r` potential tail affects
the affine field space but not the compact source regularity. It cannot
upgrade the hydrodynamic map beyond `C^1`, and the target does not claim that
it does.

**Unit A verdict: established as stated.** The global hydrodynamic/tag map is
`C^2`, the eliminated Maxwell insertion is `C^1`, and the derivative is
identity minus compact on the declared source space.

## Unit B — projected soft-cell complement

The signs in the uncharged derivative are exact:

    A_epsilon=I-a_epsilon G_1,
    partial_mu F_0=a_epsilon,
    partial_c F_0=a_epsilon r^2/2.

Symmetry of `G_1` in `L^2(dnu)` makes the cutoff potential derivative the
left functional paired with the right density
`q_t=a_epsilon phi_t`. Its normalization denominator is positive, so
`P=q_t ell_t` and `Q=I-P` are bounded finite-rank projections on the scaled
source space. The domain and range complements are the two copies of `QX`
with the stated scaled source norm; the scaled potential norm is an auxiliary
compactness norm, not an unidentified operator codomain.

In a normalized bad sequence, the source equation gives local compactness of
the pulled-back Green potentials. Escape toward the axis or physical exterior
is excluded by the fixed support, exact Green decay, and the identity part of
`A_epsilon`. The bounded local limit solves

    (-Delta-p U_+^(p-1))phi=0.

Cao Lemma 3.8 leaves the two translations; axial evenness removes the axial
one and `ell_t=0` removes the radial one. The exterior maximum-principle row
then removes the remaining normalized mass. This proves the uniform estimate
on `QX`. Since `A_epsilon` is identity minus compact and deleting one domain
and one range direction preserves index zero, the compressed operator is
onto as well as one-to-one.

**Unit B verdict: established as stated.** Equation (17) is the exact
projected-complement HSE. It is not invertibility of the unrestricted full
fixed-parameter `A_epsilon`.

## Unit C — allowed-growth exhaustion and physical matrices

The bordered blow-up correctly enlarges the bounded class. Subtracting the
constant and linear parameter columns allows logarithmic `m=0` and affine
plus decaying `m=1` exterior behavior. Regularity at the center leaves one
interior radial solution in each angular sector. Matching gives the scaling
mode

    v_a=2 U/(p-1)+y dot grad U

in `m=0` and radial translation `partial_1 U` in the even `m=1` sector;
axial translation is odd. For `m>=2`, the permitted exterior part decays and
is bounded, so Lemma 3.8 excludes an additional mode. Constants, independent
affine rows, and source-generated alternatives are subsequently tested by
the parameter columns rather than silently discarded.

Differentiating the Lane--Emden scaling and translating its compact density
gives

    integral q_a=2 Lambda_p/(p-1)=M_a,
    integral q_t=0,
    integral y_1 q_t=-Lambda_p.

Keeping `dnu=(R+s y_1)s^2dy` and
`I_z=pi rho_m integral r^2 zeta dnu`, then making the translation cell exactly
mass zero, yields the independently rederived leading physical rows

    B_I: diag(M_a,-2 pi rho_m R s Lambda_p)+lower order,
    B_R: diag(M_a,-s Lambda_p)+lower order.

The first line is consistent with
`delta I_z=2 pi rho_m R kappa delta R`; the sign difference is the chosen
orientation of `q_t=partial_1(U_+^p)`. Both determinants are nonzero for every
sufficiently small fixed core. The `B_I` and `B_R` matrices remain distinct;
neither inherits the other's finite row.

**Unit C verdict: established as stated.** The full allowed-growth cell list
and both physical finite matrices close independently.

## Unit D — bordered injectivity and surjectivity

Adding two scalar columns and two scalar rows to `I-compact` leaves each
bordered operator Fredholm of index zero. The contradiction normalization

    ||h||_Xepsilon/H_epsilon^p
      +|eta_0|/H_epsilon+|eta_1|/H_epsilon=1

uses the correct source scale: the raw density is
`epsilon_core^(-2)H_epsilon^p`, while the `epsilon_core^2` in the scaled norm
restores `H_epsilon^p`. Thus it includes, rather than loses, sequences whose
parameter coordinates dominate.

The circulation and centered impulse or mean-radius rows remove the scaling
and radial-translation density cells from Unit C. Once those coefficients
vanish, projection of the rescaled first row on `1` and `y_1` gives the
positive diagonal affine matrix

    diag(integral p U_+^(p-1),
         integral p U_+^(p-1)y_1^2).

It removes both normalized parameter columns; the finite change from
`(eta_0,eta_1)` back to `(delta_mu,delta_c)` has determinant `R s !=0` at
each fixed core. This proves injectivity in both profile- and
parameter-dominated regimes. Index zero then gives surjectivity and a bounded
inverse. No fixed-`(mu,c)` inverse is assumed.

**Unit D verdict: established as stated** for both `B_I` and `B_R` on their
respective target rows. `B_I` fixes circulation and physical impulse for the
charged branch; `B_R` fixes circulation and the frozen exact mean radius for
the local external-core path.

## Unit E — finite-window charged continuation

The speed asymptotic grows with `L=log(1/epsilon_core)`, so a fixed Maxwell
speed does not permit an all-small-core assertion. The target instead assumes
the explicit nonempty window

    L_border <= L <= L_EM

and obtains the uniform ellipticity margin
`1-c^2/c_EM^2 >= 1-(1-eta)^2`. This is the correct hypothesis, not a derived
universal hierarchy.

On that window, Unit A gives a `C^1` reduced map and Unit D gives its exact
circulation--impulse bordered derivative. Factoring the Maxwell potentials by
the signed coupling makes the fluid correction proportional to
`tau=g^2`; extending the algebraic parameter through zero permits the usual
IFT, after which either sign of `g` reconstructs the fields. The estimates

    zeta_g-zeta_0, c_g-c_0, chi_g-chi_0=O(g^2),
    E_g,B_g=O(g)

follow. The tag integral is held fixed in each variation, while the physical
charge across the branch is `Q_g=g C_chi`; no fixed nonzero charge across
different `g` is asserted.

The separate `B_R` inverse gives only a local continuous uncharged
`epsilon_core` path at fixed circulation and exact mean radius. Its
value-level scale law does not prove the quantitative thickness coverage or
differentiated transversality excluded with active `0083`.

**Unit E verdict: established as stated** in the declared finite-`L`,
uniformly subluminal transported-tag Euler--Maxwell extension. It is not a
bare-Euler solution, an infinite-frequency limit, stability, or a selected
charge/action theorem.

## Unit F — conditional leading Schur jet

The value-level source identities consistently give

    c=kappa L/(4 pi R)+O(1),
    mu=3 kappa R L/(8 pi)+O(1),
    I_z=pi rho_m kappa R^2+O(epsilon_core^2).

Differentiating only their displayed leading terms gives

    det D_(kappa,R)(mu,c)=-3 kappa L^2/(16 pi^2 R),
    det D_(kappa,R)(kappa,I_z)=2 pi rho_m kappa R,

and hence the formal/source-derived leading triangular coefficient

    S_lead=-32 pi^3 rho_m R^2/(3 L^2).

The module, tests, and exact verifier correctly implement only this algebra,
including the `2*pi*rho_m` impulse normalization and source-column signs.
They do not prove operator convergence or differentiated remainders.

An `O(1)` value remainder cannot be differentiated. Exact finite-core
nonvanishing therefore still needs uniform `C^1` control of the mass,
Pohozaev, scale, center, interface, variance, and Green remainders. Relating
an unrestricted full `A_epsilon` kernel to the sharp Schur matrix additionally
needs the stated LS tangent identity. The local `epsilon`-to-thickness path
likewise remains conditional on a nonzero derivative or a direct two-sided
coverage estimate.

**Unit F verdict: the leading Schur jet is established as conditional exact
algebra; the sharp finite-core sign, differentiated determinant theorem, LS
bridge, and quantitative local path remain blocked on the named
constructions.**

## Oracle scope and strongest verdict

The focused API and verifier are exposing algebraic support. They derive the
leading parameter and moment Jacobians, conditional Schur product, source
column signs, Lane--Emden density moments, physical two-cell determinant, and
positive affine matrix. They cannot establish the Banach mapping, compactness,
allowed-growth exhaustion, Fredholm range, or charged IFT; those are supplied
by the analytic construction reviewed above. The initial failed receipt was
an implementation-level structural-comparison failure and is correctly kept
as provenance rather than scientific counterevidence.

The strongest supported result is substantial: on the declared compact
even `C^(2,alpha)` source space and for every sufficiently small Cao core in
the nonempty finite subluminal window, the projected `Q A Q` HSE and both
physical bordered derivatives are index-zero isomorphisms. The `B_I` inverse
gives an exact local axisymmetric no-swirl charged relative-equilibrium branch
of the transported-tag Euler--Maxwell extension, with `O(g^2)` fluid/tag/speed
changes and `O(g)` finite-energy Maxwell fields. The `B_R` inverse gives the
separate local exact fixed-mean-radius uncharged path. The leading Schur jet is
correct but its exact finite-core sign is not promoted.

The next construction is not another complement review. It is either the
full fixed-parameter/LS tangent bridge plus differentiated Schur remainder if
that determinant becomes load-bearing, or the already identified downstream
charged-carrier achievements: joint constrained Hessian and
three-dimensional spectral/coercive status, quantitative thickness coverage,
and a demonstrated finite harmonic interval below the subluminal ceiling.
Nothing here licenses stability, action or charge selection, analyzer/Born/
reset dynamics, P4/P5 completion, or an electron, neutrino, particle, or
relativistic conclusion.
