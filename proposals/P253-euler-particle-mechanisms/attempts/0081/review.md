# P253/0081 independent review of final 0077

## Integrity, independence, and evidence boundary

This review executed against README SHA-256
`da46c2ce88de5ab439546d806d7f90ed7c10a601bf87b28b36095e1853a6fc1a`.
Central schema activation preceded target opening; `activation-schema.exit`
contained exactly `0` and stdout reported `WORKFLOW VALID`. Every frozen
`0077` hash matched at opening. `particle-balance-review` authored and
implemented none of `0077`. Active `0079` and `0080` remained unopened and
excluded, and no production numerical run was performed.

The one bounded correction changed evidence attribution only. Its receipt is
`0077/review-correction-0081.md`, SHA-256
`25fd574cdda1fde721a0d01e3d281f4ba248d7118cb91c2be33826fd4d732418`.
The corrected `validation.md` and `completion-receipt.md` hashes are
`470673f0643e79296b95bfe019010133f7a39c75eb66c230c8d21f239814b0a4`
and
`7a205fbf00d85b3011d7d5d5a27db58087531a7df63256cf18a6129abec3e4d9`.
The refreshed manifest has SHA-256
`1420af95049d9bf7418bcca15dfcf00e3875f262d6703de1726015064aff6f9d`
and verifies every listed file. The equations, route verdicts, and verifier
were unchanged; the verifier was appropriately not rerun.

## Source applicability

The primary source is Cao--Lai--Qin--Zhan--Zou,
arXiv:2206.10165v2, locally cached at SHA-256
`6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`.
Direct inspection supports the target's source ledger:

- Proposition 1.4 supplies the exact thin polynomial vortex-ring family with
  fixed circulation and prescribed translating-speed scaling.
- Equations (3.1)--(3.4) supply the exact streamfunction equation, axis and
  infinity conditions, circulation row, and concentration class.
- Lemma 3.8 identifies only the two translations in the bounded kernel of the
  limiting Lane--Emden linearization.
- Lemma 3.9 supplies a concentrated a-priori contradiction estimate under its
  own support and norm hypotheses. It is not the declared charged-map
  complement isomorphism.
- Theorem 1.6 is uniqueness within the source's concentrated uncharged class.
  It does not imply derivative surjectivity for a new Euler--Maxwell map.
- Equation (3.36) is an auxiliary asymptotic parameter system. Corollary 3.12
  and Proposition 3.13 relate exact carrier parameters to it with nonzero
  errors; it is neither the exact physical Schur matrix nor evidence that its
  determinant is nonzero.

The source is therefore used within its actual domain. It supplies no
charged branch, Maxwell theorem, tag stabilizer, persistence theorem, or
physical Schur determinant.

## Unit A — exact-period fixed-tail compatibility

Let

    u=U_{-2}+R,             omega=Omega_{-3}+curl R,

with the stated derivative-weighted bounds, uniformly for one period. The
Euler vorticity equation gives the exact endpoint identity

    omega(T)-omega(0)+integral_0^T [u,omega] dt=0.

Exact periodicity cancels the first term. The homogeneous bracket
`[U_{-2},Omega_{-3}]` has degree `-6`; each bracket containing `R` or
`curl R` is `O(r^(-6-gamma))`. Hence the degree-minus-six coefficient is

    T [U_{-2},Omega_{-3}]=0.

The pressure and density introduce no missing row because this is the curl of
unforced incompressible Euler. For a compactly supported tag coupling the
exterior equation is likewise unforced. The result is not a claim about a
noncompact exterior force.

For a relative endpoint, the proof does not derive decay of the endpoint
difference: it correctly takes
`omega(T)-omega(0)=o(r^-6)`, with the corresponding angular/derivative
control, as the additional compatibility hypothesis. A generic translation
or rotation need not satisfy it. Thus the result excludes a faster bounded
exact-period dressing of one fixed inertial `r^-2` tail with nonzero leading
vorticity residual; it does not exclude time-dependent tails, rotating
asymptotic sectors, or all Euler periodic carriers.

**Unit A verdict: established as stated.**

## Unit B — Lorenz-gauge Maxwell and no-swirl reduction

With

    E=-grad phi-partial_t A,  B=curl A,  y=z-ct,

the vector identity

    e_z cross curl A=grad A_z-partial_z A

fixes the sign in

    E+c e_z cross B=-grad(phi-c A_z).

Changing the plus sign to a minus sign leaves uncancelled
`partial_z A` terms, so this is an exposing sign check rather than a
convention-neutral rewriting. The Lorenz condition and the conserved
traveling current give

    L_c phi=g chi/epsilon_EM,
    L_c A=mu_EM g chi u,
    L_c=-Delta+(c^2/c_EM^2)partial_zz.

For `|c|<c_EM`, this is the chosen decaying comoving elliptic Green inverse,
not a retarded time-evolution inverse. Its fundamental solution
`[4*pi*sqrt(a_c r^2+z^2)]^-1` has the correct anisotropic normalization.
Current conservation plus decay removes the residual Lorenz-gauge harmonic
row; the additive scalar gauge constant is fixed by the decay convention.

Direct cylindrical reconstruction confirms all signs and factors in the
no-swirl reduction. In particular,

    E+u cross B=-grad Phi-H grad P,

and its toroidal curl contains the indispensable electric term
`-g grad chi cross grad Phi`. Dividing the physical toroidal curl by `r`
gives

    zeta=h(P)+(g/rho_m)[chi'(P)Phi-chi(P)H].

The radial Ampere primitive has coefficient
`1/mu_EM-epsilon_EM*c^2=a_c/mu_EM`; combining its radial derivative with
Gauss yields the scalar equation (17). Axis regularity and the whole-space
decay row fix the componentwise primitive constant. The Bernoulli identity is

    p+rho_m|W|^2/2+g chi(P)Phi=b(P),  b'=-rho_m h,

with no missing density or pressure sign. Equations (16)--(17), together with
the declared decay and continuity rows, reconstruct the full axisymmetric
Maxwell system. They define an exact reduced steady PDE; they do not establish
a charged solution of it.

**Unit B verdict: established as stated** for the prescribed-current
comoving Green problem and the exact Euler--Maxwell reduction.

## Unit C — action-volume tag and stabilizer

On a compact regular nested band, the enclosed volume action is a smooth
function of the streamline value. Normalizing angular periods gives the
physical volume form `dI d beta d theta`; changing that normalization only
rescales `I` and the profile. A smooth compact profile in this band, extended
by zero across neighborhoods of the center, critical levels, free boundary,
axis-connected exterior, and infinity, is a genuine smooth spatial scalar,
is material under `W`, has finite integral, and can be normalized to one.

Writing the physical toroidal vorticity as
`omega_0=zeta_0(I) partial_theta`, the exact coordinate components of
`[xi,omega_0]=0` are

    -zeta_0 partial_theta xi^I=0,
    -zeta_0 partial_theta xi^beta=0,
    xi^I zeta_0'-zeta_0 partial_theta xi^theta=0.

The first two divisions require `zeta_0 != 0`. Period integration of the last
row then gives `xi^I zeta_0'=0`, so eliminating the radial/action component
also requires `zeta_0' != 0`. Only on their common nonzero regular band does
`xi^I=0`, which implies `xi dot grad chi_0=0`. Tangential `beta` and `theta`
stabilizer freedom remains. No division, locking, or tag persistence is
claimed at a center, critical/flat level, support edge, or exterior region.

**Unit C verdict: established as stated** as a local regular-band stabilizer
theorem. It is not a global diffeomorphism lock or an all-time charge-profile
persistence theorem.

## Unit D — coupling parity and affine exterior domain

The exact equations are equivariant under

    g -> -g,  (phi,A,E,B,H) -> -(phi,A,E,B,H),

with the fluid, speed, streamfunction, and tag unchanged. Factoring
`phi=g hat_phi` and `A=g hat_A` makes the prescribed-current Maxwell problem
independent of the sign of `g`, while its Lorentz insertion into the fluid map
is `g` times an odd field and is therefore a smooth function of
`tau=g^2`. Thus there is no order-`g` fluid row. Odd/even parity of an actual
local branch follows only after the Unit-E isomorphism gives local uniqueness;
the target correctly states its branch estimates conditionally.

For a smooth compact tag/current, the anisotropic Newton convolution gives
potentials with a nonzero `1/r` monopole. The potentials are not generally in
ordinary inhomogeneous `L2`; subtracting the charge/current monopole places
the remainder in the declared decaying weighted class. Their derivatives
have `r^-2` tails and finite exterior field energy. Compact source support
supplies the required finite moments and core regularity. The charge row is
exactly `Q_g=g integral chi=g` after normalization, rather than a fitted
constant. Smooth dependence of this exterior solve inside the full moving
free-boundary Banach map is not separately proved and remains part of HSE.

**Unit D verdict: established at its exact equation, fixed-source Green, and
conditional-parity scope.** No affine whole-space nonlinear branch is earned
by force counting alone.

## Unit E — conjunctive branch boundary

The target writes the correct formal nonlinear rows and the correct sliced
Schur complement. It also expressly does not declare the exact
axis-weighted/affine domain and codomain, prove range closure and index zero,
or establish the differentiability of `P -> I_P -> chi_P`. Cao's limiting
two-translation kernel and Lemma 3.9 identify a plausible proof route but do
not supply that theorem. This missing infinite-dimensional result is exactly
HSE.

Even conditional on HSE, the physical circulation/impulse constraints leave
the independent matrix

    S = D_(mu,c)(kappa,I_z)-R (L_epsilon^perp)^(-1) B.

The leading impulse response
`2*pi*rho_m*r_*kappa*delta r` correctly detects the radial cell, but it does
not evaluate the full finite-core matrix or prove its determinant nonzero.
Neither source uniqueness nor auxiliary equation (3.36) supplies this second
construction. The even forcing/odd axial-translation pairing removes that
one symmetry row only; it does not imply HSE or the radial Schur sign.

The implicit-function theorem requires the conjunction

    HSE and det(S_epsilon) != 0.

Both remain missing. Consequently the exact charged Cao branch, its
`O(g^2)` fluid correction, and persistence are conditional implications, not
existence results.

**Unit E verdict: blocked on both HSE and the physical two-by-two Schur
determinant.** This is a named route frontier, not evidence against the
carrier or the gauge extension.

## Oracle scope and final verdict

The corrected evidence ledger now matches the code. The verifier checks the
Lorentz theta-curl, the source-side `K'=chi` primitive chain rule, scalar
Gauss/Ampere elimination, Bernoulli coefficient, averaged stabilizer
implication, period/degree ledger, leading impulse coefficient, and
`O(g)*O(g)` backreaction order. It does not assert the complete
traveling-potential identity or the left-hand Ampere primitive; those were
checked analytically above. It also cannot prove affine Banach mapping, HSE,
or the Schur determinant.

The strongest supported result is the conjunction of the exact fixed-tail
periodic compatibility theorem, the explicit subluminal axisymmetric
Euler--Maxwell reduction, and regular-band tag locking, together with the
exact `g -> -g` reduced-map symmetry and fixed-source Coulomb domain. The
charged stationary branch remains open on two independent analytic
constructions. Nothing here establishes persistence, P2/P3/P5, action or
charge quantization, a particle, an electron, a neutrino, or a relativistic
full-field theory.

The next dependency is exact and executable: construct HSE on fully declared
axis-weighted affine/decaying free-boundary spaces, then compute the physical
finite-core circulation--impulse Schur determinant on that same domain.
