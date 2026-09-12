# P253/0106 independent final review of P253/0105

Reviewer: `particle-balance-review`

Target owner: `root`

Activated review README SHA-256:
`a038c3287010b21b4a24a92785689c6682a23261055f6cb52a24e869250246d0`

Frozen target manifest SHA-256:
`5d25796e2cd5b16c9e6c9503c222fb1016a9828cad5a82334a577ab2f854667e`

Post-correction target manifest SHA-256:
`94d27a679d2844c45d2aceedbb30fdd250a476b1a9eae4956174387072c24819`

## Boundary and provenance

This was a content-blind, independent, non-author/non-implementer review of
final P253/0105. The activated README and all frozen target selectors matched,
and `activation-schema.exit` contained exactly `0`. P253/0095 and reserved
P253/0104 were not opened or used. P2, P4, P5, P6, particle and parent
conclusions were not adjudicated.

The mixed-base manifest has sixty entries. At substantive-review opening all
fifty-nine immutable attempt, API, test and receipt entries verified. Its sole
then-mutable central selector,
`proposals/P253-euler-particle-mechanisms/proposal.yaml`, changed from the
historical manifest hash
`bceab0bc3aebda18f00024a85b17c602466f9cffe4349c505f7af0f4f50b56fb`
to
`aac4015dd4170a1ec55621f0613cfceb10492c8b49044c71c754381f68a4b2a2`
when P253/0106 was registered. This expected central-registry drift changes no
0105 scientific artifact or predicate.

After the frozen review artifacts were first written, the shared worktree
changed four claim-bearing 0105 files. Their frozen/post-correction hashes are:

- `derivation.md`: `cef7ec8a...` / `db32c1c4...`;
- `result.yaml`: `e2dc4b44...` / `1033bbe1...`;
- `source-audit.md`: `e9a9ab0e...` / `e187d718...`; and
- `validation.md`: `188c0277...` / `3d20ff26...`.

Their modification times follow the first 0106 review/verdict writes by about
two seconds. They were initially recorded as unconsumed post-review drift. The
coordinator then supplied `0106-bounded-correction-receipt.md` at SHA-256
`e33526e3a1f13fb2a87023072ac3c74d125c7f8482433f69183ca60b6ce174cb`
for exactly one correction-only check. That check inspected only the receipt
and affected Unit-D passages. They now state the same one-normal-derivative
loss independently found by this review, explicitly deny a same-order bounded
inverse, and keep circulation, impulse, center, tag-distribution and all other
finite rows in Unit E. The correction fully closes the precision finding
without changing any equation, API, test, verifier or route verdict. All
sixty-five refreshed manifest rows verify. No second substantive pass or
oracle rerun occurred.

The final exact-v2 receipt records eleven derived algebraic checks with empty
stderr and exit zero. The focused-v2 receipt separately records twenty-one
passing API tests. Static-v3 and repository-validation-v2 record artifact and
workflow agreement. These unchanged receipts were reused without an evidence
rerun. The static-v2 line-wrapping scan failure remains transparent chronology
and has no scientific role.

## A — regular full-core cancellation and compact edge

On a positive full-core component let

    chi=h(P)=zeta/lambda_0,
    -(P_ss+P_s/s)=R^2 lambda_0 h(P).

The matched radial Poisson equations use the same regular integration
constant and give

    Phi_0'=P_s/(R^2 lambda_0 epsilon_EM).

For the first-speed dipole, exact cancellation of
`S_1=h_P v+h Phi_0'/(c_EM^2 R)` makes `h_P` nonzero for `s>0` and, with
`G=h/h_P`, gives

    v=-P_s G/(R^3 lambda_0 epsilon_EM c_EM^2).

Differentiating the profile equation yields
`L_1 P_s=-R^2 lambda_0 h_P P_s`. Substitution into the dipole equation
`L_1v=hP_s/(R epsilon_EM c_EM^2)` cancels every physical scale and leaves

    G_PP P_s^2+G_P(3P_ss+P_s/s)=0,
    partial_s(s P_s^3 G_P)=0.

The mass-density factor belongs to the common nonzero Grad--Shafranov
prefactor multiplying this response; it does not enter the zero condition.
There is no omitted radius, permittivity or wave-speed factor in the reduced
calculus.

At a smooth nondegenerate center, `P_s=-k s+O(s^3)` with `k` nonzero. A
nonzero first-integral constant gives `G_P=O(s^-4)`, then `G=O(s^-2)` and
`v=O(s^-1)`, contradicting the regular dipole condition. Hence `G` is
constant and

    h_P=a h,    h=A exp(aP).

This is a componentwise theorem under `h>0`, `P_s!=0` away from the center and
the matched regular dipole normalization; it is not a theorem through zeros or
a strict tag cutoff. A nontrivial exponential has `h(0)=A!=0`, so it cannot be
continued as the smooth vanishing source at a finite compact Cao free edge.
The compact full-core fixed-ratio cancellation route is therefore refuted by
this mechanism, while strict-band transition moments remain outside the
classification.

**Unit A verdict: established as stated at regular positive full-core scope,
including the scoped compact-edge refutation.**

## B — exact power-law Lane--Emden response

For `zeta=C P^p`, `chi=zeta/lambda_0`, `p>1`, `P_s<0`, and `P(a)=0`, define

    M_j(s)=integral_0^s t P(t)^j dt.

The signed source is

    f=C P^p P_s/(lambda_0 epsilon_EM c_EM^2 R).

Direct integration and the regular/decaying `L_1` Green formula give

    integral_s^a f=-C P^(p+1)/[(p+1)lambda_0 epsilon_EM c_EM^2 R],
    integral_0^s t^2 f
      =C[s^2 P^(p+1)-2M_(p+1)]
        /[(p+1)lambda_0 epsilon_EM c_EM^2 R],
    v=C M_(p+1)/[(p+1)lambda_0 epsilon_EM c_EM^2 R s].

Together with `Phi_0'=-C M_p/(lambda_0 epsilon_EM s)`, this produces

    S_1=C^2 P^(p-1)D(s)
      /(lambda_0^2 epsilon_EM c_EM^2 R s),
    D=pM_(p+1)/(p+1)-P M_p.

All signs and physical factors agree. The square of `lambda_0` makes the
response sign independent of its orientation. Smooth-center expansion gives

    D=-P(0)^(p+1)s^2/[2(p+1)]+O(s^4),

whereas `D` tends to the positive value
`pM_(p+1)(a)/(p+1)` at the edge. Thus `S_1` is negative near the center,
positive near the edge, tends to zero from above at the edge, and has at least
one interior zero. On the analytic positive Lane--Emden core it cannot vanish
*identically* on any open sub-band; this is not a pointwise nonvanishing or
uniform lower-bound claim.

The hypotheses ensure convergence of all finite-radius moments. The formula
does not include the existing strict-band cutoff or its transition
derivatives. The verifier's `P=1-s^2`, `p=6` calculation is correctly labeled
an algebraic sensitivity regression, not a Lane--Emden solution or physical
profile proof.

**Unit B verdict: established as stated for the declared full-core power-law
Lane--Emden profile, with the strict-band response still open.**

## C — fixed-leaf mean-flux obstruction

For a compact area-preserving cross-sectional displacement in the straight
column,

    partial_s(s xi_s)+partial_alpha xi_alpha=0

implies `partial_s(s <xi_s>_alpha)=0`. Center regularity or compact support
sets the constant to zero, and therefore every radial base satisfies

    <delta zeta>_alpha=<-zeta_s xi_s>_alpha=0.

If the tag is a fixed function of vorticity on the same material leaf, its
radial mean vanishes as well.

For a finite torus the exact physical divergence equation is

    partial_r(r xi_r)+partial_z(r xi_z)=0.

The divergence theorem on the region bounded by a regular closed `P` contour
gives `integral_C r xi dot n ds=0`. Taking
`n=grad P/|grad P|` fixes the displayed sign and yields

    integral_C r delta zeta/|grad P| ds=0.

The zero is orientation independent. This is the exact cylindrical-Jacobian
statement, not an asymptotic straight-column replacement. Consequently the
two nonzero radial inner/outer source controls used by the unconstrained
annular matrix are not first-order tangents of this fixed material leaf. The
obstruction does not cover nonradial zero-mean controls, other leaves, or the
unconstructed full steady operator.

**Unit C verdict: established as stated; the radial A/B fixed-leaf lift is
refuted by the weighted contour-flux cokernel.**

## D — raw finite-Cao material tangent

On a compact regular band write the physical meridional measure as

    r dr dz=J(I,alpha)dI d alpha,

with `chi=chi(I)` and `|chi'(I)|>=c_chi>0`. A divergence-free displacement
satisfies

    partial_I(J xi^I)+partial_alpha(J xi^alpha)=0,
    delta chi=-chi' xi^I.

Necessity follows by integrating the divergence equation around a contour:
compact support in `I` sets the constant weighted normal flux to zero, hence

    integral J(I,alpha) T(I,alpha)d alpha=0

for `T=delta chi`. Conversely, set

    xi^I=-T/chi',
    J xi^alpha=-partial_alpha,0^(-1) partial_I(J xi^I),

where the inverse is the periodic zero-mean primitive. The weighted-zero-mean
condition makes the right side periodic; smooth compact support in `I` makes
both components compactly supported. This explicitly realizes
`delta chi=T` and, when `chi=zeta/lambda_0`,
`delta zeta=lambda_0 T`.

The primitive differentiates once in the normal variable. On every compact
regular sub-band it therefore has the correctly typed tame estimate

    ||xi||_(C^k) <= C_k(J,J^(-1),chi'^(-1)) ||T||_(C^(k+1)),

and analogously `H^(s+1)` to `H^s`, with the constant controlled by the
corresponding chart/Jacobian coefficient bounds and `c_chi^(-1)`. P253/0105
claims smooth raw Fréchet surjectivity, not an untyped same-regularity bounded
inverse. No circulation, impulse, center, tag-distribution or other finite-row
restoration is needed for this raw tangent theorem; those rows belong to Unit
E.

**Unit D verdict: established as stated as a raw local material-tangent
characterization, with one normal derivative of loss and no steady-lift
inference.**

## E — full `Y_band x R^2` range

P253/0105 correctly declares

    L_full=(L_resp,L_AB):X_leaf -> Y_band x R^2,
    T_0=(S_Cao,m_0),

and exact cancellation as `T_0 in Ran L_full`. For the stated dual pairing,
annihilation by every element of `ker L_full^*` characterizes membership in
the closure of the range, not literal range membership. A closed-range
estimate, Fredholm theorem or bounded right inverse is needed to bridge that
distinction. Rank two of the finite `L_AB` block is sufficient only for that
finite target, and rank one may suffice for a particular `m_0`; neither fact
decides the full contour-function component.

The source and target norms, common-domain differentiated steady
Euler--Maxwell/free-boundary map, induced Green solves, finite rows, closed
range and pairing with the actual `T_0` are not constructed. The target does
not claim otherwise. This missing theorem has no bearing on Units A--D.

**Unit E verdict: blocked at construction and range analysis of the full
common-domain steady/free-boundary operator.**

## F — evidence, API and exclusions

The importable helpers exactly encode the finite weighted projector, power-law
response, center coefficient, cancellation residual/first integral, contour
flux pairing and local radial displacement. Their docstrings consistently
exclude continuum mapping properties and a steady/free-boundary lift. The
exact verifier derives the eleven advertised algebraic predicates; it does
not prove Lane--Emden existence, analytic regularity, the periodic primitive's
continuum estimates or `L_full`. The focused suite is correctly treated as API
regression only.

The README, derivation, result, source audit and validation agree that the
existing strict-band response, charged redesigned-tag branch, full steady
range, curvature response, compact defect, P2/P4/P5/P6 and particle
conclusions remain open. No stale exact-surjectivity or physical-evidence claim
was found. The single bounded target correction passed and records the typed
derivative-loss scope already used in the independent Unit-D verdict; all
other scientific verdicts are unchanged.

**Unit F verdict: established as an accurate algebraic/evidence package with
explicit scientific boundaries.**

## Combined verdict and next construction

P253/0105 establishes the regular full-core exponential cancellation
classification and its compact-edge obstruction, the exact sign-changing
power-law Lane--Emden response, the fixed-leaf radial mean-flux obstruction,
and the complete smooth weighted-zero-mean characterization of raw local
material tangents on regular finite-Cao bands. These results are mutually
compatible and remain useful independently of the blocked full operator.

The next construction is the common-domain linearized steady
Euler--Maxwell/free-boundary map on that explicit zero-mean tangent class,
including all finite rows, followed by a closed-range/right-inverse theorem and
pairing of the actual strict-band target with its adjoint cokernel. The strict
tag-cutoff moments should be evaluated in the same physical normalization.
This is an active supplier dependency, not P6 or parent completion.
