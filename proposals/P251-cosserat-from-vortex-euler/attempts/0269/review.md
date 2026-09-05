# Independent review of proposed C-CST-018

Reviewer: Herdr geometry-review pane `w3:p4`, fresh non-author.  The reviewed
promotion statement is `0265/claim-statement.md` at
`f94d4b0e1d83fe545e65fa69a319a067ce6f3e07f2835dda1e1532dba245a046`.
Its joint proof, assembly/chart input, and literal-current consumer are,
respectively,
`76db984528a835f660a1cede898257d4b332d5b3893f2a3e464a9b748ab855bb`,
`ffd50f9e6752abd43059d4b2014a950b69a0044cfb2d1f4858eacb34ef4e3b49`,
and
`f33ba31b0dcfbd9e6c97ad0f9c98ad804055cb59a560fc5b84c004abe8d9d45b`.
The 0269 activation receipt reports 270 accepted claims and exit zero.  No
registry or release change is part of this review.

## Verdict and strongest supported claim

`route_verdict: established as stated at the compact same-field prepared
finite-window continuum scope`

The proof supports the complete statement in `0265/claim-statement.md` with
its explicit quantifiers and exclusions.  There is one fixed smooth
stationary constant-density ensemble made from the exact periodic disjoint
assembly of compact finite-radius Euler rings.  On that same ensemble, for
each sufficiently small nonzero macroscopic `K`, compact time window, finite
time/observation derivative inventory, retained initial-state amplitude,
and positive accuracy request, there is one finite smooth linear Euler/Lin
preparation whose full physical observation, inherited KKS/Jacobi forms,
and literal material currents have the stated incompressible isotropic
Cosserat limit through `o(|K|^2)`.  The field, cell, tag law, core annuli,
measured coefficients, and whole-state symmetry law are fixed before this
preparation diagonal.  The source map may depend on `K`, the window,
derivative inventory, and accuracy, and its acoustic norm may grow like
`|K|^-2`.

This is a prepared linear-response/second-variation theorem.  It is not an
unrestricted Euler dispersion theorem, a single preparation uniform over
all windows or switching scales, an inverse at `K=0`, an all-time statement,
nonlinear finite-amplitude stability, a universal modulus, or equality of
freely prescribed surface tractions.

## One field and one ordered construction

The existential order is consistent.  The 0263 field, independently
reviewed in 0268, supplies one finite `R` compact ring with a smooth flat
exterior and an exact constant-curl inner region.  `assembly-and-chart.md`
periodizes disjoint copies in one fixed cell `L`, keeps positive constant
fluid density in the quiescent complement, and chooses both the observed
tag annulus and a separated invariant auxiliary annulus inside that same
constant-curl core.  Uniform translations and whole-state O(3), reflection,
and time-reversal operations act on the field, cell, tag, source, observer,
and input coordinates together.  They therefore make one homogeneous
isotropic stationary law rather than an average of conclusions from
different backgrounds.

The response transfer uses the exact volume action-angle chart of this
compact field.  Its acoustic frequency `Omega_1` has a strict derivative on
the fixed core annulus, while the optical frequency
`Omega_1-Omega_2` differs by `O(R^-1)`; one sufficiently large fixed `R`
therefore gives a common nonempty band.  The `n=-1` optical denominators and
the mean-zero `n=0` acoustic denominators are separated from zero on the
same annulus.  The constant-curl identity, cohomological inverse, tag
harmonics, and pressure parametrix are consequently instantiated on this
field, not merely copied from the isolated 0211 ring.

The 0250 two-polarization determinant is supported at precisely this
transfer scope.  Its exact finite-carrier rows are

    det Gamma_opt
      =theta_u f_1(f_2-f_1)(G_u S_perp-S_u G_perp)
      =theta_u f_1(f_2-f_1)(G_u R_perp-R_u G_perp).       (1)

The two tag fractions are distinct and positive.  The tangential
polarization has a leading material correction `R_perp`, so the full spin
row is independent rather than being replaced by `S=-i omega G`.  The
leading normalized determinant is nonzero in the straight Bessel core;
smooth compact-core/chart dependence preserves its strict margin for one
large finite `R`, and sufficiently large finite carrier controls the lower
pressure and normal-leakage terms.  Both real quadratures are retained.
The 0257 O(3) reconstruction is applied only after the per-realization
covariance measurement and has the exact condition-one Gram `2I/3`.

The optical and acoustic columns are then rephased cellwise as reviewed in
0267.  Their physical source phase is constant within each ring.  Periodic
image pressure and the ray-wise mean are not discarded; their smooth
high-action pairings enter the later finite error list.  This preserves the
local determinant and histories without asserting exact pressure
decoupling.

## Full finite-`K` inverse and the `|K|^-2` cost

In acoustic/optical order the normalized leading-carrier map is

    M_inf(K,omega)=[[h_T(K,omega) I_T, B_H(K,omega)],
                    [0,                    B_O(omega)]].  (2)

Here `B_O` is the actual full `(theta,G,S)` optical block and `B_H` is kept,
not assumed small.  The lower-left row is zero for the reviewed cellwise
source law: the local acoustic angle and `G` vanish, the reflected law
cancels its intrinsic axial spin, and the only moved-center spin term is a
third orientation moment.  The 0267 centering repair changes only the
reflection-odd chiral local tensor and leaves this conclusion and the even
oblique acoustic gain unchanged.

The complete compact material Fourier integral is even in `K` after the
whole-law average.  On a fixed ray (uniformly over the compact orientation
law),

    h_T(K,omega)=|K|^2 h_2(omega)+O(|K|^4),
    inf_omega |h_2(omega)|>0.                             (3)

The full oblique coefficient is the reviewed
`(4B_zxx+2B_zzz)/15`, not the earlier axis-only expression.  Since `B_O`
has a uniform inverse on a smaller common band, exact triangular inversion
gives

    ||M_inf(K,omega)^-1|| <= C |K|^-2.                   (4)

Solving `B_O` first, subtracting its literal `B_H` output, and then solving
the acoustic row keeps every mixed observation.  In particular an `O(K)`
optical hybrid contribution may require an `O(|K|^-1)` acoustic source and a
unit acoustic history may require `O(|K|^-2)`; both lie within (4).

There is no hidden interchange of the `K` and carrier limits.  For fixed
nonzero `K` and the frozen finite spectral synthesis, let the actual
finite-carrier observation map be
`M_N=M_inf+E_N`.  The arbitrary-order pressure/material diagonal permits

    ||E_N|| < [2 ||M_inf^-1||]^-1,

after all synthesis and source-normalization factors have been included.
Then the Neumann identity gives an actual finite-carrier inverse with norm
at most `2C|K|^-2`.  Equivalently, applying (2)'s inverse first and retaining
`E_N` as an output error gives the same prepared result.  To obtain a final
`o(|K|^2)` observation error after (4), the raw finite-carrier and exact-flow
errors are chosen smaller than `o(|K|^4)` times the remaining finite
synthesis constants.  The recorded arbitrary-order WKB, image-pressure,
material-jet, and sparse-tag estimates permit this choice.  Thus the phrase
“full finite-`K` gain” means full in macroscopic `K`, with all oblique and
mixed rows retained; it does not claim a carrier-independent bounded inverse
or a value at `K=0`.

After these baseline sources are frozen, their possibly large
`|K|^-2` norms and all form cross entries are finite inputs to the 0262
normalizer.  Its scale `L` may depend on those costs.  The separately
supported auxiliary carrier is chosen last, and the reviewed periodic
supplier makes its observation/current tails smaller than the remaining
budget while matching both full forms exactly.  This is why no relation such
as `|K|=N_s^-D` is needed and why large source norms do not become negative
material probabilities.

## Actual measured coefficients

The claimed inertia is the literal cell-density average of the tagged
rotational inertia tensor.  For `r=x-X`,

    I_tag=rho integral f chi (|r|^2 I-r tensor r) dx.

Haar averaging uses
`E[R(r tensor r)R^T]=|r|^2 I/3`; hence

    E[I_tag]/L^3
      =[2rho E[f]/(3L^3)] integral chi |x-X|^2 dx I
      =j I.                                               (5)

This verifies the factor `2/3`, the tag-fraction expectation, centering,
cell-volume normalization, dimensions, finiteness, and strict positivity.
For a rigid tagged angular velocity it is also exactly the coefficient in
the averaged material angular current.  The independent 0250 `G/S` controls
then prescribe the actual full-spin target `S_int=j Phi_t+o(|K|^2)`, and the
0262 pullback matches the canonical mass/phase block with this same `j`.
The integral is therefore identified in both the prepared current and
prepared action; neither identification is inferred from tube volume alone.

Likewise whole-state Haar averaging gives

    a=[1/(3L^3)] integral_cell |u_per|^2 dx>0,            (6)

as the actual isotropic velocity second moment.  The tag gyration length

    ell_tag^2=integral chi|x-X|^2 dx/integral chi dx>0

is fraction-independent.  Choosing a positive actual frequency `nu` in the
common band and

    cT=j nu^2/(4rho)+nu^2 ell_tag^2,
    cL=nu^2 ell_tag^2                                    (7)

gives

    mu=rho a>0,                 alpha=j nu^2/4>0,
    gammaT=j nu^2 ell_tag^2>0, gammaL=j nu^2 ell_tag^2>0. (8)

These are actual, fixed, dimensionally consistent inputs to the prepared
target.  The finite-window controls and form normalization realize that
target.  Equations (5)-(8) do not purport to derive an unprepared Euler
dispersion relation from background moments.

## Action, chart connection, and literal current

The microscopic action coordinate may be chosen to be the physical
velocity-flux action from 0258.  Its nonzero return-map twist and the volume
action-angle Jacobian are retained on the same regular compact core; no
periodic noncontractible tube or analytic continuation across the flat
support boundary is substituted.  The normalizer's auxiliary sources occupy
the separated invariant annulus of that same field and retain the full Bloch
Leray pressure and all baseline/auxiliary cross forms.

Let `E` be the actual initial observation map after adding the auxiliary
columns.  The 0262 affine equations match

    H_phys=E* H_can E,              Omega_phys=E* Omega_can E,  (9)

with

    H_can=diag(K2,M0),
    Omega_can=[[0,M0],[-M0,0]],     M0=diag(rho I_T,jI).

The same `E` occurs in both equations.  It is not replaced by its limiting
map, and its observed auxiliary tails are retained.  The controlled history
map has a bounded inverse in physical state variables even though the source
inverse grows as `|K|^-2`.

For the evolving source family, `Z(t)=B(t)E(t)^-1` and the stationary Euler
Hamiltonian identity `H=-Omega(.,L.)` give

    Omega(Z,Z_t)=-H_Z-Omega_Z E_t E^-1.                  (10)

The finite-window history bounds and (9) imply
`H_Z=H_can+o(|K|^2)`, `Omega_Z=Omega_can+o(|K|^2)`, and
`E_tE^-1=A_can+o(|K|^2)`.  Since
`Omega_can A_can=-H_can`, the two leading terms in (10) cancel.  Fast
unobserved auxiliary motion therefore does not leave an omitted finite chart
connection.  This conclusion uses both the actual history and both matched
forms; action matching alone would not prove it.

Finally, 0262 transfers the already reviewed universal material identities
without changing their representative.  With

    J_H,t=div F_full,
    S_full,t=div N_full-ax F_full,
    Q_ij=q(t)epsilon_ijk U_t,k,

the improved variables are

    S_int=S_full-div Q,       N_int=N_full-Q_t,
    J_int=J_H,                F_int=F_full.               (11)

The full optical/acoustic output targets, rather than (11) by itself, give
`S_int=j Phi_t+o(|K|^2)`.  The source inventory includes the lower-endpoint
current, `q_t` memory, initial angular charge, pressure, continuous ambient
fluid, convective transport, and moving-boundary reaction.  Combining these
literal balances with the cubic branch mismatch yields

    P_T div(F_int-F_can)=o(|K|^2),
    div(N_int-N_can)-ax(F_int-F_can)=o(|K|^2).            (12)

Equation (12) is an operator statement across the retained initial-state
space because the physical observed-state map is boundedly invertible.  Its
integration-by-parts consequence is complete periodic/bulk constitutive
virtual-work equivalence.  On a cut body the surface term remains explicit;
no pointwise canonical traction equality is claimed.

## Minimum correction and promotion boundary

No false or missing load-bearing step was found in the frozen C-CST-018
statement.  The minimum scientific correction is `none`.  The strongest
license is exactly the prepared, nonzero-`K`, compact-time,
linear-response/second-variation claim above.  In particular, “one
preparation” means one finite linear map for the entire retained initial
amplitude space at a fixed stage of the declared diagonal, not a universal
map uniform in accuracy, `K`, or time window.

C-CST-018 is genuinely additive to C-CST-017: it uses compact Euclidean
closed cores and the new same-field response/action/current implication,
while C-CST-017 retains its periodic noncontractible axial-core scope.  This
review supplies the required individual scientific assessment; it does not
itself edit the registry, render documentation, create a release, replay
downstream consumers, or declare issue200 complete before those separate
promotion transactions pass.
