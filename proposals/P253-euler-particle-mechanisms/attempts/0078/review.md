# P253/0078 independent final review of corrected 0074

## Integrity, independence, and scope

This review used the corrected anti-narrowing contract at README SHA-256
`9ea8dd8d8ec7bd8008e3bda31108bd53a675069b6739c85986ea2314e39972b9`.
Root replayed the schema after that correction; `activation-schema.exit`
contained exactly `0` and stdout reported `WORKFLOW VALID` before any target
body was opened. All corrected target hashes matched at opening.

The earlier activation against README SHA
`0c6d59503ead129a921066aff10176db9b4fefc139d88763c6ba96e612744423`
is superseded postactivation provenance only. No `0074` body was opened under
that broader wording. The pre-correction `d8936320...` / `91a0e36a...`
package and intermediate digests supplied no evidence.

`particle-balance-review` authored or implemented none of `0074`. Active
`0077` and all equal-frequency analyzer/mixing work remained unopened and
excluded. Candidate supplier APIs listed content-blind in the contract were
not named as proof inputs by the final `0074` source audit and therefore
remained unopened. No production numerical run was performed.

The independently reviewed `0066/0073` regular-cell Weyl sequence and
essential inclusion are treated only as compatibility input. They were not
reconstructed or made a new acceptance gate.

## Source applicability

The target accurately delimits its primary sources.

- Cao--Lai--Qin--Zhan--Zou, arXiv:2206.10165v2, cached at SHA-256
  `6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`,
  supplies the translating no-swirl equation, scaled `C^1` profile limit,
  outer free-boundary description, and Lane--Emden nondegeneracy. Direct
  inspection confirms Lemma 3.8 states that the bounded kernel of
  `-Delta-pU_+^(p-1)` is exactly the two translation derivatives. It does not
  supply the full-vector fiber, common graph, converse essential spectrum, or
  Cao Riesz transfer; `0074` derives those.
- Gallay--Smets, arXiv:1805.05064v3, cached at SHA-256
  `081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`,
  uses the full radial half-line enstrophy space. Proposition 2.1 gives a
  transport-plus-compact split for a straight column and `k != 0`; its passive
  exterior and bounded-operator setting are not imported into the compact
  Cao DA closure. Remark 2.5 supports the existence of axisymmetric imaginary
  Kelvin modes but not their curved-ring transfer.
- The reviewed `0048/0053` input supplies the simple `m=0` frequency sequence
  and positive constrained-energy/Krein sign. The final correction properly
  withdraws any identification of its Sturm amplitude with the absolute
  physical orbit Hessian or KKS norm.

No cited source is used beyond its stated operator, domain, or regularity.

## Verdict unit 1 — exact fiber and global Hodge identities

### Differential and physical normalization

The physical harmonic convention

    v=(2*pi)^(-1/2) v_n(r,z) exp(i n theta)

gives `||v||^2=integral |v_n|^2 r dr dz`. Under
`r=R_delta+a_delta y_1`, `z=z_delta+a_delta y_2`, the measure is exactly
`R_delta a_delta^2 q_delta dy`, with `q_delta=1+delta y_1`; the half-density
factor is therefore correct.

Direct cylindrical calculation confirms the scaled gradient, divergence and
curl rows. In particular, the connection term
`+delta q_delta^(-1)v_theta` in the third curl component is required by
`D C=0`. The combinations `v_r+i v_theta`, `v_r-i v_theta`, and `v_z` carry
the exact angular masses `(k+delta)^2`, `(k-delta)^2`, and `k^2`. Hence all
`delta*n_delta` contributions remain in the leading operator when
`k_delta=delta n_delta=O(1)`.

The weighted adjoint identity `D=-G^*` in measure `q_delta dy` gives the stated
Leray formula. The Piola divergence and curl pullbacks have the correct
Jacobian/cofactor orientation. The support edge is correctly treated as a
coefficient interface of the whole-space Hodge solve, not as an artificial
elliptic boundary.

### Compact-source global comparison

Projecting the physical Newton kernel and applying the half-density unitary
leaves `sqrt(q q')/(4*pi)`: the scale factor `R/a=delta^-1` cancels the
`dalpha=delta dt` substitution. The local kernel tends to
`(2*pi)^(-1)K_0(k|y-y'|)`. The proof upgrades this pointwise fact using the
localized resolvent identity, uniform collar coercivity, and separated smooth
commutator kernel, yielding scalar `H^3 -> H^5` and curl/Hodge
`H^3 -> H^4` convergence.

The global tail is not omitted. Integer-`n` periodic endpoint terms cancel;
regular helical axis powers and exact exterior decaying solutions control the
remote axis and infinity. The conclusion is correctly restricted to compact
DA vorticity or pressure sources, not arbitrary exterior forcing.

The final correction closes the only low-frequency vulnerability. For a
compact DA curl, Stokes gives zero vector mean. At general axial frequency,
Fourier divergence yields

    qhat_j(0,k)=-k partial_(xi_j)qhat_z(0,k),

and fixed support controls the first moment. Thus
`|qhat(xi,k)| <= C(|xi|+|k|)||q||_(H3)`, which cancels the precise zero-shift
helical `K_0` monopole. No uncontrolled curl-potential norm is used.

The positive-core inversion proves the fixed-nonzero-fiber ambient DA closure
equals divergence-free `H_0^3` compact vorticity with the declared finite
rows. It never divides by the vanishing edge coefficient. The `p>=6` boundary
jet gives the three zero traces needed for the `H_0^3` extension.

**Unit 1 verdict: established as stated** for exact fiber/Piola/Hodge algebra
and compact-source global comparison. It supplies no arbitrary-exterior
resolvent or nonlinear orbit chart.

## Verdict unit 2 — common-domain action-flow/Piola convergence

A raw Hanzawa perturbation would indeed be unbounded on the anisotropic
transport graph: the displayed `exp(iNI)` sequence makes a transverse
`delta a partial_I` term grow like `delta N`. The final construction does not
use that invalid map.

The physical meridional area form and Hamiltonian identity give canonical
coordinates

    mu_delta=d beta_delta wedge dI_delta,
    W_delta=omega_delta(I_delta) partial_beta_delta.

The parameter-dependent isochore Morse construction extends through the
elliptic center, the normalized flow-time chart covers the regular annulus,
and the nonzero boundary normal derivative supplies the edge chart. Phase
changes on overlaps depend only on `I`, so gluing preserves both the area form
and flow alignment. The same map is extended through a two-sided collar and
to the axis/exterior, then used for both transport and the metric Hodge solve.
No second, merely compatible coordinate extension is substituted.

The regularity ledger is sufficient. For `p>=6`, `(t_+)^p` is `C^5`; nested
interior Schauder bootstrapping gives the `C^(7,alpha)` profile and
`C^(6,alpha)` velocity control needed for a uniformly `C^5` chart with
`C^4` coefficient convergence.

In common action coordinates, naturality of the bracket gives exact principal
alignment. The difference is

    (omega_delta-omega_0) partial_beta q
      -q^I partial_I(omega_delta-omega_0) partial_beta,

with no transverse `partial_I q`. Since `omega` is bounded away from zero,
the column transport graph controls `partial_beta q`; the displayed
difference therefore maps `D_col -> X` with norm tending to zero. The global
Hodge, curvature, carrier, and finite-row pieces are bounded `X -> X` and
converge at the stated rates. This yields the two-sided graph isomorphisms and

    ||(U_delta^-1 A_delta U_delta-A_col(k_delta))q||_X
      <=rho_delta ||q||_(D_col),  rho_delta ->0.

The common-domain bounded-perturbation theorem gives uniform graph bounds.
On a fixed resolvent contour, the relative product has norm below one, so the
resolvent identity yields `X -> D_col` convergence and norm convergence of
the Riesz projectors. This is leading full-operator graph convergence, not a
second-order `HJ2` claim.

The reviewed `0066/0073` singular sequence is compatible with the same
regular-cell chart/domain; no new reconstruction is required.

**Unit 2 verdict: established as stated.** The action-flow/Piola intertwiner,
two-sided graph bounds, `D_col -> X` smallness, and graph-relative
resolvent/projector convergence are supported on the declared fixed massive
fiber.

## Verdict unit 3 — fixed-fiber essential equality and gap

In volume action-angle variables, the streamline-`m` transport block is a
scalar multiplier `-i m omega(I)` plus a square-zero shear. Outside

    Sigma_T={0} union closure{-i m omega(I):m !=0},

the periodic characteristic/monodromy inverse has uniformly controlled
denominators and first three `I` derivatives. The proof separately handles
the elliptic center in physical Cartesian flow coordinates and proves that
the edge flow preserves all three `H_0^3` traces. It therefore maps the
ambient compact DA space into the closed transport domain and supplies the
reverse Fredholm containment.

Together with the already reviewed compatible Weyl inclusion, this gives
transport essential-spectrum equality under the declared Fredholm convention.
The full Hodge/stretching block is compact on the fixed compact-support
`H_0^3` space: the corrected low-frequency row gives a global one-derivative
Hodge gain, multiplication returns to the coefficient support, and Rellich is
compact. Finite rows do not change the Fredholm essential set.

Hence, in the one fixed nonzero toroidal fiber `n_delta`, the actual carrier
has

    sigma_ess(A)={0} union closure{-i m omega(I):m !=0}

and the punctured essential gap `0<|Im lambda|<omega_min`. This is an
essential-spectrum gap only; discrete Kelvin eigenvalues may occupy it. No
uniform statement over other toroidal fibers is needed or inferred.

**Unit 3 verdict: established as stated** at fixed `n_delta` and fixed
carrier, including reverse containment/equality and the punctured essential
gap.

## Verdict unit 4 — every streamline m in the fixed fiber

For large `|m|`, the triangular transport inverse is `O(|m|^-1)` on the
physical `H^3` block. The full div--curl estimate is uniform in `m`: the
helical masses contain `(m+1)^2`, `(m-1)^2`, and `m^2`, while physical
covariant derivatives absorb the apparent stretching factors. The corrected
compact-curl zero moment controls the exceptional zero-shift `m=+/-1`,
`k=0` channel. Thus a Neumann factorization removes every sufficiently large
`|m|` from a common neighborhood of zero.

For each finite `|m|>=2`, the `k=0` zero equation reduces exactly to

    (-Delta-pU_+^(p-1))phi=0.

Cao's bounded-kernel theorem leaves only translations, which live in
`m=+/-1`; common-domain continuity keeps the finitely many other blocks
zero-free for small positive `k`.

The translation blocks are not incorrectly quotiented away. Their exact
`K_1` exterior Dirichlet-to-Neumann expansion produces a
`k^2 log k` Feshbach determinant with nonzero leading coefficient, so both
translations move off zero at every sufficiently small `k>0`. A finite
minimum after the uniform large-`m` estimate gives one genuine
`epsilon_*>0` excluding zero and a punctured neighborhood from all `m !=0`
blocks at the chosen `k_0`.

**Unit 4 verdict: zero and a punctured neighborhood are excluded as stated**
for every nonzero streamline harmonic `m` in the limiting column used by the
fixed-`n_delta` transfer. This is not a global assertion over all toroidal
fibers.

## Verdict unit 5 — two distinct-frequency positive-Krein modes

The reviewed axisymmetric Sturm problem supplies a simple sequence of
nonzero frequencies tending to zero. Two distinct large fixed indices can be
chosen below `epsilon_*/3`; simplicity separates them from each other, while
Unit 4 separates them from every nonzero-`m` point and essential value.
Disjoint contours therefore have complex rank one in the full column, not
merely in a scalar subblock.

Common-domain graph convergence at `k_delta -> k_0` makes the relative
resolvent norm less than one on each contour. The exact resolvent identity
preserves both ranks and transfers the two isolated modes to every
sufficiently thin Cao ring at
`n_delta=nearest(k_0/delta)`.

The Krein conclusion is correctly scoped. The reviewed supplier proves the
positive constrained-energy sign, and a nonzero-frequency column eigenvector
reconstructs a smooth compact DA displacement. Under the common physical
Piola map, carrier, displacement representatives, Jacobian, and normalized
fiber factors converge in the physical KKS integral. On each finite
rank-one range the nonzero sign is therefore locally constant. Definite
Hamiltonian simplicity keeps the transferred pair on the imaginary axis.

The identity

    |Omega_KKS(e_J,conjugate(e_J))|=1/|sigma_J|

is exact only after normalizing by the actual orbit Hessian. The target
properly leaves open the conversion from the source Sturm amplitude to the
absolute density/Fourier/real-complex physical normalization. The two modes
also retain distinct frequencies. They are not an equal-frequency invariant
doublet and do not supply analyzer mixing axes.

**Unit 5 verdict: two simple, isolated, distinct-frequency positive-Krein Cao
modes are transferred as stated.** Absolute physical KKS normalization,
frequency matching, and noncommuting deformation couplings remain open.

## Oracle evidence and final verdict

The final exact verifier independently derives the scalar differential
complex, shifted helical masses, metric expansions, half-density Newton
coefficient, and cylindrical transport/stretching connection signs; a wrong
curl connection sign is explicitly detected. The initial failed execution
was an implementation-only SymPy substitution that did not test a scientific
identity. The corrected and expanded symbolic receipts have exit `0`. These
checks support Unit 1 algebra only, exactly as validation states; Units 2--5
rest on the analytic proof above.

No bounded correction is required. All five final corrected claim units are
established at their declared scopes. The strongest theorem is a genuine
fixed-carrier, fixed-massive-fiber full-vector graph and Fredholm transfer with
two distinct simple positive-Krein modes. It is a powerful linear supplier,
not a nonlinear carrier, equal-frequency analyzer, action-quantized state, or
particle result.

The next construction is the explicitly retained one: calculate the absolute
physical orbit-Hessian/KKS conversion, tune an exact equal-frequency
resonance while preserving disjoint contours, and derive two noncommuting
physical deformation couplings. Those achievements belong outside this
review and may not use excluded `0077` evidence retroactively.
