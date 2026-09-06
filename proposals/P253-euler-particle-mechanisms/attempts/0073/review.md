# P253/0073 independent final review of 0066

## Review boundary and collision provenance

This is the fixed, non-author and non-implementer review preregistered by
`0073/README.md` at SHA-256
`43653f748a86d09c9ab0f7c286fa52fc470e36082d5ddb9a9e8e8cf31ec02277`.
The reviewer owns only this review and its verdict record. The target is
attempt `0066`, owned by `particle-foundations`; root-owned `0072` was never
opened and supplied no definition, evidence, correction, or conclusion.

A shared-worktree collision temporarily replaced the frozen review README and
then left it absent during an interrupted restoration. The authoritative
chronology is `0073/collision-repair.md`, SHA-256
`d02cb4c1f21377daf838708de8066f75adae5f3e6f5c6922162207e34ba25c7c`.
It records the files opened under the apparent first activation, the pause
before any verdict, exact restoration, and the final schema replay. The
interrupted first receipt has SHA-256
`5bdd305e63a627956551c9346719f0fca2ca917120a79fff75ac5424ec95ca04`
and is explicitly superseded. Neither invalid activation supplied scientific
authority or altered the frozen target hashes.

After exact README restoration, the authoritative activation command, stdout,
empty stderr, and exit have SHA-256 values
`0b8bcc78ad3e326535d62e436c7e5f623fef1e2c3c8d7c12ef1a792f07e8e1c8`,
`d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`,
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`,
and
`9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.
The exit file contains exactly `0` and stdout is `WORKFLOW VALID`.

The accepted authority base is release `v0.183.0`. The target-local activation
receipts and unchanged README retain their frozen hashes. One bounded proof
completion was requested during the substantive pass. The owner supplied
`0066/0073-bounded-proof-completion-receipt.md`, SHA-256
`1151da645862db766143ef9ed71a244678257f9f2d3eb600d00a75cc732022b4`.
The correction-only check verified the four requested constructions and their
direct summaries; no second substantive pass or oracle rerun was performed.

Final corrected claim-bearing hashes are:

| Artifact | Final SHA-256 |
|---|---|
| `0066/derivation.md` | `7bebcc7c51c93c6c28d5736730b9c1a07deb155c440942433a8ac55f54be1507` |
| `0066/source-audit.md` | `330ea0954f10d54f06e383407e484f2182ca817719f8da471da2c6ba0f3056d0` |
| `0066/result.yaml` | `ceedfa810234491fa368f10310dacf35f956637a5a88042b8337d8983acad4d1` |
| `0066/validation.md` | `4341e462b5de231e0537b5563f79da242a2edd50e51b3bdf323e2badd543868b` |

## Source applicability

Cao--Lai--Qin--Zhan--Zou, arXiv:2206.10165v2, is pinned at SHA-256
`6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`.
Equations (1.3)--(1.6) supply the axisymmetric no-swirl translating equation;
Proposition 1.4 supplies the polynomial compact ring; Proposition 3.2 supplies
scaled `C^1` convergence; and Lemma A.2, especially (A.5)--(A.8), supplies the
single outer boundary and its nonzero normal derivative. The source does not
state a nonaxisymmetric DA graph, essential-spectrum theorem, curved Kelvin
mode, contour transfer, or nonlinear branch. The additional regularity and
Weyl construction reviewed below are calculations in `0066`, not imported
Cao theorems.

Gallay--Smets, arXiv:1805.05064v3, is pinned at SHA-256
`081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`.
Its Proposition 2.1 proves the multiplication/compact split in its straight-
column enstrophy space. Equations (2.17)--(2.18) and Remark 2.5 supply the
axisymmetric nonzero-axial-wave-number pencil, energy identity, and imaginary
Kelvin modes. They do not transfer their space, compactness, full spectrum, or
isolated modes to a curved Cao ring. The generalized compact-core
Sturm--Liouville asymptotic is separately retained from `0048` at the scope of
the independent `0053` review; the correction now defines its Liouville
length explicitly.

The corrected `0062`/independent `0063` input is used only for the cylindrical
Hodge block, positive-core DA map, KKS density, and the `s=4`, `p>=6` collar
counterexample. It supplies no global inverse, limiting absorption, `q_*`, or
`V_*`. No source named above supplies the high-toroidal Piola/Hodge transfer
or whole-column spectral isolation.

## 1. Interface bootstrap, boundary jet, and regular-cell foliation

For the physical defining function

    P=G_1 zeta-(W/2)r^2 log(1/epsilon)-mu,
    zeta=epsilon^(-2)(P_+)^p,

the axisymmetric operator gives

    -Delta P+r^(-1)partial_r P=r^2 epsilon^(-2)(P_+)^p.

The translating term is harmonic for this operator. Lemma A.2 gives the
nonzero normal derivative on the actual thin-ring free boundary; fixed-
carrier elliptic regularity and normal Taylor expansion then give

    P=kappa_P(s)d+O(d^2),
    zeta=kappa_zeta(s)d^p+O(d^(p+1)),

with positive coefficients and two-sided comparison along each fixed compact
boundary. This is a fixed sufficiently thin carrier statement, not a uniform
raw physical scaling theorem for all family members.

The bounded correction supplies the previously compressed interface
bootstrap. On nested fixed two-dimensional balls, source `C^1` convergence
and the exact rescaled equation give a uniformly bounded full right-hand side.
Interior `W^{2,q}` estimates for `q>2` yield uniform `C^{1,gamma}` with
`gamma=1-2/q`. For `alpha<gamma`, interpolation controls the positive-part
term in `C^{0,alpha}`, and the complete drift
`delta(1+delta*y_1)^(-1) partial_1 w_delta` tends to zero in the same norm.
Interior Schauder therefore gives `C^{2,alpha}` convergence across the
positive-part interface; no `C^1` shortcut is used.

The limiting center Hessian is negative definite. Quantitative persistence
near the center and a uniform lower bound for `|grad U|` on its compact
complement exclude secondary critical points. Together with the single
regular outer boundary, this gives one elliptic center and regular closed
positive levels without a separatrix. The center Hessian and the outer
gradient, cylindrical weight, and curve-length bounds give positive finite
limits at both ends of the period interval. Continuity then yields, for each
fixed sufficiently thin carrier,

    0<omega_min<=omega(I)<=omega_max<infinity.

**Verdict:** established at this fixed-carrier/rescaled scope. No uniform
all-family physical-time bound or arbitrary-Cao-member foliation follows.

## 2. Fredholm convention and singular-sequence bridge

The correction now defines

    X_n^3=closure in H^3(R^3) of compact, fixed-n DA curls

with the declared finite rows, and

    D(T_n)={q in X_n^3:T_n q in X_n^3},
    D(A_n)=D(T_n),  A_n=T_n+K_n.

The equality of domains follows from the bounded compact `K_n` result checked
below. The essential-spectrum convention is explicit:

    sigma_ess^+(A_n)={lambda:A_n-lambda is not left Fredholm}.

For the normalized packets, `q_j` converges weakly to zero in `X_n^3` and
`(A_n-lambda)q_j` converges strongly to zero. Hence
`A_n q_j=lambda q_j+o(1)`, so the sequence is weakly null in the graph domain
as well. If a bounded graph-domain left regularizer existed with

    R(A_n-lambda)=I-C,

where `C` is compact on the graph domain, then both `R(A_n-lambda)q_j` and
`Cq_j` would tend strongly to zero. This would force `q_j` to zero in graph
norm, contradicting its unit `H^3` norm. Thus `A_n-lambda` is not left
Fredholm, and in particular is not Fredholm.

**Verdict:** the construction proves membership in the upper-semi-Fredholm
essential spectrum and therefore in the standard non-Fredholm essential
spectrum. It does not prove equality with the full spectrum, a converse
direct-integral theorem, or absence of interface/separatrix essential
spectrum.

## 3. Fixed-nonzero Hodge compactness

For one fixed integer `n!=0`, all DA vorticities have common compact support
in `K_0`. Divergence freedom and compact zero extension give zero spatial
integral, while the fixed support and `H^3` bound control first moments.
Consequently `q_hat(k)=O(|k|)` at the origin. This removes the low-frequency
singularity of the decaying whole-space Biot--Savart multiplier, while its
high-frequency symbol gains one derivative:

    B_n:H^3(K_0)->H^4(R^3).

Restriction to a fixed neighborhood of `K_0`, followed by
`H^4 compactly embedded in H^3`, is compact. In
`[B_n q,zeta partial_theta]`, the apparent derivative of `B_n q` in the
second bracket leg is only the fixed multiplier `i*n` plus cylindrical
connection terms. It is therefore zero order in this fixed fiber. The `C^4`
carrier coefficients and local Rellich compactness prove that
`K_n:X_n^3->X_n^3` is compact, with the global exterior velocity and decaying
harmonic row retained.

**Verdict:** established for each fixed nonzero `n`. It is neither a uniform
large-`n` estimate nor an `n=0` compactness theorem.

## 4. Canonical cell chart and constrained Weyl packets

On a regular meridional annulus, the corrected action is the enclosed
`r dr dz` area divided by `2*pi`; the travel-time angle `beta` is normalized
to period `2*pi`. With the compatible orientation,

    r dr dz=dI d beta,
    dV=dI d beta d theta.

Thus this is a canonical volume action-angle chart, not an arbitrary
streamline label. Both `m` and `n` are integers. For fixed `n!=0`, regular
`I_0`, `epsilon_j->0`, `N_j->infinity`, and
`N_j*epsilon_j->infinity`, the compact packet polarization

    a_(m,n)=partial_beta-(m/n)partial_theta

is tangent to `I=constant`, exactly divergence free, and kills the unipotent
transport shear. Since `omega_0=zeta(I)partial_theta`, it satisfies exactly

    C_0 xi_j=i*n*zeta*xi_j.

After `H^3` normalization, its amplitude is of order
`(N_j^3 sqrt(epsilon_j))^(-1)` and

    ||(T_n+i*m*omega(I_0))q_j||_H3
      <=C(epsilon_j+N_j^(-1)).

The packets are smooth compact elements of `D(A_n)` and are graph-weakly
null. Any independent continuous finite rows can therefore be removed using
fixed smooth profiles in the same fixed-`n` DA graph sector; the correction is
`o(1)` in graph norm. Compactness of `K_n` then gives, with the harmless sign
absorbed by `m in Z`,

    i*m*omega(I_0) in sigma_ess^+(A_n)

for every regular cell, regular `I_0`, integer `m`, and fixed nonzero `n`.

**Verdict:** the regular-cell constrained DA Weyl/Fredholm essential-spectrum
inclusion is established. This is the strongest new operator theorem in
`0066`.

## 5. Nonlinear weighted quotient versus tangent norm

The corrected `0062` collar sequence retains a displacement quotient bounded
below while its ambient `H^3` vorticity and `H^3` generator-graph norms tend
to zero. It therefore refutes both ordinary unweighted lower bounds at
`s=4`, `p>=6`, fixed `|n|>=2`. It does not refute the positive-core inverse,
the weighted-orbit tangent norm, or a source-specific sandwiched trace.

`0066` defines a separated weighted tangent quotient and proves the exact
smooth-path identity

    omega_1-omega_0=curl integral_0^1(u_t cross omega_t)dt,
    B(omega_1-omega_0)=P_L integral_0^1(u_t cross omega_t)dt.

It also records composition/inversion in the ambient `H^5`
volume-preserving group. It does not prove representative-independent descent
of those operations and Hodge/pressure maps to a completed nonlinear quotient.
The reverse Hardy/interface estimate needed for a two-sided nonlinear chart
is absent and is explicitly named.

**Verdict:** tangent weighted norm and smooth-path Hodge identity established;
ordinary unweighted equivalence refuted at the stated collar scope; nonlinear
weighted quotient/chart remains open.

## 6. Full-vector high-n Piola/Hodge transfer

The scaling `n_delta=nearest(k_0/delta)` correctly keeps
`k_delta=n_delta*delta` near a fixed positive `k_0` and avoids the zero-mass
scalar Hodge singularity. But scalar coercivity and fixed-domain profile
convergence do not control the physical vector operator when
`delta*n_delta=O(1)`. `0066` expressly leaves unproved:

- the full Piola transform of all vector components, curl, divergence, Hodge,
  Leray pressure, density, Jacobian, and cofactors;
- every `delta*n_delta` contribution;
- the moving interface, axis, collar, exterior, and harmonic maps; and
- the common-domain `D(A_col)->X` graph-relative estimate.

Equations (34)--(37) are targets or conditional consequences, not theorems.

**Verdict:** positive-mass scaling and scalar coercivity are retained exact
atoms; the full-vector Piola/Hodge graph transfer is conditional.

## 7. Limiting modes versus whole-column isolation

For the limiting compact column, the axisymmetric pencil becomes

    H_k u=lambda Phi u,
    -(r u')'+(r^(-1)+k^2 r)u=lambda r Phi u.

The correction defines
`L_Phi=integral_0^a sqrt(Phi(r))dr`. The reviewed compact-core
Sturm--Liouville law gives

    sqrt(lambda_j)=pi*j/L_Phi+O(1),
    sigma_j(k)=k/sqrt(lambda_j)->0,

with simple modes inside this `m=0` sector. The Rayleigh identity supplies
positive wave energy there. This establishes the limiting sectoral supply of
low positive-energy Kelvin frequencies.

It does not isolate either proposed circle from the entire constrained
full-column spectrum. A common analytic domain for all sectors, large-`|m|`
resolvent control, exclusion of finitely many remaining sectoral crossings,
and the full constraint/pressure/exterior rows are absent. Without that
isolation and the Unit 6 graph transfer, Riesz rank preservation cannot be
invoked for the Cao ring.

**Verdict:** limiting `m=0` Sturm--Liouville modes and their sectoral positive
energy are established. Whole-column isolation, Cao eigenmodes, their Krein
normalization, the doublet, and equations (38)--(39) are conditional.

## 8. Fixed-domain regularity, uniqueness, and `omega_min`

The nested-ball argument proves the required `C^{2,alpha}` profile convergence
on a fixed domain for the sufficiently thin family. Negative-definite center
Hessian persistence gives a unique center locally, the off-center gradient
bound excludes any other critical point, and the single outer component makes
the positive core one global regular cell. The endpoint period calculation
then proves a strictly positive finite `omega_min` for each such fixed
carrier. This is profile/cell uniqueness, not an inverse-function theorem for
a three-dimensional nonaxisymmetric branch.

**Verdict:** established for the fixed-domain thin-family profile and its
global positive regular cell. No nonlinear branch uniqueness follows.

## Additional exact atoms and route boundaries

The explicit interior fixed-`|n|>=2` divergence-free displacement has the
correct cylindrical component and KKS calculation

    Omega_KKS(q_test,conjugate(q_test))
      =4*pi*i*rho_0 integral r^2 zeta f g dr dz !=0.

It is a physical DA symplectic seed, not an eigenmode. Translation seeds are
isotropic zero modes, tilted rotations are generalized zero chains ending in
translations, and axial rotation stabilizes the ring; this refutes only the
symmetry-seed route to a nonzero-frequency internal state. The toroidal
selection rule `n=(k-2j)n_0` is exact, while a finite cross-sectional `m`
lattice is not. The Hamiltonian/KKS Green identity is established on smooth
compact DA tangents, with extension to the completed weighted/interface graph
still dependent on Unit 5.

No target-local API, verifier, test, captured oracle output, or production
numerical evidence exists in the frozen `0066` inventory. The evidence is
analytic. The passing schema and earlier supplier tests cannot replace the
continuum regularity or operator proofs, and no unchanged oracle needed a
rerun for this correction.

## Final verdict and next dependency

**Established:** for every sufficiently thin fixed `p>=6` Cao carrier in the
declared `H^3` vorticity scope, `0066` establishes the nonzero boundary jet,
one-center regular foliation and finite positive frequency bounds; the
ordinary unweighted collar lower-bound failure; the weighted tangent and
smooth-path Hodge identities; a compact physical nonzero-KKS DA seed;
symmetry and toroidal-selection results; and, most importantly, the
fixed-nonzero regular-cell constrained DA transport-band inclusion in the
upper-semi-Fredholm and standard non-Fredholm essential spectra.

**Still conditional:** the converse constrained direct-integral and
boundary/interface spectrum, hence the claimed open essential gap; the
nonlinear weighted quotient chart; the full-vector high-`n` Piola/Hodge
`D->X` estimate; isolation of the limiting `m=0` modes from the whole column;
Riesz transfer to Cao; positive Cao Krein modes; the two-mode doublet and its
noncommuting deformations; and equations (38)--(39). No nonlinear branch,
stability neighborhood, P2/LP2, particle, quantum, charge, or relativistic
conclusion is licensed.

The next construction is the converse constrained DA direct-integral plus
boundary/interface spectral theorem. In parallel after that domain is fixed,
the massive route needs the complete vector Piola/Hodge ledger at
`n_delta~delta^(-1)` and whole-column common-domain isolation. Only those
results can activate a contour transfer or nonlinear two-mode range equation.

The single bounded correction passes and fully closes the findings of this
review. The route verdict is **established at the exact regular-cell
essential-inclusion scope**, with the stronger gap/doublet/transfer route
honestly open rather than refuted.
