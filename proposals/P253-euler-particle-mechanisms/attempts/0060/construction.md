# Exact neutral-cell moment cancellation and band boundary

## 1. Exact stationary composite, and why its pressure cannot couple cells

Let `(u_*,p_*)` be one accepted smooth stationary compact-velocity Euler
component at physical density `rho_0`. Write
`pi_*=(p_*-p_infinity)/rho_0` for kinematic pressure normalized to vanish on
the connected exterior. Put six rotated copies at well-separated signed
coordinate-axis centers. Their velocity supports, including flat collars,
are disjoint. On the support of one component every other velocity vanishes
and every other `pi` has zero gradient; outside all supports every velocity
and pressure gradient vanishes. Hence

    u_cell=sum_a u_a,
    pi_cell=sum_a pi_a,
    p_cell=p_infinity+rho_0 pi_cell                        (1)

is an exact stationary Euler field with positive additive energy and tags.
This convention avoids summing the ambient pressure constant six times.

There is a stronger identity than the preregistered generic quadrupole
cancellation. Outside the compact velocity support, Euler gives
`grad p_*=0`. Hence the Newton potential of
`partial_i partial_j(u_i u_j)` is constant there and every exterior pressure
multipole vanishes. At the first level this follows from the virial identity

    integral u_i u_k dx=-delta_ik integral pi dx
      =-(delta_ik/rho_0) integral (p-p_infinity) dx,      (2)

obtained by multiplying stationary momentum balance by `x_k` and integrating.
Thus `M_ik` is isotropic and its contraction with the harmonic Newton Hessian
is zero. Higher stationary identities cancel the higher moments as well.
Translating and rotating a field with identically constant exterior pressure
cannot restore a tail.

**Static-pressure route verdict:** refuted as an intercell coupling mechanism.
The exact compact assembly is positive carrier geometry, but its stationary
pressure field supplies no intercell force rather than an `r^-5` interaction.

## 2. Six sites, the `C4` stabilizer, and exact translated moments

Let `O` be the 24-element orientation-preserving octahedral group and let
`H=C4` stabilize the base signed-axis site and its axisymmetric carrier. The
stationary cell has six sites indexed by `O/H`, not 24 disjoint sites. For a
compact-core DA seed `eta_*=-L_xi Omega_*`, define its stabilizer projection

    eta_bar=(1/4) sum_(h in H) h_* eta_*.                 (3)

A generic nonzero seed does not imply `eta_bar!=0`: if
`eta_*=(I-h_*)zeta` for a generator `h` of `H`, the sum telescopes to zero.
If a **nonzero** `eta_bar` is exhibited, then its six coset images have
disjoint supports and give a nonzero product-leaf tangent. Equivalently their
sum is one quarter of the full 24-element average of `eta_*`. A second valid
construction is a genuinely free 24-site orbit with trivial stabilizer and
24 disjoint supports. Neither construction uses signed algebraic weights.

The moment statement includes translations. For local seed moments

    z_a=integral eta_a dy,
    q_ba=integral y_b eta_a dy,
    s_bca=integral y_b y_c eta_a dy,

the copy at `X_g` with rotation matrix `R_g` has

    Z_i^g=R_ia z_a,
    Q_ji^g=X_gj Z_i^g+R_jb R_ia q_ba,
    S_jki^g=X_gj X_gk Z_i^g
       +X_gj R_kb R_ia q_ba+X_gk R_jb R_ia q_ba
       +R_jb R_kc R_ia s_bca.                            (4)

For compact curl vorticity `z=0`; integration by parts makes `q_ba`
antisymmetric, while `s_bca=s_cba` in its two spatial indices. With
`X_g=R_g X_*`, the full group average in (4) acts on the exact translated
tensors. Its averaged rank-two tensor is isotropic and antisymmetric, hence
zero. Its averaged rank-three tensor is symmetric in `j,k`; although proper
octahedral symmetry admits the pseudoscalar `epsilon_jki`, that symmetry
excludes it, so the average is zero. For a nonzero `H`-invariant seed the
six-coset sum is one quarter of this full average and has the same
cancellations.

The exact group projectors used in that last step are

    sum_R R a=0,
    sum_R R M R^T=8 tr(M) I,
    sum_R R_(i a)R_(j b)R_(k c)T_(abc)=0                (5)

for every vector `a`, rank-two `M`, and tensor `T` symmetric in its first two
indices. Thus, conditional on one of the nonzero tangent constructions above,
the Hodge dipole and quadrupole cancel. A rank-four cubic invariant is allowed
and can supply the octupole, so the robust upper bound is

    |v_O(x)|=O(|x|^-5),       |grad v_O(x)|=O(|x|^-6).    (6)

Equation (6) is an upper bound; an additional seed moment cancellation may
make it faster. It holds for the trivial octahedral tangent representation.
Generic symmetry-breaking tangents can restore the `r^-3` dipole and remain
in the physical complement.

For the direct-product statement, let `K_g` be each vorticity core and choose
`xi_g` divergence free, compactly supported in its interior, and zero on a
flat collar. Then

    eta_g=-L_(xi_g) Omega_g=curl(xi_g cross Omega_g),
    eta_cell=sum_g eta_g.                                (7)

Each component separately retains its material-circulation rows
`delta Gamma_(g,a)=0`; any solid-torus harmonic coefficients are retained and
fixed rather than discarded; and a centered slice imposes the three
component translation rows
`Omega_KKS,g(eta_g,-L_(e_j)Omega_g)=0`. The stabilizer row is
`h_*eta_g=eta_g` in the six-coset construction. These rows, support/collar
conditions, and the displacement `xi_g` are imposed component by component
before taking the direct sum. Rigid rotations commute with curl and Lie
derivative. The KKS form is a direct sum because its integrand is supported in
the disjoint cores, but the kinetic Hessian, Hodge velocity, and pressure have
global cross terms and are not asserted to split. Thus (7) is an exact
product-leaf tangent, not an invariant product dynamics. Its perturbation
energy is finite because the conditional `r^-5` tail is square integrable.

**Moment route verdict:** the translated octahedral tensor projector and its
conditional `r^-5`/`r^-6` implication are established. The six-site nonzero
DA tangent is blocked by the missing nonzero `C4`-invariant projected seed; it
remains the neutral-cell candidate, not yet an invariant Euler normal mode.

## 3. Lattice summability and the conditional band theorem

The pointwise statement is exact: an `r^-5` field has finite scalar zeroth and
first lattice moments in three dimensions. A Bloch theorem requires more. If
normalized invariant cell modes have a full pressure/Leray coupling satisfying

    ||J(R)||_(D(A_0)->X)+||partial_R J(R)||_(D(A_0)->X)
       <=C (1+|R|)^-5,                                  (8)

uniformly in cell orientation, where `D(A_0)` is the physical graph domain,
then

    sum_R ||J(R)||<infinity,
    sum_R |R| ||J(R)||<infinity.                          (9)

The comparison integrands are `r^-3` and `r^-2`. Therefore

    A(k)=A_0+sum_(R!=0) exp(i k dot R)J(R)               (10)

is continuously differentiable in `k` as an operator `D(A_0)->X`, with
`partial_k A=sum iR exp(i k dot R)J(R)`. For any simple isolated positive KKS
band, a uniform graph-resolvent estimate then gives a differentiable Riesz
family and bounded group velocity. Pointwise decay of one Hodge velocity does
not prove (8): the normalized mode, pressure/Leray derivative, interface and
operator-domain estimates are all additional inputs.

Those inputs are not currently available. C-CST-018 gives prepared
finite-window response rather than an autonomous invariant cell mode. The
static compact translations have no coupling. A generic DA perturbation
leaks from the trivial octahedral representation once a Bloch phase breaks
the within-cell symmetry unless the full projected/complement operator is
controlled.

**Bloch route verdict:** blocked by the named invariant-mode, graph-domain
`r^-5` coupling, and complement construction. The conditional moment theorem
makes the scalar lattice threshold plausible and checkable; it does not create
a `C^1` Bloch operator or spectral band.

## 4. Isolated scattering alternative

Two separated symmetric composites have a finite first-moment linear Hodge
coupling in the protected tangent sector. But no reviewed result controls the
nonlinear persistence of that sector during encounter, and symmetry-breaking
deformation can restore lower multipoles. The isolated-scattering route is
therefore blocked by a same-family persistence/modulation theorem with the
dangerous moment coordinates retained.

## 5. Strongest result and next construction

The useful exact result is the translated octahedral moment projector: for a
nonzero trivial-sector compact DA tangent, it removes the Hodge dipole and
quadrupole using only physical rotated copies, giving an `r^-5` velocity upper
bound and the scalar finite first spatial moment needed by a strong
effective-band estimate. The attempted static pressure interaction is exactly
absent, a consequence of compact stationarity rather than a numerical
cancellation. Constructing the required nonzero `C4`-invariant tangent on the
six-site cell remains open rather than being inferred from a generic seed.

The next construction must first exhibit a nonzero `C4`-invariant compact DA
tangent (or use the free 24-site alternative), then choose an actual positive
internal mode and prove the full periodic Euler/Leray operator-domain estimate
(8), an isolated invariant trivial-sector Riesz family, and uniform leakage
control in `k`. Its Bloch derivative then follows from (9)--(10). This remains
an effective band with algebraic full-field leakage, not a Lorentz cone,
action selection, detector, exchange rule, electron, or neutrino.
