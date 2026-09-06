# Exact neutral-cell moment cancellation and band boundary

## 1. Exact stationary composite, and why its pressure cannot couple cells

Let `(u_*,p_*)` be one accepted smooth stationary compact-velocity Euler
component. Put six rotated copies at well-separated signed coordinate-axis
centers. Because their velocity supports are disjoint, every cross-advection
term vanishes. Each other pressure is constant on a given support, so

    u_cell=sum_a u_a,       p_cell=sum_a p_a              (1)

is an exact stationary Euler field with positive additive energy and tags.

There is a stronger identity than the preregistered generic quadrupole
cancellation. Outside the compact velocity support, Euler gives
`grad p_*=0`. Hence the Newton potential of
`partial_i partial_j(u_i u_j)` is constant there and every exterior pressure
multipole vanishes. At the first level this follows from the virial identity

    integral u_i u_k dx=-delta_ik integral p dx,          (2)

obtained by multiplying stationary momentum balance by `x_k` and integrating.
Thus `M_ik` is isotropic and its contraction with the harmonic Newton Hessian
is zero. Higher stationary identities cancel the higher moments as well.
Translating and rotating a field with identically constant exterior pressure
cannot restore a tail.

**Static-pressure route verdict:** refuted as an intercell coupling mechanism.
The exact compact assembly is positive carrier geometry, but its stationary
pressure field supplies a flat band rather than an `r^-5` interaction.

## 2. Physical octahedral projection of a DA tangent

Let `eta_*=-L_xi Omega_*` be any nonzero compact-core dynamically accessible
vorticity tangent. Its Hodge velocity `v_*=curl(-Delta)^-1 eta_*` need not be
compact. Place and rotate the tangent with the same six-component cell, or
equivalently average it over the 24-element proper octahedral group `O`:

    eta_O=sum_(R in O) R_* eta_*.                         (3)

Because the component supports are disjoint, (3) is a nonzero tangent to the
product coadjoint leaf whenever the seed is nonzero. It contains no signed
algebraic copy. The exact group projectors give

    sum_R R a=0,
    sum_R R M R^T=8 tr(M) I,
    sum_R R_(i a)R_(j b)R_(k c)T_(abc)=0                (4)

for every vector `a`, rank-two `M`, and tensor `T` symmetric in its first two
indices.  The compact divergence-free vorticity first moment is
antisymmetric and trace free, so its rank-two orbit sum vanishes.  The next
vorticity moment is symmetric in its two spatial indices and has the
rank-three form in (4), so it also vanishes.  Thus the Hodge dipole and
quadrupole cancel. A rank-four cubic invariant is allowed and can supply the
octupole, so the robust bound is

    |v_O(x)|=O(|x|^-5),       |grad v_O(x)|=O(|x|^-6).    (5)

Equation (5) is an upper bound; an additional seed moment cancellation may
make it faster. It holds for the trivial octahedral tangent representation.
Generic symmetry-breaking tangents can restore the `r^-3` dipole and remain
in the physical complement.

The DA property is exact: rigid rotations commute with curl, Lie derivative,
and the Euclidean Hodge operator, while the sum over disjoint components is
the direct product leaf tangent. The perturbation energy is finite because
`r^-5` is square integrable.

**Moment route verdict:** established for the cellwise octahedral DA tangent.
It is a physical cancellation theorem, not yet an invariant Euler normal
mode.

## 3. Lattice summability and the conditional band theorem

If a normalized invariant cell mode has the tail (5), its leading linearized
coupling to a distant compact core is `J(R)=O(|R|^-5)`. On a three-dimensional
lattice,

    sum_R ||J(R)||<infinity,
    sum_R |R| ||J(R)||<infinity,                          (6)

because the comparison integrands are `r^-3` and `r^-2`. Therefore

    A(k)=A_0+sum_(R!=0) exp(i k dot R)J(R)                (7)

is continuously differentiable in `k`, with
`partial_k A=sum iR exp(i k dot R)J(R)`. For any simple isolated positive
KKS band, standard Riesz differentiation then gives a bounded group velocity.
This conclusion is exact conditional operator algebra once the invariant
mode, isolation, and full Euler coupling bound are supplied.

Those inputs are not currently available. C-CST-018 gives prepared
finite-window response rather than an autonomous invariant cell mode. The
static compact translations have no coupling. A generic DA perturbation
leaks from the trivial octahedral representation once a Bloch phase breaks
the within-cell symmetry unless the full projected/complement operator is
controlled.

**Bloch route verdict:** blocked by the named invariant-mode and complement
construction. The moment theorem makes the required lattice kernel
summability plausible and checkable; it does not create a spectral band.

## 4. Isolated scattering alternative

Two separated symmetric composites have a finite first-moment linear Hodge
coupling in the protected tangent sector. But no reviewed result controls the
nonlinear persistence of that sector during encounter, and symmetry-breaking
deformation can restore lower multipoles. The isolated-scattering route is
therefore blocked by a same-family persistence/modulation theorem with the
dangerous moment coordinates retained.

## 5. Strongest result and next construction

The useful exact result is the octahedral dynamically accessible moment
projector: it removes the Hodge dipole and quadrupole using only physical
rotated copies, giving an `r^-5` velocity bound and the finite first
spatial moment needed by a strong effective-band estimate. The attempted
static pressure interaction is exactly absent, a consequence of compact
stationarity rather than a numerical cancellation.

The next construction must choose an actual positive internal mode of one
compact component, form its physical octahedral orbit, and prove that the
full periodic Euler/Leray operator has an isolated invariant trivial-sector
Riesz family with leakage controlled uniformly in `k`. Its Bloch derivative
then follows from (6)--(7). This remains an effective band with algebraic
full-field leakage, not a Lorentz cone, action selection, detector, exchange
rule, electron, or neutrino.
