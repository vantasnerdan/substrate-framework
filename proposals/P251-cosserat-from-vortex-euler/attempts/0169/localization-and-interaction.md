# What exact pressure localization can transfer to the current construction

## 1. The full Euler identity

Use pressure divided by density and a smooth stationary Euler pair
(u,p), with div u=0 and Du u=-grad p. For a smooth scalar cutoff chi,
set U=chi(p)u and choose P by P'(p)=chi(p)^2. Direct chain rules give

    div U=chi'(p) A p,
    DU U+grad P=chi(p)chi'(p)(A p)u,  A=u.grad.

Thus the advected-pressure condition is the exact license for arbitrary
pressure modulation; it is not a consequence of stationarity or of a
globally defined Bernoulli function. Where grad p is nonzero and
chi' is active, the first residual alone exposes a false localization.

For the actual triangular field of0161/0167, including ALL generalized
Bernoulli lifts W=sqrt(C+lambda²psi²), the physical pressure is unchanged
and independent of the axial coordinate. Hence

    u.grad p=v.grad p=-(1/2)v.grad |v|².

The lambda²psi² contribution differentiates to zero because v.grad psi=0.
The remaining expression is nonzero for the actual three-wave field, as
the exact point evaluation in verify.py establishes. It is analytic and
does not depend on C. Raising the axial speed therefore cannot create
the advected-pressure license for this geometry.

More generally, on this connected analytic cell any nonconstant smooth
chi on the attained pressure interval has chi' nonzero on an open
subinterval. Its preimage is nonempty and open. If chi'(p)A p vanished
there, analyticity would force A p to vanish on the entire cell,
contradicting the exact evaluation. Thus the pressure-cutoff class does
not localize this field at any choice of C. This is a coverage argument
for that scalar modulation route, not a no-go for steady-Euler gluing.

## 2. Repair by changing the actual stationary field

The localizable axisymmetric candidates impose an additional equation
alongside stationary Grad-Shafranov balance. In the source orientation,

    u=(psi_z e_r-psi_r e_z+F(psi)e_phi)/r,
    -Delta_*psi=(F²/2)' +r² P_plasma'(psi),
    |grad psi|²+F²=2r² A_speed(psi).

The hydrodynamic pressure is -P_plasma-A_speed. This second constraint
changes the actual core jets; it cannot be obtained by calling the
original Bernoulli lift generalized force-free. The existing global
solutions leave room for investigating an annular optical parcel in
their smooth shell, but do not retain the triangular elliptic-core or
whole-cell acoustic estimate. Those are newly required constructions.

The generic local torus of0136 remains a separate unlocalized candidate.
Its Grad-Shafranov solution supplies exact Euler and a closed tube but
does not satisfy the extra speed condition by default. A genuine global
elliptic extension or a newly constructed compatible localizable profile
is needed for that transfer. The constant-curl EPS Runge theorem is not
changed or used beyond its actual equation.

## 3. Exact disjoint assembly and its actual static interaction

Suppose a compact smooth Euler template and its pressure are given in
a bounded support region, with pressure constant on the exterior.
Choose separated copies by translations and orthogonal transformations,
and subtract each exterior pressure constant. Their sum is an exact
smooth Euler field wherever the supports remain disjoint: every cross
velocity product and cross derivative vanishes pointwise. A uniformly
separated periodic placement gives a bounded smooth periodic field.
Whole-field uniform cell phase and Haar rotation yield a stationary,
isotropic finite-energy-density law. This is an actual nonlinear
superposition license from disjoint supports, not linear Beltrami Runge.

Its complete kinetic energy is a sum of the individual energies.
Independent rigid rotations and translations that keep supports disjoint
leave every summand invariant. Therefore the mixed static Hessian for
these whole-template positions/angles is exactly zero. The Biot-Savart
cross energy agrees: each template's actual velocity is zero on the
others, so retaining the complete inverse-curl tails does not restore
an interaction that has already cancelled in this construction.

This refutes a positive static whole-template locking coefficient from
the disjoint family. It does not refute an internal shape-polarization
oscillator, an induced-flow dynamical current, or overlapping velocity
fields. Any proposed coupling from those mechanisms must use their own
actual variations rather than the independently rotated-template ones.

## 4. Next executable route

The scalar cutoff of the existing triangular field is refuted by its
nonadvected pressure. The disjoint global assembly is established, while
its whole-template static locking route is refuted by exact additive
energy. Neither verdict propagates to the global geometry obligation.

The continuation candidates are actual localizable-shell internal
polarization, smooth vorticity-supported rings with noncompact velocity
interaction, and the unchanged EPS constant-curl whole-field construction
with generic-direction rather than axial-only acoustic response. The
source theorem for a discontinuous ring is not enough for the smooth
route; explicit regularity and profile continuation remain to be derived.
0163/0164/0166 meanwhile continue their live response constructions.
