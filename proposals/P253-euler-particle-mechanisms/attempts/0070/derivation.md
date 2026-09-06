# Can the electromagnetic U(1) sector emerge from Euler?

## 1. Clebsch gauge is a redundancy, not a field sector

On a regular local chart write the Euler velocity one-form as

    u_flat=d phi+alpha d beta.                            (1)

For any smooth one-variable function `F`, the transformation

    phi'=phi-F(beta),   alpha'=alpha+F'(beta),   beta'=beta (2)

leaves (1) unchanged.  It also leaves the Clebsch kinetic one-form
`phi_t+alpha beta_t` unchanged, because the two added terms cancel.  Its
infinitesimal generator therefore satisfies

    delta u= -F'(beta) grad beta+F'(beta) grad beta=0.    (3)

Consequently its pulled-back kinetic energy, physical Euler KKS pairing, and
every observable depending only on the Euler state vanish along this
direction.  This is an exact local gauge redundancy of coordinates.  It is
not an independent field with curvature energy or a charged source.

The physical two-form `d u_flat=d alpha wedge d beta` is vorticity.  Euler
advects it; it does not satisfy the independent Maxwell displacement-current
equation.  Nonzero helicity additionally obstructs one global regular
Clebsch pair, as independently reviewed in 0011.  Passing to multiple charts
repairs representation coverage but does not turn a kernel direction into a
new physical canonical pair.

**Route A verdict:** the Clebsch/relabeling transformation is established as
an exact gauge redundancy and refuted as an exact derivation of the 0068
Maxwell sector.  This verdict does not cover a new collective connection with
an independently derived physical curvature and action.

## 2. The physical principal symbols do not match

Linearize constant-density incompressible Euler about rest.  For a Fourier
mode with `k!=0`,

    partial_t v=-i k pi,       k dot v=0.                 (4)

Applying the transverse projector `P_T(k)` gives `partial_t v=0` on the two
physical velocity polarizations.  A material tag about a constant reference
also obeys `partial_t delta chi=0`.  About a uniform relative background `U`,

    (partial_t+U dot grad)v=-grad pi,
    (partial_t+U dot grad)delta chi=0,                    (5)

so every physical root is the convective frequency

    omega=U dot k.                                       (6)

By contrast, the source-free 0068 gauge field satisfies

    E_t=c_EM^2 curl B,        B_t=-curl E,                (7)

with `k dot E=k dot B=0`.  Its four-dimensional transverse state has two
polarizations on each temporal branch

    omega=+c_EM |k|,    omega=-c_EM |k|.                 (8)

Writing `lambda` for the time-generator spectral parameter, the physical
characteristic factors are

    chi_Euler=(lambda+i U dot k)^2,
    chi_Maxwell=(lambda^2+c_EM^2 |k|^2)^2.               (9)

At rest the Euler root is zero; at uniform flow it is a repeated one-sided
convective root.  It is neither the isotropic two-sided cone (8) nor a Gauss
constraint on an independent electric field.

Let `S(k)` be any time-local, translation-invariant finite-order differential
or elliptic pseudodifferential change of physical variables whose principal
symbol is invertible on the declared quotient.  Its transformed generator is
`S A S^(-1)` at principal order, so its characteristic polynomial and roots
are unchanged.  Equations (6)--(9) therefore refute an exact local invertible
Euler-to-Maxwell conjugacy at rest or uniform background.

The same high-frequency conclusion holds locally about a smooth background:
the principal Euler characteristic is the repeated material root
`omega=u_0(x) dot k`; strain and Hodge projection enter the amplitude/lower
order system.  This statement does not classify low-frequency Bloch bands,
homogenized structured media, time-nonlocal maps, singular limits, or an
enlargement by new degrees of freedom.

**Route B verdict:** established as a characteristic-set obstruction for the
frozen class of local invertible physical field redefinitions.  It is not a
no-go for structured-background collective limits or foundation extensions.

## 3. Structured backgrounds separate propagation from charge

The accepted `C-CST-018` response contains a prepared massless transverse
sector.  It already evades the rest-state statement by using a structured
Euler background and a finite-window collective response.  It does not evade
the charge representation test.  For a local homogeneous rotation scalar
`q`, covariance gives at nonzero `k`

    f_i(k)=i k_i F(|k|^2)q(k),                            (10)

and the transverse response obeys

    P^T_ij(k)f_j(k)=0.                                   (11)

An internal vector produces a nonzero source but the exact 0067 response is
orientation dependent and has the attractive positive-elastic-field sign.
A triplet removes the tensor only after adding two sectors and a common-frame
lock, while preserving the attractive sign.  Thus the current accepted
collective sector supplies neither the independent Gauss variable nor the
electric scalar source.

A successful structured-background route must construct, on one autonomous
state space, both a longitudinal constrained Green response and a transverse
propagating pair, and it must derive the reciprocal source and carrier force
from one action.  No accepted or author-stage construction presently supplies
those four rows.  This is a named missing construction, not evidence that no
nonlinear or nonlocal Euler collective limit can do so.

**Route C verdict:** refuted for direct reuse of the accepted prepared
transverse sector; blocked for a new structured-background homogenization by
the absent autonomous longitudinal/Gauss sector and same-carrier source law.

## 4. What the 0068 extension adds

Attempt 0068 adds precisely the structures absent above:

    rho_q=g chi,       J=g chi u,
    div(epsilon E)=rho_q,
    (1/mu)curl B-epsilon E_t=J.                          (12)

The current follows from an existing transported Euler tag, but `(E,B)`, its
positive field energy, and `epsilon,mu` are new.  The coupling `g` is also a
new normalization.  A common material-map plus gauge action then supplies the
Lorentz force and electric-sign Coulomb energy.  The construction is a
coherent minimal extension because deleting the independent gauge state
returns the characteristic and projector obstructions above.

Multiplying the complete classical action by a constant changes its
symplectic/action scale without changing the equations.  Independently,
`epsilon*mu` fixes the gauge speed.  Even the extension therefore selects
neither `hbar` nor a charge quantum.  It supplies a classical charge sector,
not an electron.

**Route D verdict:** established as a minimal-delta ledger for the routes
covered here.  The exact added objects are one constrained U(1) field, its
positive action, and three independent constants.  Adoption remains a user
foundation decision; bare-Euler carrier and doublet routes continue in 0066.

## 5. Continuation toward the particle objective

The campaign now has two sharply separated positive lines:

1. **Bare Euler:** construct the 0066 persistent carrier and positive physical
   doublet.  Charge still requires a nonlinear/configuration-dependent
   topological source or a new longitudinal collective sector.
2. **Minimal charged-fluid extension:** couple the reviewed transported tag to
   the 0068 U(1) field on the same persistent carrier.  Then derive rather
   than assign the carrier charge integral, magnetic moment, and action scale.

Neither line yet supplies spin half, fermionic exchange, Born probabilities,
a Lorentz one-particle representation, or the neutrino chiral/flavour sector.
Those remain active parent obligations rather than defects in the exact
route-scoped results above.
