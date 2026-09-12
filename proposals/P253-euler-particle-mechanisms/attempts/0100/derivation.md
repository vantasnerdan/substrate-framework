# Ertel axial current and the same-carrier lock

## 1. Exact local material invariant

Let `u` be a classical incompressible Euler velocity, let
`omega=curl u`, and let `chi` be a transported true scalar.  With
`D_t=partial_t+u dot grad`,

    D_t omega=(omega dot grad)u,
    D_t chi=0,
    D_t grad chi=-(grad u)^T grad chi.                    (1)

Define

    q_E=omega dot grad chi.                               (2)

In indices, the two stretching terms cancel after relabeling `i,j`:

    D_t q_E
      =omega_j partial_j u_i partial_i chi
       -omega_i partial_i u_j partial_j chi
      =0.                                                 (3)

Since `div u=0`, this is the exact continuity law

    partial_t q_E+div(q_E u)=0.                           (4)

For a true scalar `chi`, `grad chi` is polar and `omega` is axial.  Thus
`q_E` is a pseudoscalar and `q_E u` is an axial spatial current.  If `chi`
is instead a pseudoscalar, `grad chi` is axial, so the density is scalar and
the spatial current polar.  Equation (4) is a Galilean continuity-current
tuple; it is not a Lorentz chiral four-current.

## 2. Exact forced Ertel current

For constant mass density `rho_m`, add a body acceleration `f`:

    D_t u=-grad p/rho_m+f.                                (F1)

Taking curl changes the first equation in (1) to

    D_t omega=(omega dot grad)u+curl f.                   (F2)

The stretching cancellation in (3) then leaves

    D_t q_E=(curl f) dot grad chi
           =div(f cross grad chi).                        (F3)

The exact forced continuity law is therefore

    partial_t q_E+div(q_E u-f cross grad chi)=0.           (F4)

The sign follows from
`div(f cross grad chi)=grad chi dot curl f`.

In the reviewed Euler--Maxwell extension, write the electromagnetic material
charge tag as `chi_c(P)` to distinguish it from later Ertel labels:

    f=(g chi_c/rho_m)(E+u cross B).                       (F5)

On the steady charged Cao branch,

    E+u cross B=-grad Phi-H grad P,
    grad chi_c=chi_c'(P)grad P.                           (F6)

Although the regular axisymmetric tag has `q_E=0`, its corrected axial
spatial flux need not vanish:

    j_F=-f cross grad chi_c
       =(g/rho_m)chi_c chi_c'(P)grad Phi cross grad P
       =(g/(2rho_m))grad Phi cross grad(chi_c^2)
       =curl[(g Phi/(2rho_m))grad(chi_c^2)].               (F7)

Thus `div j_F=0` pointwise for constant `g,rho_m`. It is a toroidal
magnetization/superpotential current made from two poloidal gradients. Its
density and global Ertel charge remain zero, and its integrated transport
content is a boundary row. This closes one explicit superpotential member of
the wider first-derivative current class left open by 0097, without exhausting
that class or creating a Lorentz weak current.

## 3. The global charge is an exact boundary row

The vorticity identity `div omega=0` gives

    q_E=div(chi omega),
    integral_D q_E dx=integral_(boundary D) chi omega dot n dS.  (5)

Consequently `integral_R3 q_E dx=0` for smooth compactly supported
`chi omega`, and also under decay that kills the sphere flux.  The same
conclusion holds on a closed compact boundaryless manifold, or a noncompact
boundaryless domain with adequate infinity decay, when `chi omega` is
globally smooth and single-valued: topology alone does not invalidate Stokes.

A nonzero global Ertel charge must therefore enter through an actual boundary
flux, a puncture or distributional defect, or a non-global/multivalued label
for which `chi omega` is not a globally defined smooth vector field.  The
local conservation law (4) remains useful, but it supplies no independent
nonzero whole-space charge in the regular single-valued class.

## 4. Constitutive vector/axial lock and monodromy

First test the local relation in unforced Euler:

    q_E=lambda chi.                                       (6)

Because both `q_E` and `chi` are transported in Route A0, an initially imposed
relation is preserved precisely when

    (D_t lambda)chi=0.                                    (7)

Thus an unforced variable `lambda` must itself be advected wherever `chi` is nonzero.
For true scalar `chi`, covariance types `lambda` as a pseudoscalar if (6) is
to relate the pseudoscalar `q_E` to `chi` without hiding a parity-breaking
background.

Along a vorticity line `dx/dtau=omega(x)`, (6) reads

    d chi/dtau=lambda chi,
    chi(tau)=chi(0) exp(integral_0^tau lambda(s) ds).      (8)

For a closed line of period `T` carrying nonzero single-valued `chi`,

    exp(integral_0^T lambda d tau)=1.                     (9)

For real `lambda`, (9) is `integral_0^T lambda d tau=0`.
A real nonzero constant lock is therefore impossible on such a line, while a
variable lock with zero line integral is not excluded. The monodromy statement
is kinematic at each fixed time and therefore also constrains forced locks.

Under forcing, `q_E` is no longer materially transported. Substituting (6)
into (F3) instead gives the distinct exact condition

    chi D_t lambda=(curl f) dot grad chi
                  =div(f cross grad chi).                (9a)

Consequently forcing can in principle select or sustain a variable lock. The
unforced `D_t lambda=0` condition must not be imposed on that route. This
attempt establishes (9a), but it does not construct a force-selected lock
that also satisfies the closed-line monodromy, parity, action, and endpoint
interaction rows.

## 5. Exact test on the charged Cao carrier

In the Cao convention used by the reviewed carrier packages,

    omega=r zeta(r,z) e_theta=zeta(r,z) partial_theta.    (10)

For the reviewed transported material tag `chi=F(I)` on the regular tagged
band, axisymmetry gives `partial_theta chi=0`.  Therefore

    q_E=omega dot grad chi=zeta partial_theta chi=0.       (11)

No center, edge, or cylindrical metric factor changes (11); the tag is
supported inside the regular band and both objects are smooth there.  The
Ertel route consequently gives no axial density on this same axisymmetric
no-swirl charged carrier.  This is a geometry-specific refutation, not a
no-go for nonaxisymmetric carriers, multivalued labels, or defects.

## 6. Exact defect-domain Cao continuation

The refutation (11) uses a globally real axisymmetric tag. The failure points
to a distinct object already present in the geometry. On
`M=R^3` minus the symmetry axis, the azimuthal phase is a circle-valued
geometric label. Its closed one-form

    d theta=(-y dx+x dy)/(x^2+y^2)                       (12)

is smooth and global on `M` but is not exact there. At the translating
axisymmetric no-swirl Cao equilibrium, every local lift obeys

    D_t theta=u dot grad theta=0.                         (13)

The local forced Ertel calculation therefore patches to the global density

    q_theta=d theta(omega)=d theta(zeta partial_theta)=zeta. (14)

The vorticity support is a compact toroidal core separated from the removed
axis, so (14) is smooth and integrable even though the phase has a defect.
With the exact Cao circulation convention

    kappa=integral zeta(r,z) r dr dz,

one obtains the exact equilibrium integral

    Q_theta=integral_M q_theta dx=2 pi kappa.             (15)

This does not contradict (5): there is no global real `chi` for which
`chi omega` supplies the exact-divergence primitive. The removed axis and
non-exact `d theta` are the precise topological data used by the construction.

For this forced equilibrium, the conserved local flux is
`q_theta u-f cross dtheta`. Both `zeta` and the electromagnetic force are
supported away from the removed axis and infinity, so the inner-axis and
outer boundary fluxes vanish at the base. In an axisymmetric no-swirl
evolution the fixed geometric `theta` remains material and the same boundary
argument conserves (15), provided those support separations persist.

For a general nonaxisymmetric or swirling perturbation, however,

    D_t theta=u_theta/r.                                  (16)

so the fixed geometric azimuth is not a transported material label. A
dynamical extension must instead include an independently advected
circle-valued phase `Theta`, initialized as `theta`, whose defect line and
domain move with the flow; its total is conserved only after the corresponding
moving inner-boundary and outer-flux rows are proved. Thus (15) is
unconditionally an equilibrium integral, not an all-perturbation invariant.

Equations (12)--(15) establish a nonzero same-carrier Ertel density at
defect-domain equilibrium scope. It is not the electromagnetic charge density
`g chi_c(P)`. Its magnitude remains the continuous Euler circulation. The
chosen axis/circle phase does not transform as the true scalar assumed in
Section 1 under unrestricted `O(3)`, supplies no Lorentz chiral projector, and
has no derived symplectic/action coefficient or endpoint interaction. It is
therefore a useful P6 supplier rather than a neutrino or charge/action
selection mechanism.

## 7. Failure-derived positive continuations

The strongest next candidates are materially distinct:

1. construct a persistent nonaxisymmetric carrier with a smooth transported
   true scalar and nonzero local `omega dot grad chi`, while retaining an
   electron-linked interaction current and the shared P4 action;
2. promote the exact azimuthal defect-domain current (12)--(15) by deriving
   an advected circle phase and moving defect domain, finite-energy material
   action, parity/rotation representation, endpoint interaction, and relation
   to the electron current;
3. replace the prohibited real exponential lock on closed lines by a compact
   material connection with derived holonomy and endpoint interaction frame;
4. classify multilabel first-derivative currents together with the still-open
   0097 invariant-tensor/superpotential sector.

Each candidate must live on an independently persistent P2 carrier and must
join the shared P4 spin/action representation before it can contribute to P6.
The exact local axial current (4), by itself, supplies neither neutrino flavor
oscillation nor a Lorentz chiral interaction.
