# Physical frame response and the invariant first-gradient coupling

## 1. Direct material variation, not an assumed polar frame

For a finite positive material tag with central second mass tensor I,
an affine material displacement eta=a+E r changes the tensor by

    delta I=E I+I E^T.                                        (1)

This follows by differentiating both factors in the defining material
integral. The central translation cancels exactly. For a planar tag whose
principal axes initially diagonalize I=diag(A,B), A!=B, the actual
central-quadrupole angle is half arg(Ixx-Iyy+2i Ixy). Thus

    delta theta=(A E_yx+B E_xy)/(A-B)
               =omega+(A+B)/(A-B) e_xy,
    omega=(E_yx-E_xy)/2, e_xy=(E_yx+E_xy)/2.                  (2)

The coefficient of an actual rigid rotation is exactly one. Symmetric
strain generally also rotates the measured angle. Averaging the complete
tag and affine preparation over the in-plane reference angle cancels
that linear symmetric-strain row, while preserving omega. This is a
kinematic statement at the same instant, not a theorem that actual Euler
cell displacements are affine throughout an acoustic interval.

For a non-affine actual displacement, the additional row is exactly

    delta theta_extra=Im[2 integral rho w (x+iy)
                              (eta_extra,x+i eta_extra,y)]/(2Q),
    Q=integral rho w (x+iy)^2 real and nonzero.               (3)

It is bounded by ||eta_extra||_(L2(rho w))
times (integral rho w r^2)^(1/2)/|Q|. Small quadrupole cancellation
therefore amplifies non-affine response; small support alone is not a
relative-angle error bound. The same formula can be differentiated in
time and preparation parameters after the moving material map is retained.

In0144 the actual first-gradient field is
eta=exp(i epsilon kappa.x)[U+i epsilon chi+...]. Its local first
gradient contains BOTH U tensor kappa and grad chi. In0139's actual
initial compact core translation, eta is constant on the selected core.
Consequently (1) gives delta theta=0 there, including its non-affine
correction, although macroscopic curl U may be nonzero. The implication
"actual tag angle equals curl U/2 for every acoustic preparation" is
refuted by this initial-data example. The affine preparation remains an
admissible different route; its subsequent cell history must be derived.

## 2. Passive orientation has an independent material-label row

On an axisymmetric column u=V(r)e_theta+W(r)e_z, let
eta=g(r)e_theta with smooth compact radial support away from the axis,
independent of theta,z,t. Both eta and u are divergence-free and
[u,eta]=0: in cylindrical coordinate components both are linear
combinations of partial_theta and partial_z with coefficients depending
only on r. Hence the exact Lin velocity w=eta_t+[u,eta] is zero.
Its Euler energy and vorticity perturbations vanish identically.
Nevertheless its material quadrupole angle changes by the nonzero
weighted g(r)/r rotation of a nonaxisymmetric passive tag.

This is a material relabeling with a distinguishable passive tag, not
an optical Euler mode. It proves that unrestricted absolute passive-tag
angle cannot by itself carry a positive Euler locking energy. It does
not refute0142's explicit prepared two-dimensional Euler history/action:
that history fixes its initial label data and contains genuine velocity.
Nor does it refute accepted C-CST-009's declared Cauchy--Born action.
The next physical-frame route can use a core-shape director tied to
vorticity, or retain the extra label row explicitly.

## 3. The observable coupling survives the frame-coordinate issue

Retain0117's actual current identity. Under its complete-response isotropy
hypothesis, the optical shape-rate first moment vanishes. If a prepared
physical optical angle q has derived spin S=j qdot at the retained order,
the hybrid centroid-plus-ambient displacement U and point-fluid mean X
obey

    U=X-c curl q+O(grad^2), c=j/(2rho).                       (4)

The total density rho includes the ambient fluid. This is an observed
current map, not an added dynamical interaction. Actual time connections
or spin mismatch in0142/0147 enter (4) with their established error;
they cannot be discarded in an exact autonomous claim.

To identify the algebraic consequence once a common actual action is
licensed, take one transverse helicity component, curl=h (h real), and
an independently derived diagonal quadratic action in X,q with
mass diag(rho,j) and potential diag(mu h^2,kappa+d h^2).
The conditional premise is stated here:0144's exact separation of axial
Fourier sectors does not alone prove this autonomous action or isotropy.
Keep a general scalar frame response b by writing

    U=X-c h q, Phi=q+b h X,
    D=1+bc h^2,
    (X,q)=D^-1(U+c h Phi, Phi-b h U).                       (5)

For D!=0 this is an invertible field chart, without any change in physics.
Pull back BOTH full kinetic and potential matrices. In particular

    M_UU=(rho+j b^2 h^2)/D^2,
    M_UPhi=(rho c-jb)h/D^2,
    M_PhiPhi=(j+rho c^2 h^2)/D^2,
    K_UPhi=[mu c h^3-bh(kappa+d h^2)]/D^2.                 (6)

Thus b=1/2 removes the leading mixed kinetic term and gives the familiar
locking potential kappa(Phi-curl U/2)^2/2 at this spatial order.
For any b, however, at the zero-wave-number optical frequency
omega0^2=kappa/j,

    K_UPhi-omega0^2 M_UPhi=-(kappa/2)h+O(h^3).             (7)

The first-gradient dynamical coupling is independent of b. Directly,
the underlying pure optical history X=0 gives U/Phi=-j h/(2rho),
also independently of b. Setting b=0 does not eliminate the physical
centroid response; it moves the coupling entirely into mixed inertia.
Dropping that inertia would erase the observed current and falsify (7).
The exact determinant pencil is changed only by det(T)^2, so no field
renaming creates an eigenfrequency, positive stiffness or new branch.

Acoustic spin can add a second-gradient term to (4), and actual optical
shape-rate/current corrections add their licensed rows. They affect the
second-gradient diagonal coefficients and must be kept for a full
continuum match; they do not alter the displayed first-gradient invariant
when their stated order holds. No coefficient in (6) is fitted.

## 4. Route scope and continuation

Established: direct physical affine/non-affine tag-angle variation; the
distinct passive-label row; and the complete conditional field-chart
identity with a frame-independent observed optical/translation first jet.
Refuted route: universal identification of an arbitrary actual acoustic
tag with the polar macro frame, and positive unrestricted passive-label
locking, by the concrete zero-angle/zero-Euler-energy mechanisms above.
Neither verdict propagates to the parent objective.

The actual current coupling survives and remains a positive component of
the construction. The next open achievements are an actual vorticity-tied
core director or a controlled affine acoustic preparation, the common
acoustic/optical action including its retained state, and its isotropic
and stationary EPS realization. These are constructions, not new review
requirements imposed on accepted008..010 or frozen0142.
