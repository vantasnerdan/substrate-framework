# Physical coupled normalization and the small-inertia margin

This is the stationary reference of the actual moving action, not a claim
that its time-dependent remainder vanishes. All coefficients below are
obtained after coherent action averaging as in `coherent-action.md`.
In particular B_T is the resulting equation-ratio curvature, not an
average of individual squared frequencies. Its actual spatial supplier
and clock-error bounds remain explicit.

## 1. Both quadratic forms undergo the physical field map

Write the transverse reference action in actual mean X and registered
optical angle q, with h the real curl-helicity eigenvalue, h²=k²:

    M0=diag(rho+m_a2 k², j+j2 k²),
    K0=diag(rho a k², j omega²+(j B_T+omega² j2)k²).

The circular-reference physical map from0164 is

    U=(1-I k²/(4rho))X-j h q/(2rho), Phi=q+h X/2.       (1)

I is the literal acoustic tagged polar moment density. It is NOT the
optical phase inertia j. The other static acoustic shape filters change
the scalar gradient mass below; they do not change the first-curl map.
Pull both forms back with the inverse of(1), retaining their second jets:

    M_UU=rho+m_U k², m_U=m_a2+I/2-j/4,
    M_U,Phi=0 through k²,
    M_Phi,Phi=j+(j2-j²/(4rho))k²;
    K_UU=(rho a+j omega²/4)k²,
    K_U,Phi=-j omega² h/2,
    K_Phi,Phi=j omega²
                 +[j B_T+omega² j2-j² omega²/(2rho)]k². (2)

Thus dropping gradient inertia while retaining the potential would give
the wrong curvature. Let D=diag(rho,j), M2 the gradient-mass matrix in(2).
The explicit derivative observation map

    (U,Phi)_physical=(Id-D^(-1)M2 k²/2)(U_N,Phi_N)       (3)

gives mass D through this order, and changes the potential as well. It is
invertible for sufficiently small k at every selected finite j>0.
Equations(1),(3), and the literal current rows specify what U_N,Phi_N mean;
no unobserved canonical angle is substituted for the registered tag.

The resulting standard coefficients are

    mu=rho a, alpha=j omega²/4,
    C_T=j B_T-j² omega²/(4rho), C_L=j B_L.              (4)

The transverse pencil is precisely

    K_N=[[(mu+alpha)k², -2alpha h],
         [-2alpha h, 4alpha+C_T k²]], M_N=diag(rho,j).

Its acoustic speed squared is a. Its optical gap is omega² and optical
curvature is C_T/j+alpha/rho=B_T, as required by the invertible original
map. The longitudinal spin has curvature B_L; incompressibility removes
longitudinal displacement, not longitudinal spin.

Positive B_T therefore permits, but does not alone imply, positive C_T:

    0<j<4rho B_T/omega².                               (5)

The finite-cell construction makes this a genuine choice of actual
nonzero packet density. It does not supply or fit an inertia. Neither I
nor j2 creates another positivity condition at this spatial order.

## 2. Local curvature energy and its explicit boundary improvement

For a purely axial leading scalar curvature B, fourth-Haar averaging gives
B_T=B/5 and B_L=3B/5. Before normalizing the curl mass, its raw local energy
is jB[(tr G)²+2||sym G||²]/10, with G_ij=partial_j Phi_i.
The curl mass normalization contributes -j² omega²||skew G||²/(4rho).
That raw representative need not be pointwise positive on arbitrary G.

For periodic macro fields, or with its surface flux retained, add the
explicit null Lagrangian eta[(div Phi)²-tr(G²)]. The identity is

    (div Phi)²-partial_i Phi_j partial_j Phi_i
      =div[Phi div Phi-(Phi dot grad)Phi].              (6)

It changes neither C_T nor C_L. Choose the equivalent coefficients

    c_s=c_a=C_T/2, c_tr=(C_L-C_T)/2.

These have positive pointwise curvature energy when C_T>0 and
3C_L>2C_T, which holds for the axial ratio and is stable under small
properly scaled transverse errors. The boundary couple traction changes
by the flux in(6); this is recorded, not a new microscopic energy.

## 3. Moving physical maps retain their connections

For actual time-dependent I, shape filters, optical clocks and integrated
spin, replace(1) by the complete Q(t,K) rows in `coherent-action.md`.
The derivative map(3) is then time dependent as well: its time derivative
and second derivative contribute to the mixed current and potential.
The exact phase pullback includes these terms. Formula(4) names only the
stationary reference; it does not erase them. Translation momentum is
the measured Euler mean with its hybrid current, and angular momentum
includes the literal tag spin, initial displacement dipole and shape
current. Equality to the normalized canonical momenta is a derived
reference relation plus these displayed errors, not an identity by name.
