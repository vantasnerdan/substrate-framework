# The fixed-lambda leading jet on a general thin tube

This local derivative calculation complements the independently completed
axisymmetric persistence construction. It does not supply the missing
fixed-lambda arbitrary-knot return twist.

On a Bishop-frame chart write

    X(s,y)=gamma(s)+epsilon*y_i E_i(s),
    t'=kappa_i E_i,  E_i'=-kappa_i t,
    B=1-epsilon*kappa·y,
    ds_phys^2=B^2 ds^2+epsilon^2(dy_1^2+dy_2^2).

Take the oriented frame (t,E1,E2), let v=v^s partial_s+w_i partial_yi,
and f=B^2 v^s. Its metric dual is f ds+epsilon^2 w_i dyi. The exact
curl and divergence equations are

    partial_1 w_2-partial_2 w_1=lambda*f/B,
    partial_2 f-epsilon^2 partial_s w_2=lambda*B*epsilon^2*w_1,
    epsilon^2 partial_s w_1-partial_1 f=lambda*B*epsilon^2*w_2,
    partial_s(f/B)+div_y(Bw)=0.                         (1)

The unit-disk boundary condition is w·y=0. EPS 6.8, for a fixed finite
bound on |lambda|, and the harmonic estimates give a uniform C^k bound
on w in these coordinates. Thus (1) gives grad_y f=O(epsilon^2).
The integrated divergence equation gives a cross-section-independent
flux integral f/B dy. Since integral_D B^(-1)dy=pi+O(epsilon^2), and the
prescribed harmonic component fixes the leading flux to pi under the
source's unit axial normalization, one obtains

    f=1+O(epsilon^2).

The transverse equations then give curl_y w=lambda+O(epsilon),
div_y w=O(epsilon), and w·y=0. The disk div-curl boundary problem has no
nonzero harmonic tangent vector field: its homogeneous kernel vanishes
by the simply connected Neumann harmonic-potential argument. Its fixed-
domain elliptic estimate therefore yields

    w=(lambda/2) J y+O(epsilon).                        (2)

All remainders hold in each fixed finite derivative norm, with constants
depending on the chosen geometry and eigenvalue bound. Bishop frames
need only be local; their normal holonomy is retained in global return
maps. Equations (1)-(2) do not assume a global zero-holonomy frame.

The physical field is u=(f/B)t+epsilon*w_i E_i. At the reference core its
leading velocity-gradient matrix, ordered (t,E1,E2), is

    Du = [[0, kappa^T], [kappa, (lambda/2)J]]+O(epsilon).

The physical transverse rotation is indeed lambda/2; it is not merely
the torsion of a moving coordinate system. Subtracting the Bishop frame
connection gives the particle coordinate matrix

    B_particle = [[0, 2*kappa^T], [0, (lambda/2)J]]+O(epsilon).

In the EPS Frenet convention the transverse block is
(lambda/2+tau_EPS)J. The leading particle return angle is therefore
lambda*L/2+integral tau_EPS ds. An elliptic strict-margin choice gives a
periodic core by the usual return-map implicit function theorem. It does
not itself make the boundary circle KAM robust.

For the adjoint covector normalized by k·u=1 the leading normal equation
is k_h'=(lambda/2+tau_EPS)J k_h-2*kappa. The geometric coupling is retained:
one cannot simply set k parallel to u on a curved tube. The true Euler
material amplitude subtracts lambda*k cross a/|k|^2 from the particle
matrix. For a large-eigenvalue geometric regime the oscillatory periodic
covector equation is the appropriate next estimate; a bound based only
on the norm of its forcing would lose the useful inverse-lambda scale.
The completed axisymmetric candidate makes this equation constant and
establishes that scale directly, while also proving its boundary twist.

`route_verdict: established` for equations (1)-(2) and the leading jet.
`evidence_scope: FIXED_EIGENVALUE_LOCAL_THIN_TUBE_ASYMPTOTICS`.
The boundary-twist/global arbitrary-knot assertion is not included.
