# Direct EPS physical-angle theorem

## Object, conventions, and imports

Use an actual field supplied by Enciso–Peralta-Salas, *Existence of knotted
vortex tubes in steady Euler flows*, [1210.6271](https://arxiv.org/abs/1210.6271),
Theorem 1.1 and its local-to-global construction (Theorem 8.3). Thus
`curl u0=lambda u0`, `omega=curl u0`, and an invariant material solid torus
`D` exists. Choose an interior point of the tube where `omega != 0`.
This proof operates directly on that R³ field with the R³ Leray operator P;
it does not infer a compact-domain Green operator from local approximation.

For the estimates below fix a length unit and write coordinates, lambda,
and k in that unit. All inequalities are dimensionless inequalities in those
coordinates. The final H and B are the physical integrals, with physical
generators; they have units energy and action. No length unit is a modulus.

Compact smooth divergence-free generators act by pushing forward vorticity:

    v_xi = P(xi cross omega),
    delta omega = curl(xi cross omega),
    H(xi,eta) = rho integral[v_xi.v_eta - v_xi.curl(v_eta)/lambda],
    Omega(xi,eta) = rho integral omega.(xi cross eta).

These are the physical orbit energy Hessian and KKS form, not the unconstrained
material Jacobi Hessian. Indeed `delta² omega=curl(xi cross delta omega)`;
the second energy variation is
`rho integral[|v|²+u0.(xi cross curl v)]` and `omega=lambda u0`
reduces it to H. The Euler generator X=u satisfies `i_X Omega=dE`, since
`u.(eta cross omega)=omega.(u cross eta)`. This fixes the sign.
Compact generators make these Hessian integrals finite even though the EPS
background need not have finite total kinetic energy.

## 1. A physical core angle, rather than a relabeling coordinate

Orient a coordinate axis so `omega_z >= wstar > 0` on a closed ball inside
D. Inside a smaller ball choose a point x0 and a unit axis e perpendicular
to omega(x0). With r=x-x0 and a smooth bump c equal to one near x0, define

    xi_R = curl[-c(r) |r|² e/2].

It is compact, divergence free, and equals `e cross r` near x0. Consequently

    xi_R(x0)=0,
    delta omega_R(x0) = (omega.grad)xi_R = e cross omega(x0) != 0.

Its amplitude q is exactly the infinitesimal physical tilt angle of the
vorticity direction at the tube core point. A circular cross-section does
not erase this observable. The claim is this local direction/rotation jet,
not rigid rotation of every point of a finite cross-section. The compact
extension connects this core motion to fluid in the same tube.

## 2. A compact negative-helicity cage and its exact KKS pairing

Choose a nonzero smooth real bump phi in an annulus disjoint from xi_R,
still inside the positive-omega_z ball. For signed k with k/lambda>0 set

    p1=(cos(kz), sin(kz),0), p2=(-sin(kz),cos(kz),0),
    xi_i = -curl(phi p_i)/k = phi p_i - grad(phi) cross p_i/k.

Both generators are compact and exactly divergence free, hence preserve D
and leave its boundary unchanged. Their cross product is the exact identity

    xi_1 cross xi_2 = phi² ez + (phi/k) J grad_perp(phi)
                     + phi_z grad(phi)/k²,
    J(vx,vy,0)=(-vy,vx,0).

Using `curl omega=lambda omega` and compact support gives

    B/rho = (1-lambda/(2k)) B0 + T/k²,
    B0 = integral omega_z phi² > 0,
    T = integral phi_z omega.grad(phi).

Let `Tstar=||omega||_infinity ||phi_z||_2 ||grad phi||_2`, with the supremum
on the cage support. For K=|k| satisfying `K>=|lambda|` and
`K²>4 Tstar/B0`, the exact finite-k pairing obeys `B>rho B0/4`.
For `eta_1=xi_1+xi_R`, `eta_2=xi_2`, B is unchanged exactly: the added cross
product vanishes pointwise by disjoint supports. Projection tails do not
enter KKS; this is a local generator pairing.

## 3. Full finite-k energy estimate, including Leray and cutoff effects

Write the exact force fields as

    xi_i cross omega = Re[(F0_i + F1_i/k) exp(ikz)],

where F0_i and F1_i are fixed compact smooth complex amplitudes obtained
by expanding the displayed curl formula. They do not depend on k.
Let Pi=I-ez ez^T. The principal projected fields are

    v0_1 = a (sin(kz),-cos(kz),0),
    v0_2 = a (cos(kz), sin(kz),0),  a=phi omega_z.

Put A=integral a²>0 and G=||grad a||_2. For every real t=(t1,t2),
`v0=t1 v0_1+t2 v0_2` satisfies

    ||v0||_2²=A |t|²,
    integral v0.curl(v0)=-k A |t|²,
    ||curl(v0)||_2 <= (K sqrt(A)+G)|t|.

The middle identity is exact: the amplitude-derivative term is orthogonal
to v0 pointwise. Thus principal H is `rho(1+k/lambda) A I`.

Here is a usable finite bound on the discarded terms. For every Fourier
frequency eta, the Euclidean Leray symbol obeys

    ||P(eta+k ez)-Pi||_op <= 4 |eta|/K.

For |eta|<=K/2 this follows by comparing the directions of eta+k ez and
k ez; for |eta|>K/2 use the unit bound for two orthogonal projections.
Plancherel therefore bounds the projected F0 error by
`4 ||grad F0||_2/K`. The F1 term is bounded by `||F1||_2/K`.
Also curl P=curl, so the curl of the F0 error is the curl of
`(I-Pi)F0 exp(ikz)`; its fast z derivative vanishes because this vector
is vertical. The F1 curl has a bounded fast part and a `1/K` slow part.
The following deliberately generous constants consequently suffice for K>=1:

    D = sum_i [4 ||grad F0_i||_2 + ||F1_i||_2],
    E = sum_i [2 ||grad F0_i||_2 + ||F1_i||_2 + ||curl F1_i||_2].

For `v=t1 P(xi_1 cross omega)+t2 P(xi_2 cross omega)` and w=v-v0,

    ||w||_2 <= D |t|/K,   ||curl w||_2 <= E |t|.

Expanding both terms of H and using Cauchy–Schwarz proves

    |H_cage(t,t)-rho(1+K/|lambda|) A |t|²| <= rho C0 |t|²,
    C0 = 2 sqrt(A) D + D²
         + [D sqrt(A)+D G+sqrt(A) E+D E]/|lambda|.

No asymptotic coefficient is substituted for the exact finite-k Hessian.
The O(1) remainder has the explicit constant C0. Complex-amplitude norms
bound the corresponding real-part norms, so no uncounted averaging factor
appears in D or E.

## 4. Attach the physical core jet without losing positivity

Let `vR=P(xi_R cross omega)` and

    HR=H(xi_R,xi_R),
    CR=||vR||_2 + ||curl vR||_2/|lambda|.

All are finite. Curl is self-adjoint for these decaying smooth tangents,
so a cross entry can be bounded without differentiating the fast cage:

    |H(xi_R,t1 xi_1+t2 xi_2)|
       <= rho CR (sqrt(A)+D) |t|,    K>=1.

For the actual physical pair eta_1=xi_1+xi_R, eta_2=xi_2, this gives

    H_eta(t,t) >= rho [(1+K/|lambda|) A - Ctotal] |t|²,
    Ctotal=C0+|HR|/rho+2 CR(sqrt(A)+D).

Choose any finite K with

    K>=1,  K>=|lambda|,
    K²>4 Tstar/B0,
    K>|lambda| Ctotal/A,
    k=sign(lambda) K.

Then the FULL exact two-coordinate Hessian is positive definite and
`Omega(eta_1,eta_2)>0`. Every threshold is computed from the declared
background, cutoffs, and Leray projection. None uses a desired modulus,
gap, or Cosserat comparison value. A continuum carrier on R³ has no integer
constraint; integer carriers on an axial circle can be chosen above the
same finite bound after their period is fixed.

## 5. The exact positive angle action and scope

Evaluate the finite-k integrals `hij=H(eta_i,eta_j)`, `B=Omega(eta_1,eta_2)`.
The quadratic orbit action in the declared two-generator ensemble is

    L2 = B s qdot - (h11 q² + 2 h12 q s + h22 s²)/2.

Eliminate s, not an independently supplied mechanical mass:

    s=(B qdot-h12 q)/h22,
    L_angle = I qdot²/2 - Kangle q²/2 + a total derivative,
    I=B²/h22>0,
    Kangle=det(H)/h22>0,
    omega_angle²=det(H)/B²>0.

The cage and core react through the same Euler energy and KKS form. At x0,
the cage generators vanish and q gives the physical core tilt calculated
above; s is its same-fluid conjugate. There is no imposed strain reservoir
or postulated angle modulus. Exact integral coefficients, rather than a
numerical value for an unspecified EPS field, are the theorem's output.

**Route verdict: established as stated.** Every actual smooth constant-lambda
EPS tube with a nonzero interior vorticity point supports this explicit
compact two-generator positive physical-angle Euler orbit action. This is
an exact constrained quadratic-action statement, consistent with the
parent's prescribed microscopic/affine ensemble; it is not a claim that
arbitrary microscopic Euler initial data remain in a two-dimensional
invariant subspace. It establishes the positive microscopic angle-action
obligation directly on the EPS field. The common-angle/translation action
and spatial coarse-graining coefficients remain distinct parent
constructions, not silently supplied by this local theorem.
