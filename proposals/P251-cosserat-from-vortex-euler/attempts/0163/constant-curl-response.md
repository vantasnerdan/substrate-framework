# Full current and action of the constant-curl triangular lift

Use the unchanged planar field and whole-cell conventions of 0161, and
set W=-lambda psi=zeta/lambda, u=(v,W). Then curl u=lambda u and
p=-|u|^2/2. The physical axial Fourier convention is exp(ikz). This
calculation concerns the complete periodic horizontal cell, not a tube
subparcel or a frozen ambient fluid.

## 1. Exact Euler variables and mean

Let eta=curl_h w_h, b=w_z, r=b-eta/lambda, and retain the full Hodge map

    w_h=m+K eta+i k d b,  K=-J grad G, d=grad G,
    G=(-Delta_h)^(-1),   <eta>=<b>=0.

Pressure pi is the complete solution of

    (-Delta_h+k^2)pi
      =2 sum_ij (partial_i v_j)(partial_j w_h,i)
        +2 i k grad W.w_h.                               (1)

Its mean is -2<W b>; this vanishes in the excited vector representation,
but is retained in the general formula. The full horizontal mean is

    m_t=-i k <W w_h+v b>.                                (2)

The two exact local evolution identities are

    eta_t=-A eta-w_h.grad zeta
      -i k[W eta+(grad W cross w_h)_z-zeta b],
    r_t+A r=-i k S,
    S=pi+v.w_h+W b+W r,       A=v.grad.                  (3)

Here grad W=lambda Jv, so
(grad W cross w_h)_z/lambda=-v.w_h. At k=0 the horizontal
eta equation is precisely 0161's Euler equation and r is passively
advected. Thus the mean-zero group in the norm
||eta||_2+||r||_2 is bounded in the same vector representation: the
first part is the proved whole-cell Arnold/dual-row group, the second
is unitary transport. The pressure and Hodge terms give a bounded
O(k) perturbation on this fixed space. This uses neither an ambient
spectral gap nor a discarded circulation variable.

The first-shell Green identity gives, for ANY eta,

    <W K eta>=-<v eta>/lambda,
    <W d b>=-J<v b>/lambda.

Consequently (2) is exactly

    m_t=-i k <v r>-(k^2/lambda)J<v b>.                    (4)

For the complete translation T X=(-Dv X,-grad W.X), its two original
terms cancel:

    <W T_h X+v T_z X>=-<X.grad(Wv)>=0.                   (5)

One term alone would therefore give a false O(k) conclusion. Equations
(3)-(4), rather than that isolated term, determine the remaining sign.

## 2. A genuine slow tangent of the lifted Euler operator

Put q=sqrt(lambda^2+k^2). For each of the six first-shell modes with
horizontal wave vector beta, use its positive-helicity polarization for
n=(beta,k). If e1=J beta/lambda and
e2=(lambda e_z-k beta/lambda)/q, then e1+i e2 has curl eigenvalue q.
The sideband translation T_k X has coefficient

    (Psi lambda/2)(beta.X)(e1+i e2).                     (6)

At k=0 this is the actual coefficient of -X.grad u; the six-mode
definition handles opposite beta with their corresponding polarizations.
It has zero cell mean, is exactly solenoidal, and

    curl_k T_k=q T_k,  T_k=T+O(k),
    (T_k)_z=(lambda/q)T_z,
    r(T_k)=(1-q/lambda)(T_k)_z=O(k^2).                   (7)

The FULL linear Euler generator about a Beltrami field is

    L_k w=P_k[u cross (curl_k w-lambda w)].               (8)

Hence L_k T_k=(q-lambda)P_k(u cross T_k)=O(k^2) in the
fixed response space. This is an actual divergence/pressure-corrected
slow tangent, not a postulated cell oscillator. Its linear Bernoulli
pressure and S in (3) are bounded and analytic near k=0 in the selected
representation, with S(T)=0 and S(T_k)=O(k). The harmonic pressure term
in (1) creates no pole; the relevant scalar is in the vector sector.

Let X_t=m and write the actual solution w=m+T_k X+y, using the Hodge
identification for the zero-mean remainder y. The bounded k=0 group,
L_k-L_0=O(k), L_k T_k=O(k^2), and L_k m-T_k m=O(k)m
give, for epsilon=|k| and x=epsilon X,

    sup_{0<=t<=T/epsilon} (|x|+|m|+||y||)
       <=C_T(|x(0)|+|m(0)|+||y(0)||).                  (9)

As in 0161, Gronwall is applied AFTER controlling the fast groups; its
exponent is C T. Translation amplitudes O(1/k) are explicitly scaled.

## 3. The complete current fixes the actual acoustic sign

Use the same bounded material-cell coordinate R of 0161, A R=v.
Equations (3)-(4) imply the exact improved current identity

    p_c=m+i k<R r>,
    (p_c)_t=k^2<R S>-(k^2/lambda)J<v b>.                 (10)

For the full translation,

    <v T_z X>=lambda C_v JX.

Together with (7)-(9), this gives on slow time tau=epsilon t

    x_tau=p_c+O_T(epsilon),
    (p_c)_tau=+C_v x+O_T(epsilon),
    p_c-m=O_T(epsilon).                                (11)

All pressure, passive-r and fast-mode terms remain in the controlled
remainder of (10). The sign is POSITIVE, unlike the planar array.
In particular the actual common-V initial data have x(0)=0,y(0)=0 and

    m(t)=cosh(|k| sqrt(C_v)t)V0+O_T(|k|)|V0|,
    0<=t<=T/|k|.                                       (12)

Prepared displacement data w(0)=T_k X0 give the corresponding sinh
phase. Constants use the fixed cell units; the dimensionless error is
O(|k|/lambda). This is a controlled growing Euler response, not a
numerically inferred eigenvalue or an assertion about all times.

## 4. Independent same-material-action sign

Specify the displacement phase by actual initial material position
xi0=(X0,0)exp(ikz) and velocity w0=T_k X0. Its circulation data are those
of this preparation; they are not asserted to equal a fixed-Kelvin leaf.
The common-V phase has xi0=0 and pi0=rho(V0,0). Lin reconstruction gives

    pi_D=rho D_t xi0=rho(T_k-T)X0,
    B xi0=i k W(X0,0),  B=u.grad.

The complete Jacobi Hamiltonian is

    H_J=||pi||^2/(2rho)-Re<pi,B xi>
          +<xi.Hess(p_phys)xi>/2.

The periodic Hessian average is zero for xi0. The exact six-mode
calculation from (6) yields

    ||(T_k-T)X||^2=2lambda^2(1-lambda/q) X.C_v X,
    Re<(T_k-T)X,i k W(X,0)>=lambda k^2 X.C_v X/q.

Thus the actual phase energy is

    H_phase=rho |V0|^2/2-rho lambda(q-lambda)X0.C_v X0.   (13)

The symplectic mass pairing is exactly rho X0.V0; all additional
horizontal momentum rows have zero mean. There is no mixed initial
energy term. Equation (13) has leading stiffness -rho k^2 C_v,
independently confirming (11). The full moving-frame connection is
retained when using the observation chart (x,p_c), as in 0161. Neither
energy nor current is repaired by changing the sign of W alone: the
opposite-handed constant-curl lift has the corresponding same even sign.

## 5. Route verdict and failure-generated continuation

The actual constant-curl response and its action sign are established.
As a candidate for a POSITIVE axial acoustic stiffness, this route is
refuted by (11)-(13), with its mechanism named: axial-flow and
polarization/pressure reaction reverse the planar restoring sign.
This does not refute its local optical mechanism, generalized force-free
arrays, other directions, or the parent coupled-continuum objective.

In particular the generic-direction, whole-field SO(3) ensemble remains a
registered constant-curl candidate. The axial result alone does not fix
that tensor. The actual one-wave second-jet example in 0151 demonstrates
why signs of separate orientations cannot decide an averaged physical
response. Such a continuation must compute the shared action/current and
response, not average proposed frequencies.

The next registered repair keeps the same planar array and takes
W_C=sqrt(C+lambda^2 psi^2), C>0. This is an exact generalized-force-free
Euler field. In the specified axial frame moving at sqrt(C), its
h_C=W_C-sqrt(C) and all fixed spatial derivatives are O(C^(-1/2)).
A finite-C ordered comparison with the established planar acoustic
histories can therefore retain every pressure/shear term and choose C
last at each nonzero k. That continuation is not a fixed-C all-k
homogenization theorem and does not alter constant-curl canon.
