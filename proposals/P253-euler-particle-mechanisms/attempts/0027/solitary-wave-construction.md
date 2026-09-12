# Same-fluid exterior reduction and logarithmically scaled Euler solitary wave

This is a direct candidate proof, to be reviewed at its full operator and
nonlinear reconstruction scope. The algebraic inputs in exterior-construction.md
passed structure-check.exit=0. No Sun cylinder theorem is silently invoked.

## 1. Freeze one smooth background and the physical claim

Use fixed background length and velocity units ell_0,U_0. All coordinates,
speeds, streamfunctions and operators below are nondimensional; physical
units are restored by x_phys=ell_0 x, u_phys=U_0 u and
t_phys=(ell_0/U_0)t. Thus log L_mu in(3) has a dimensionless argument.
These units set no fitted particle constant and are held fixed as mu varies.

Take w(r^2) smooth, nonnegative, constant and strictly positive near r=0,
zero for r>=R_0, with a nontrivial smooth transition. Fix R>R_0 with a margin.
The column has U=V(r)e_theta, L=rV=int_0^r s w(s^2)ds and positive circulation.
There is no external wall, density interface, gravity or surface tension.

The claim sought here is a nontrivial smooth axisymmetric traveling Euler
solution U+v(r,z-ct), for c in a one-sided interval above its critical speed,
with v in H^m(R3) for every finite m and finite absolutely defined excess
kinetic energy above U. It is an exact solitary excitation, not a nonlinear
orbital-stability theorem. Persistence of its exact traveling profile does
not imply stability of an open set of neighboring fields.

## 2. Finite-interval operator with the exact exterior attached

Let H=L2((0,R),dr/r), V the completion in int f'^2/r of smooth f with
f(0)=0. Trace at R is continuous on V. At c=c_0, define

    A_0=-partial_rr+r^-1 partial_r-Q_c0,

with the natural condition f'(R)=0. Its closed quadratic form on V is
semidefinite with a simple normalized positive kernel f_0. The compact
embedding V->H gives discrete spectrum on this finite interval; its
restriction to f_0^perp has a positive gap. This finite-interval statement
is compatible with, and different from, the whole-space threshold resonance.

For c>c_0 close, A_c has first eigenvalue lambda(c)>0 with

    lambda'(c_0)=2 int Phi f_0^2 dr/r / c_0^3>0.

Use mu=lambda(c) as parameter. Let f_c be its positive H-normalized first
eigenfunction. The remaining eigenvalues stay uniformly separated from zero.
For physical axial Fourier wave number k, attach the exact exterior:

    L_c(k)=A_c+k^2 I+T_R(k) b tensor b,
    b(f)=f(R), T_R(k)=|k|K_0(|k|R)/(R K_1(|k|R)).      (1)

This is a form operator V->V*. It is positive and satisfies
L_c(k)>=(mu+k^2)I in H. The nonlocal boundary term is not dropped.

An exact rank-one inverse makes its singular component explicit. Write
P_c=f_c tensor f_c, Q_c=I-P_c (Q_c here is a projection, distinct from the
radial coefficient Q_c(r)),

    R_q(k)=[Q_c(A_c+k^2)Q_c]^-1 on Q_c H,
    q_k=R_q(k)Q_c b, d_k=b(q_k), h_c=f_c(R),
    e_k=f_c-[T_R(k)h_c/(1+T_R(k)d_k)]q_k,
    m_mu(k)=mu+k^2+T_R(k)h_c^2/(1+T_R(k)d_k).

Then the full inverse, as H->V, is exactly

    L_c(k)^-1=R_q(k)-[T_R(k)/(1+T_R(k)d_k)]q_k tensor q_k
                         + e_k tensor e_k/m_mu(k).    (2)

The trace functionals are interpreted through the variational inverse, not
as L2 delta functions. In particular, d_k>=0 and the denominator is>=1.

## 3. Logarithmic scale and the scalar limit

Choose the large root L_mu>sqrt(e) of

    mu L_mu^2=h_0^2 log L_mu, h_0=f_0(R)>0.           (3)

It exists uniquely on that branch for sufficiently small mu. Set X=z/L_mu,
and K_mu(kappa)=mu L_c(kappa/L_mu)^-1. The claim needed for nonlinear
continuation is the operator-norm convergence

    K_mu -> K_0=P_0/(1+kappa^2)

in both the uniform H->V multiplier norm and its first kappa derivative. (4)

For bounded kappa, (2), the smooth eigenpair in c, and
T_R(k)=k^2[log(2/(|k|R))-EulerGamma]+o(k^2) give

    m_mu(kappa/L_mu)/mu ->1+kappa^2,
    e_(kappa/L_mu)->f_0,
    mu[R_q-T q_k tensor q_k/(1+Td_k)]->0.             (5)

The same convergences hold after one kappa derivative, including kappa=0:
kappa log|kappa| is bounded near zero and its coefficient tends to zero.
One must control the entire Fourier line, rather than only (5).

Here is the required large-frequency control. The radial form and bounded
potential give, with constants uniform in c near c_0,

    ||L_c(k)^-1||_(H->V)
       <=C[(mu+k^2)^-1/2+(mu+k^2)^-1].               (6)

For |kappa|>=L_mu^delta with any fixed 0<delta<1, (6) multiplied by mu
tends to zero uniformly, because
k^2/mu=kappa^2/(h_0^2 log L_mu). For M<=|kappa|<=L_mu^delta, k->0;
q_k is uniformly bounded in V, d_k is bounded, T_R(k)->0, and
log(L_mu/|kappa|)>=(1-delta)log L_mu. Formula(2) therefore bounds
||K_mu|| by C/(1+kappa^2)+o(1), uniformly in this middle region.
First take M large, then mu small. This proves the zeroth-derivative part of(4).

For the first derivative, the exterior minimization identity gives exactly

    T_R'(k)=2k int_R^infinity |f_k(r)|^2 dr/r,
    0<=k T_R'(k)<=2T_R(k), k>0,                     (7)

where f_k has boundary value1. Hence, in quadratic forms,

    0<=partial_k L_c(k)<=2L_c(k)/k.

Differentiating the inverse and using (6) gives

    ||partial_kappa K_mu||_(H->V)
       <=(C/|kappa|) mu[(mu+k^2)^-1/2+(mu+k^2)^-1]. (8)

This controls the far region. In the middle region, differentiate(2): the
bounded complement and small-k Bessel expansions give
C/(|kappa|(1+kappa^2))+o(1). On bounded kappa, differentiating the explicit
small-k expansions gives convergence to partial_kappa K_0. This completes(4).
The endpoint kappa=0 is treated by continuity, not division by zero in(8).

## 4. Nonlinear equation in a fixed Banach space

Choose s>1/2 sufficiently large and let

    Z^s={F: F and XF in H^s(R_X;V)},
    Y^s={G: G and XG in H^s(R_X;H)}.

Use their real even-X subspaces. The first derivative in(4) is exactly what
makes Fourier multipliers converge Z/Y with this spatial weight:
Fourier(XKG)=i K' G_hat+K Fourier(XG). Thus

    K_mu -> K_0 in bounded operators Y^s->Z^s.        (9)

The background labels F_c(psi),B_c(psi) from the direct Euler construction
are smooth in psi and c on a fixed interval containing the relevant ranges.
Because w is constant near the axis, the swirl label is linear in psi and
B_c' is constant there. Because w is zero beyond R_0, the labels are constant
on the exterior range. Extend the interior linear label smoothly past psi=0
if needed. Small solutions will remain in the physical range after the
regularity estimate below. The exact nonlinear remainder N_c(r,f) in

    L_c(D_z) f=N_c(r,f)                              (10)

is therefore quadratic, smooth, and supported in a fixed annulus away from
r=0 and r=R when ||f||_V is small. Its leading coefficient is J_c(r)/2.
On that annulus V embeds into bounded continuous functions. Sobolev product
estimates in X and Taylor's integral remainder give C1 convergence on bounded
subsets of Z^s:

    n_mu(F)=mu^-2 N_c(r,mu F)
        ->n_0(F)=J_c0(r) F^2/2,  Z^s->Y^s.           (11)

The weight is carried by one factor of F, so it does not change this estimate.
The nonlinear equation for f(r,z)=mu F(r,z/L_mu) is exactly

    F-K_mu n_mu(F)=0.                               (12)

There is no dropped finite-core term in(12).

## 5. Solve the limiting and perturbed equations

At mu=0 the range of K_0 lies in span f_0. Set F=f_0 A(X). With

    beta=(1/2)int_0^R J_c0 f_0^3 dr/r>0,

proved in exterior-construction.md, equation(12) becomes

    A-A''=beta A^2.

Its even positive homoclinic is exactly

    A_*(X)=3/(2beta) sech^2(X/2).                    (13)

The full derivative of the limiting map is identity on the radial complement,
with a triangular coupling into the critical component. The remaining scalar
operator is I-(1-partial_XX)^-1 2beta A_*. Its associated Schrödinger operator

    -partial_XX+1-3sech^2(X/2)

has one zero mode A_*', which is odd. Its even sector is invertible: rescaling
X=2y gives one-quarter of the l=3 Pöschl-Teller operator with shifted bound
values -5,0,3 and continuum starting4; the zero mode is excluded by parity.
Equivalently, reduction of order from A_*' proves the second zero solution
is not decaying, so no even zero mode exists; Fredholmness follows from the
decaying multiplication potential and the smoothing free inverse.
The inverse preserves the one spatial weight because the potential decays
exponentially and the resolvent away from its zero sector has exponentially
decaying Green functions. Thus the derivative is an isomorphism on even Z^s.

The C1 convergence (9),(11), and the implicit-function perturbation theorem
now give F_mu in even Z^s with

    F_mu -> f_0 A_* in Z^s,
    f_mu(r,z)=mu F_mu(r,z/L_mu),                       (14)

solving the complete equation(10) for all sufficiently small mu>0. Smooth
parameter differentiability at mu=0 is not required: an invertible derivative
and C1-small perturbation of the map suffice.

## 6. Physical reconstruction and localization

The exterior is reconstructed by(6) of exterior-construction.md. The weighted
Z^s control gives L1_z(V) as well as L2_z(V). Consequently all core swirl
perturbations are L1 in z and the cross kinetic term with the background is
absolutely integrable. The exterior poloidal perturbation has finite energy
by its exact minimized form(7); there the swirl label is constant, so there
is no background-swirl cross term. Thus the literal excess kinetic energy
exists, without canceling infinite background energies.

At the axis use f=r^2 phi. The equation becomes the radial-four-dimensional
part of a five-dimensional elliptic equation. The V/H solution has precisely
the locally finite H1/L2 radial weights for phi; coefficients are smooth in
r^2, and the nonlinear forcing is supported away from the axis. Interior
elliptic regularity therefore gives smooth phi and uniform O(mu) bounds on
f/r^2 and f_r/r, using arbitrarily high fixed s and physical z derivatives
L_mu^-1 partial_X. This preserves u_z=-c+f_r/r<0 for small mu, the label
ranges, and the exterior margin. It also ensures smooth Cartesian velocity
through the axis. At r=R the exact Dirichlet-to-Neumann relation matches
psi_r; elliptic transmission then supplies smoothness across that artificial
integration surface. Bernoulli reconstruction matches pressure there.

The laboratory field is obtained by adding c e_z to the translating-frame
velocity and replacing z by z-ct. It solves the full constant-density Euler
system and is a localized excitation above U. Its nontriviality follows from
(14) and the nonzero homoclinic amplitude. No prescribed time-history source,
rigid wall, external surface force or quantum constant is present.

## Status and next achievements

This is the candidate exact existence proof to audit, particularly (4),(11)
and the axis/weighted reconstruction. The scalar and operator algebra will be
checked through reusable APIs before its verdict is frozen. Even if the full
existence proof passes, nearby-state restoring control, same-family interaction,
physical action/current, spin/statistics, relativity and electron/neutrino
identification remain open. Exact persistence of one solitary wave is a useful
object; it is not the robust particle requested by the parent campaign.
