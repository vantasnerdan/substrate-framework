# Full Euler column equation, exact exterior, and the threshold mode

All calculations below are direct Euler derivations, independent of the
unavailable Sun1995 body. They concern a translating axisymmetric excitation
on the unbounded smooth vortex-column background of0022.

## 1. Exact full equation and background linearization

In a frame translating at c, set u_r=-psi_z/r, u_z=psi_r/r,
u_theta=F(psi)/r and specific pressure

    p=B(psi)-(|grad psi|^2+F(psi)^2)/(2r^2).

The cylindrical Euler equations reduce exactly to

    Delta_*psi+F F'(psi)-r^2 B'(psi)=0,
    Delta_*=partial_rr-r^-1 partial_r+partial_zz.       (1)

For a column with axial velocity W(r), swirl V(r), define
q=W-c, L(r)=rV(r), psi_0'=rq. Where q!=0, the labels are determined by

    F(psi_0(r))=L(r),
    B'(psi_0(r))=W'(r)/r+L(r)L'(r)/(r^3 q(r)).          (2)

These relations come from actual radial Euler pressure balance
p_0'=V^2/r. They do not assign a material equation of state.
Writing f=psi-psi_0, the full linear equation is

    Delta_* f+Q_c(r)f=0,
    Q_c=2LL'/(r^3 q^2)-(W''-W'/r)/q.                 (3)

For the pure swirl column W=0 and c>0,

    Q_c=Phi(r)/c^2, Phi=(L^2)'/r^3.                  (4)

For smooth nonnegative compactly supported axial vorticity, L'=r w>=0,
so Phi>=0 and has compact radial support. Assume it is not identically zero.

The exact steady functional is

    S[psi]=int [(|grad psi|^2-F(psi)^2)/(2r)+r B(psi)]drdz. (5)

Its Euler-Lagrange equation is(1), after multiplying by-r. The relative
functional subtracts its background value and first variation. This is a
steady variational functional, not automatically the physical kinetic-energy
Lyapunov function on an arbitrary Euler orbit.

## 2. Replace an artificial wall by the actual Euler exterior

Choose a matching radius R beyond the vorticity support with a positive
radial margin. On the exterior F is constant and B is constant. For a small
perturbation retaining this label range, (1) reduces exactly, not just
linearly, to Delta_* f=0 on r>=R. The decaying Fourier solution for axial
wave number k!=0 and boundary trace g_hat(k) is

    f_hat(k,r)=g_hat(k) r K_1(|k|r)/(R K_1(|k|R)),
    f_r_hat(k,R)=-|k| K_0(|k|R)/K_1(|k|R) g_hat(k).    (6)

This supplies the actual exterior Dirichlet-to-Neumann map. Matching psi
and psi_r matches all velocity components; the common F/B and Bernoulli
formula then match pressure. R is an integration surface, not a vessel wall.
The minimized exterior contribution to(5) is

    (1/2)int T_R(k)|g_hat(k)|^2 dk/(2pi),
    T_R(k)=|k| K_0(|k|R)/(R K_1(|k|R)).                (7)

The small-wave-number behavior is

    T_R(k)=k^2[log(2/(|k|R))-EulerGamma]+o(k^2),       (8)

and the large-|k| behavior is T_R(k)~|k|/R. The low-frequency logarithm
comes from the same-fluid exterior kinetic integral. Replacing it by a
constant local k^2 coefficient loses the far-field dynamics.

## 3. A positive critical speed exists, but its mode is a threshold resonance

At k=0 the exterior solution is constant in r. On[0,R] the regular critical
mode therefore solves

    f_0''-f_0'/r+(Phi/c_0^2)f_0=0,
    f_0(r)=O(r^2) at0, f_0'(R)=0.                    (9)

The largest generalized eigenvalue is

    c_0^2=sup_f int_0^R Phi f^2 dr/r / int_0^R f'^2 dr/r. (10)

The finite interval weighted energy space with f(0)=0 embeds compactly into
the weighted potential norm: near zero |f(r)|<=r ||f'||_L2(dr/r)/sqrt(2),
while away from zero this is ordinary one-dimensional compactness. Thus(10)
is attained. Sturm comparison supplies a positive simple first mode. Because
(f_0'/r)'=-Phi f_0/(c_0^2 r)<=0 and f_0'(R)=0, f_0'>=0. For a nontrivial
smooth positive core it is strictly positive through an initial interval.
The exterior constant is f_0(R)>0.

The associated h=f_0/r is not L2(rdr) at infinity: the norm diverges
logarithmically. This is a zero-frequency threshold resonance, not an isolated
L2 radial eigenfunction on R3. Eliminating the exterior before using the
finite-interval gap is essential.

An exact Rankine-core check distinguishes the boundary mechanisms. With
L=Omega r^2 inside R and constant outside, Phi=4Omega^2 in the core. The
regular solution is r J_1(2Omega r/c). Actual exterior matching requires
J_0(2Omega R/c_0)=0; a rigid impermeable wall instead requires J_1(...)=0.
The critical speeds differ. This patch example checks the exterior mechanism;
the smooth-profile eigenproblem above does not assume a Rankine interface.

## 4. Exact quadratic nonlinearity on a smooth pure-swirl column

Let D=(r q)^-1 partial_r, q=-c. The coefficient J multiplying f^2/2 in(1) is

    J=D^2(FF')-r^2 D^2 B'
     =4/q^3 [(L'^2+L L'')/r^4-3L L'/r^5]
     =2 Q_c'(r)/(r q).                               (11)

In a uniform-vorticity interior it vanishes; the transition profile supplies
the nonlinearity. Its projection onto the positive critical mode is

    beta=(1/2)int J f_0^3 dr/r
        =(1/q)int Q_c'(r) f_0^3 dr/r^2.              (12)

At c=c_0, integration by parts using(9) gives

    int Q_c' f_0^3/r^2 dr
       =-(1/2)int [f_0^3/r^5] h (3h^2-8h+8)dr<0,
    h=r f_0'/f_0.                                   (13)

All endpoint terms vanish: f_0=O(r^2) at zero and is constant outside R.
Here h>=0 and 3h^2-8h+8>0 for real h; the nontrivial core makes the inequality
strict. Since q=-c_0<0, beta>0. Thus a real smooth monotone column has a
nonzero focusing quadratic coefficient without fitting a radial profile.
Equation(13) is an analytic identity to be checked by the exact oracle before
assigning its terminal verdict.

## 5. Concrete solitary-wave continuation still to execute

For c just above c_0, the critical radial eigenvalue is small and positive.
Eliminate the radial complement on[0,R] while retaining the nonlocal exterior
term(7). The exact scalar Schur symbol should have the structure

    m_mu(k)=a mu+b k^2 log(1/|k|)+O(k^2+mu^2), a,b>0,

with a projected quadratic nonlinearity beta A^2 and controlled higher terms.
The logarithm selects axial scale sqrt(log(1/mu)/mu), rather than the naive
mu^-1/2 KdV scale. This is a derived candidate scale, not yet an exact Euler
solitary wave.

The next analytic tasks are: construct the complement inverse uniformly in
k and mu with the axis regularity and boundary operator(7); prove the scalar
inverse converges after that logarithmic rescaling in a suitable Sobolev
operator norm; solve the even translation-fixed homoclinic problem by an
actual inverse/variational argument; and bound the reconstructed solution
in C1 strongly enough to preserve the background label ranges and smooth
exterior margin. Only then does(1) provide an actual whole-space Euler wave.
Stability, physical action/current and particle identity remain separate.
