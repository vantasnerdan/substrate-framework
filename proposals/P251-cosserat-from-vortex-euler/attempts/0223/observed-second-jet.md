# Actual pressure-retaining physical second jet

This is the axial slice of the registered common-K construction, on the
same fixed C016 cell. Coefficients below are Taylor coefficients, so the
second coefficient is one half of the second derivative. No finite-wave
label formula is imported. The exact source is `axial-kelvin-operator.md`.

## 1. Explicit stationary operator coefficients

Let Q remove the normal constant, f=H^-1 h, z=s+h,
C=div(z grad psi), and C_sym as defined in the source. The exact scalar
generator has L(k)=L0+k L1+k² L2+O(k³). Its components are

    L0(h,s)=(-Th, H^-1 T(s+h)-Ts),

    (L1 q)_h=i Q[psi(s-h)+H^-1 C],
    (L1 q)_s=-i H^-1{div[psi grad(f-s)]+C_sym f},

    (L2 q)_h=-H^-1 T(s-f),
    (L2 q)_s=-H^-1 T f.                              (1)

The Q in the first line is the actual constant pressure row. H^-1 C is
already mean free. Expanding only nonzero normal modes of (H+k²)^-1
gives (1); no inverse acts on a removed Bloch harmonic. In particular
there is no inverse small axial wave number in these coefficients.

For actual initial coefficients q_j(0), the full histories are

    q0(t)=E(t)q0(0),
    q1(t)=E(t)q1(0)+int_0^t E(t-r)L1 q0(r) dr,
    q2(t)=E(t)q2(0)+int_0^t E(t-r)[L1 q1(r)+L2 q0(r)] dr,
    E(t)=exp(t L0).                                  (2)

These are actual Euler/Kelvin response columns, not a postulated local
oscillator. For each selected smooth preparation, ordinary fixed-time
linear Euler estimates and the analytic mean-free pressure multiplier
bound the remainder by C_(T,preparation) |k|³ in the finite Sobolev and
time norms needed by the observations. The constant includes every
selected initial correction. A diagonal limit requires its actual cost.

## 2. Physical observations, including the second-order flux

Keep chi_psi distinct from the product chi psi. Define linear rows

    A(s)=Q_tag^-1 int ab chi_psi Ts,
    B(h)=-i Q_tag^-1 int ab chi_psi grad(H^-1 h).grad psi,
    C(s)=-rho int (2chi+chi_psi Dpsi)s,
    D(h)=-i rho int chi_psi Rpsi H^-1 h,

    F(h,s)=rho int [chi psi r.grad s
                                  +(2Fchi-chi Dpsi)h],
    J(h)=i rho int chi psi r cross grad(H^-1 h).      (3)

Here Q_tag is the fixed reference quadrupole; it is not the mean-free
projection Q in (1). Centering is unchanged because the exact C2 sector
has zero tag centroid. The literal response coefficients are

    theta_j=A(s_j)+B(h_(j-1)),
    G_j=C(s_j)+D(h_(j-1)),
    S_j=partial_t G_j+i F(q_(j-1))+i J(h_(j-2)),     (4)

with negative-index rows zero. Equations (3)--(4) follow from the actual
transported density and moving spin, not from canonical momentum.
In particular a coherent theta_2 is insufficient to establish a matched
physical spin: the retained F(q1) and J(h0) enter S2. The symmetric
material shape moments are the same linear density functional with ab
replaced by the appropriate centered quadratic polynomial.

For clarity, the current identity has a local, integration-independent
form. Let xi denote its normal component, a0=2Fchi-chi Dpsi, and use
the actual Lin rate. Then

    S_density-partial_t G_density
       -ik rho[chi psi(r cross xi)+a0 h]
      =rho div_perp[chi v(r cross xi)+a0 xi].         (5)

The compact tag removes the boundary term. This exposes precisely where
omitting either axial transport or displacement-current motion changes
the result. The equality uses div_perp xi=-ik h, not planar
incompressibility. It also retains the initial G.

## 3. The actual compact forcing after the principal secular term

For q0=(-s,s), s_t=-Ts, set z1=h1+s1 and initially leave its preparation
free. Then

    h1_t=-Th1+2i Q(psi s),
    z1_t=-B0 T z1+i Fz s,   B0=1-H^-1,
    Fz=2 Q psi+H^-1 div[psi grad(H^-1+1)]
                                      +H^-1 C_sym H^-1. (6)

On the nonzero streamline harmonic, <psi s>=0. The principal radial
symbol of Fz is multiplication by psi, not 2psi. Indeed

    Fz s=Q(psi s)+R s,
    R s=H^-1[-grad psi.grad s+psi s
                 +div(psi grad H^-1 s)+C_sym H^-1 s], (7)

with the mean subtraction retained as in (6). On this selected harmonic
the displayed quantities have the compatible zero means; the general
formula is (6). The explicitly computed secular leading response is

    h1=2it psi s,  z1=it psi s,  s1=-it psi s,        (8)

up to actual homogeneous initial terms and the remainder forced by
(7) and H^-1 T(psi s). Formula (8) is a principal decomposition, not
an exact Doppler theorem. Substitution into (6) gives the exact residual

    (partial_t+B0 T)(it psi s)-i Fz s
                    =-it H^-1 T(psi s)-i R s.       (9)

Thus an exact pressure-retaining correction is forced by
i R s+it H^-1 T(psi s). The second summand is a different compact
operator and cannot be absorbed into a scalar frequency average.

On a fixed regular streamline annulus with coordinate E and angular
harmonic l, write s=f(E)e^(il theta_o). The only order-minus-one term
of R is factored explicitly as H^-1(a(E,theta_o) f'(E)e^(il theta_o));
all remaining terms are order-minus-two operators on f(E), with smooth
fixed-annulus coefficients. This makes the input derivative cost
visible. It does not make R s physically negligible for a narrow band:
its observation against a fixed smooth tag can have a nonzero limiting
value even when its norm is lower order than the leading field norm.

## 4. What the actual preparation must close

An initial band-center correction contributes the actual derivatives
of the source with respect to E; it is not a pressure replacement.
Additional Kelvin initial q1(0),q2(0) enter exactly through (2).
The finite-output problem is to make (4), including its F and J rows,
equal the desired physical angle/spin/current jet on the fixed window.
An approximate Sylvester graph for (9) is a sufficient competing
construction, but the entire Euler complement need not be stationary.

The original same-cell acoustic D,V data and off-tag returns contribute
their own cross observations and forms. The axial slice has
n dot (K cross U)=0 and therefore cannot establish the nonzero acoustic
rotation coefficient required by the generic coupled parent. At zero
wave, mean-free optical generators and their canonical momenta have
zero pairing with the uncorrected common uniform translation/velocity
columns. Their first and second cross jets are actual integrals, not
zero by this zeroth-order fact. Root0228's algebraic normalization can
match forms once these rows are supplied; it cannot create a missing
material-angle response.

This source establishes the exact inherited spatial response and
physical-current interface. The controlled autonomous observed second
jet and generic-K acoustic observation remain the active construction;
the full parent campaign is not closed by retaining an infinite state.
