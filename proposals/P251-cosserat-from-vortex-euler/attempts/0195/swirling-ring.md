# Global small-swirl ring and two genuine flux-action twists

## 1. Exact Euler equation and the finite-radius continuation

Fix the smooth 0186 core, its sufficiently large FINITE ring radius R,
and its actual full-poloidal positive optical contour. Write phi=psi/R;
the physical cylindrical angle is called varphi_phys to distinguish it
from this scaled streamfunction. For a real smooth law g and small real
epsilon, prescribe

    F(psi)=epsilon R g(psi/R),
    B'(psi)=-f(psi/R)/R,
    u=(-psi_z/r,F(psi)/r,psi_r/r).                    (1)

Here B=p/rho+|u|²/2 is the physical Bernoulli function and f is exactly
the fixed 0186 smooth nonnegative near-Rankine law. Direct curl and cross
product give

    omega=F'(psi)u+[-Delta*psi-F F']/r e_varphi,
    u cross omega=(Delta*psi+F F') grad psi/r².

Thus the full steady Euler equation, INCLUDING the swirl source, is

    -Delta*phi=(r²/R²)f(phi)+epsilon² g(phi)g'(phi).  (2)

It is not the no-swirl stream equation plus a free axial velocity. Define

    q_epsilon(r,phi)=f(phi)+(R²/r²)epsilon² g(phi)g'(phi).

With the exact full-space stream kernel K_R of 0186, solve

    phi=Ktilde_R q_epsilon(phi)+A-U x-U x²/(2R),
    integral q_epsilon(phi)=Gamma,
    integral x q_epsilon(phi)=0.                    (3)

Integrals are over the same fixed meridional disk; x=r-R. The two borders
now constrain the actual effective toroidal source q, not f alone. This
retains the same exterior velocity and actual translation determination.
At epsilon=0 the entire derivative of (3) in (phi,A,U) is the invertible
0186 finite-R border. The implicit-function theorem therefore constructs
the global family for |epsilon|<epsilon_R, with

    phi_epsilon-phi_0, A_epsilon-A_0, U_epsilon-U_0=O(epsilon²)

in every required fixed finite smooth norm. These bounds are at FIXED R;
their constants are not asserted uniform as R tends to infinity.

Choose g supported on a stream-value interval strictly inside the flat
positive part of f. Then all new sources are strictly interior to the
old active core. For sufficiently small epsilon, q remains nonnegative,
smooth and supported on the continued toroidal core. The full-space
extension is precisely

    phi_global=K_R q_epsilon-mu_epsilon-U_epsilon r²/(2R),
    mu_epsilon=L_R Gamma/(2pi)-A_epsilon-U_epsilon R/2.

Both mu_epsilon and U_epsilon stay positive by continuity at this fixed R.
As in 0186, the function is negative on the axis, at infinity and on the
outer meridional disk boundary; the homogeneous exterior maximum principle
excludes another positive source region. In the exterior g=0 as well.
Consequently (2) holds globally, not only in the IFT disk. A constant shift
of psi to set its axis value to zero is accompanied by the same shift of
the arguments of F and B. There is no singular axis filament.

Velocity and vorticity are C-infinity throughout R³. Vorticity has compact
toroidal support, the actual far field is -U_epsilon e_z, and the finite
rest-at-infinity energy has the steady-frame Hamiltonian
E-U_epsilon I_z, with the actual Euler impulse from 0186. No rigid rotor
energy or exterior wall has been added.

## 2. Actual straightened trajectories and the closed core

The meridional streamlines remain nested about a nondegenerate elliptic
center. Let I be their outward enclosed r dr dz area divided by 2pi and
theta the corresponding canonical angle. Fix its orientation so that

    d volume=dI dtheta dvarphi_phys,
    u_meridional=Omega(I) partial_theta,
    psi_I=Omega(I), phi_I=Omega(I)/R.               (4)

Omega is signed; near a maximum it is negative in this convention. Its
absolute value is bounded away from zero on the fixed core. Define

    H(I)=(2pi)^(-1) integral_0^(2pi) r(I,theta)^(-2)dtheta,
    V(I)=epsilon R g(phi(I)) H(I).

Solve the periodic, zero-mean equation

    Omega partial_theta b=F(psi)/r²-V(I)

and put varphi=varphi_phys-b(I,theta). The right side has exactly zero
theta average, so the solution is genuine and periodic. The coordinate
change preserves volume, and the COMPLETE particle dynamics is

    u=Omega(I) partial_theta+V(I) partial_varphi.     (5)

At I=0 the center circle has nonzero speed F(psi_c)/r_c when
epsilon g(phi_c) is nonzero. Its physical period is
2pi r_c²/|F(psi_c)|, and its transverse return multipliers are
exp(+-2pi i Omega(0)/V(0)). These are physical streamline quantities,
not the optical eigenfrequency or a relabeling clock.

On each sufficiently small surrounding torus with V>0, a varphi-return
has meridional rotation number rho_u=Omega/V. The return section carries
the ACTUAL flux form V(I)dI wedge dtheta. Its flux action (flux/2pi) is

    J_u(I)=integral_0^I V(s)ds,
    d rho_u/dJ_u=(Omega' V-Omega V')/V³.             (6)

In particular it is not the unweighted meridional area derivative.

## 3. Constructing nonzero twist without arithmetic assumptions

At epsilon=0 select a small inner interval about the elliptic center,
well separated from the tag annuli considered below. Prescribe a positive
g_c and freely choose its first two stream derivatives g_1,g_2 there.
Take g to be that quadratic polynomial in a smaller neighborhood and
smoothly cut it off before reaching the tag annuli. Positivity of g on
a smaller core neighborhood follows from g_c>0, irrespective of the
finite derivative choices. Since f is constant there and g is polynomial,
the interior elliptic equation and its solution are analytic locally.

At epsilon=0 the leading logarithmic derivative controlling (6) is

    T_u=Omega_0'/Omega_0-(g_1/g_c)phi_0'-H_0'/H_0.  (7)

The coefficient of g_1 is -phi_0'/g_c, nonzero by (4). Choose g_1
nonzero and different from its one exceptional value. The continued
quantity is T_u+O(epsilon²), so the actual streamline flux twist is
nonzero for all sufficiently small nonzero epsilon in a core neighborhood.

Vorticity has its own exact expression, not an assumed parallelism:

    omega=epsilon g'(phi) u+(f(phi)/R)partial_varphi,
    omega=alpha(I)partial_theta+beta(I)partial_varphi,
    alpha=epsilon g' Omega,
    beta=f/R+epsilon² R g' g H.                    (8)

Thus the SAME level tori and core circle are actual vortex-line tori and
a nonzero closed vortex line. Near the center beta>0. Their vortex-line
return rotation is rho_omega=alpha/beta, with flux action
J_omega'=beta. Because f is constant in this inner region,

    rho_omega'=epsilon R/f_c [g_2 phi_0' Omega_0
                                      +g_1 Omega_0']+O(epsilon³). (9)

The coefficient of g_2 is nonzero. Choose it to avoid the one exceptional
value in (9). Both actual flux-action twists are then nonzero. There is
no need to impose an unproved arithmetic property of a Bessel zero or
to numerically select a tiny spectral splitting.

At the core, rho_u(epsilon)=c_u/epsilon+O(epsilon) with c_u nonzero,
and rho_omega(epsilon)=c_omega epsilon+O(epsilon³) with c_omega nonzero.
On every sufficiently small positive interval their derivatives are
nonzero. Diophantine numbers have full measure; the inverse images under
these regular scalar maps therefore have full measure. Select epsilon
in their intersection, inside all geometric and spectral smallness bounds.
Both core returns are elliptic and nonresonant. By nonzero I derivatives,
one can likewise select ONE inner boundary action I_b for which both
rotation numbers are Diophantine. These are actual free profile/parameter
choices, not postulated frequency relations.

The integrable returns on that inner disk are analytic and preserve their
respective flux forms. Their nonzero flux-action derivatives are precisely
nonzero normal torsion. The invariant-circle persistence theorem and
flux-form straightening in [Enciso--Peralta-Salas, arXiv1210.6271,
Theorem 7.6 and the proof of 7.10](https://arxiv.org/pdf/1210.6271)
therefore apply to sufficiently close analytic divergence-free fields.
The selected core periodic orbits also persist by the nonunit transverse
return multipliers. This is geometric stability of each tube; it does not
assert that the two separately continued invariant tori coincide after an
arbitrary perturbation. It constructs neither arbitrary knots nor a
constant-curl approximant to this ring.

## 4. Tag-compatible swirl support, not an untransported painted tag

In the generic swirl case a finite-arc density initially
chi_0(I,varphi) evolves exactly as chi_0(I,varphi-V(I)t). Its angular
shape generally shears. Section 2 did not make it stationary.

For the present construction, choose the support cutoff of g strictly
INSIDE the innermost of the two positive 0186 tag annuli, leaving an open
margin to all three of their finite tag-control variations. Invariant
stream intervals persist as epsilon varies. On those complete annular
neighborhoods g and every derivative vanish, and hence the actual angular
velocity F/r² is exactly zero. There the straightening b may be zero;
the actual densities chi(I,varphi_phys) are exactly stationary because
both u.I and u.varphi_phys vanish. Their poloidal material circulation
is still retained. No boundary is fixed against a nonzero Euler normal
velocity and no off-shell tag constraint is substituted for transport.

The inner core and its twisted tori are disjoint from these outer tag
annuli, all within the SAME connected finite-core vortex ring. The tags
are not claimed to lie inside the selected inner KAM boundary. They
continue to observe the same global optical eigenmode through the full
Euler pressure. The next companion proves the positive pole and the
actual tag/current continuation on this precise geometry.
