# Intrinsic-swirl Euler carrier: moment, action, tail, and variational control

## 1. Axisymmetric conventions and physical observables

Let the constant mass density be `rho0`, and use cylindrical coordinates `(r,theta,z)`.
For the Cao--Zhan sign convention write

    u = (-Psi_z e_r + xi e_theta + Psi_r e_z)/r,
    xi = r u_theta,
    zeta = omega_theta/r = L Psi,
    dnu = r dr dz,
    L = -(1/r) partial_r((1/r) partial_r) -(1/r^2) partial_zz.

For a decaying field, `Psi=K zeta`. Direct cylindrical integration and the weak
definition of `K` give

    H = (rho0/2) integral_R3 |u|^2 dx
      = pi rho0 integral_Pi [zeta K zeta + xi^2/r^2] dnu,              (1)

    P_z = (rho0/2) integral_R3 (x cross omega)_z dx
        = pi rho0 integral_Pi r^2 zeta dnu,                            (2)

    j_z = rho0 integral_R3 (x cross u)_z dx
        = 2 pi rho0 integral_Pi xi dnu.                                (3)

Thus swirl `xi`, not merely ring circulation `integral zeta dnu`, supplies the
centered axial intrinsic angular momentum. For Cao--Zhan, (1) converges because
the lab disturbance decays as `|x|^-3`, (2) because vorticity is compact, and
(3) because the azimuthal velocity is compactly supported. The full absolute
integral `integral |x cross u|` need not converge for a dipolar poloidal tail;
no such stronger claim is used. All three integrals converge absolutely for the
compact Gavrilov velocity.

If `xi=F(Psi)` and boundary terms vanish, then

    Hel(u) = integral_R3 u dot omega dx = 4 pi integral_Pi xi zeta dnu. (4)

This follows because the meridional contribution integrates to a second copy of
the azimuthal contribution. Equation (4) is a classical helicity identity; no
quantum interpretation is attached to it.

## 2. The actual axisymmetric Euler leaf

For a smooth axisymmetric time-dependent Euler solution,

    D_t xi = 0,
    D_t zeta = (1/r^4) partial_z(xi^2),                                 (5)

where `D_t=partial_t+u_r partial_r+u_z partial_z`. The meridional flow preserves
`dnu`. Consequently, with sufficient decay,

    C_C = integral C(xi) dnu,
    C_F = integral zeta F(xi) dnu                                     (6)

are conserved for every smooth one-variable `C,F`; the source in (5) becomes a
`z` boundary term in `d C_F/dt`. In particular, circulation `integral zeta`,
angular momentum `integral xi`, and helicity `integral xi zeta` are genuine
conserved-label quantities. This is the invariant leaf relevant to stability,
not an arbitrary variation of `zeta` alone.

Normalize energy and axial impulse by

    E = H/(2 pi rho0)
      = (1/2) integral [zeta K zeta + xi^2/r^2] dnu,
    I = P_z/(2 pi rho0) = (1/2) integral r^2 zeta dnu.                  (7)

An energy--impulse--Casimir functional is

    A = E - c I + integral [C(xi)+zeta A0(xi)] dnu.                    (8)

Its first variations are

    delta_zeta A = K zeta - c r^2/2 + A0(xi),
    delta_xi A   = xi/r^2 + C'(xi)+zeta A0'(xi).                       (9)

If the moving-frame streamfunction is `psi=K zeta-c r^2/2` and
`xi=F(psi)`, choose

    A0'(xi) = -1/F'(psi),
    C'(F(psi)) = -B'(psi)/F'(psi).                                    (10)

Then (9) is exactly

    zeta = -B'(psi) + F(psi)F'(psi)/r^2,                              (11)

the Bragg--Hawthorne constitutive relation. This is an exact Hamiltonian/leaf
bridge, not just an analogy between two elliptic equations.

For perturbations `(eta,chi)=(delta zeta,delta xi)`, its Hessian is

    Q[eta,chi] = integral {eta K eta + 2 A0'(xi) eta chi
                    + [1/r^2+C''(xi)+zeta A0''(xi)] chi^2} dnu.        (12)

Nonlinear orbital control would follow from a coercive sign of (12), modulo
symmetry zero modes, on the dynamically accessible tangent space. Neither
source paper proves that result.

## 3. Cao--Zhan/Turkington generator: exact Schur-complement bridge

Cao--Zhan Theorem 1.3 has

    psi_e = K zeta_e -(U_e r^2)/2 - mu_e,
    U_e = W log(1/e),
    xi_e = H(psi_e)/e,
    zeta_e = e^-2[-B'(psi_e)+H(psi_e)H'(psi_e)/r^2].                  (13)

The Turkington generator singled out in their discussion is

    H(psi)=psi_+,              -B'(psi)=alpha 1_{psi>0}.               (14)

On the active core, (13) is the Euler critical point of the *conserved*
functional

    A_e = E - U_e I - mu_e integral zeta dnu
          +(alpha/e) integral xi dnu - e integral zeta xi dnu.        (15)

Indeed,

    delta_zeta A_e = K zeta-U_e r^2/2-mu_e-e xi,
    delta_xi A_e   = xi/r^2+alpha/e-e zeta,                            (16)

and setting these to zero gives exactly (13)--(14). Every term in (15) is an
Euler invariant from (6)--(7). Its Hessian completes a square:

    Q_e[eta,chi]
      = integral [eta K eta +chi^2/r^2-2e eta chi] dnu
      = integral [eta K eta-e^2 r^2 eta^2] dnu
        + integral (chi/r-e r eta)^2 dnu.                              (17)

Eliminating `xi` pointwise from (15) gives

    xi_* = e r^2 zeta-alpha r^2/e,

    min_xi [xi^2/(2r^2)+(alpha/e)xi-e zeta xi]
      = -(r^2/(2e^2))(e^2 zeta-alpha)^2.                               (18)

With the positivity constraint this is precisely Cao--Zhan's

    -(1/e^2) J(r,e^2 zeta),
    J(r,s)=(1/2)r^2(s-alpha)_+^2.                                     (19)

Therefore their reduced functional

    E_e(zeta)=(1/2) integral zeta K zeta
       -(U_e/2) integral r^2 zeta -(1/e^2) integral J(r,e^2 zeta)      (20)

is exactly the Schur complement of a physical conserved Euler functional for
this generator. Its second variation on the active set,

    delta^2 E_e[eta]=integral [eta K eta-e^2r^2eta^2] dnu,             (21)

is the first term of (17).

This is useful variational control, but it also identifies the unproved step.
Cao--Zhan maximize (20) over a relaxed box/mass/support class. The eliminated
`xi` direction in (18) is a minimum, and (17) has manifest positive pure-swirl
directions. Hence a reduced maximum is not automatically a full Euler
energy--Casimir extremum. If (21) has a strict negative direction, (17) is a
saddle on the unconstrained joint space. The required stability calculation is
the sign/coercive gap of (17) on the *dynamically accessible* tangent satisfying
the full family of linearized constraints (6), after removing the axial
translation zero mode, followed by the non-axisymmetric sectors. The source's
explicit statement that stability/nonstability and local uniqueness remain open
is consistent with this exact remainder.

## 4. Same Cao--Zhan family: time transfer, moment, tails, and action

The theorem's steady field `v_e` tends to `-U_e e_z`. Galilean invariance gives
the actual finite-energy lab solution

    q_e(t,r,z)=v_e(r,z-U_e t)+U_e e_z,                                 (22)

which tends to zero and translates the core at speed `U_e`. At `t=0`, its
decaying streamfunction is `K zeta_e` (additive constants are irrelevant), so
(1)--(3) apply to this same time-dependent family.

Let `r_*=kappa/(4 pi W)`. Cao--Zhan prove

    integral zeta_e dnu=kappa,
    integral r^2 zeta_e dnu -> kappa r_*^2.                            (23)

Consequently

    P_z(e) -> pi rho0 kappa r_*^2
            = rho0 kappa^3/(16 pi W^2).                               (24)

For the actual swirl generator (14),

    j_e = (2 pi rho0/e) integral (psi_e)_+ dnu >0.                     (25)

Using their core diameter `Theta(e)`, bounded rescaled support, and the
subsequential `C^1_loc` profile convergence
`Psi_e(y)=psi_e(X_e+e y)->Psi(y)`,

    j_e/e -> 2 pi rho0 r_* integral_R2 Psi_+(y)dy >0.                  (26)

Thus every fixed member has finite nonzero intrinsic angular momentum, but the
thin-ring limit drives it to zero. In the same limit,

    H_e = (rho0 kappa^2 r_*/2) log(1/e)+O(1),                          (27)

because `integral zeta K zeta = [kappa^2 r_*/(2pi)]log(1/e)+O(1)` and
the swirl-energy term is `O(1)`. The kinematic helicity obeys

    Hel(q_e)=Hel(v_e)=4pi integral xi_e zeta_e dnu,
    e Hel(q_e) -> 4pi r_* integral_R2 Psi(y) phi(y)dy >0               (28)

along the same convergent subsequence (`phi=e^2 zeta_e` in core coordinates).
The Galilean addition does not change the integrated helicity because
`integral omega_e dx=0`. Hence, along each profile-convergent subsequence with
the displayed positive limits, this family has the notable joint scaling

    j_e ~ c_j e,       Hel_e ~ c_H/e,       H_e=Theta(log e^-1),       (29)

not a finite nonzero universal action in its desingularized limit.

Let `Ivec=(1/2) integral x cross omega dx=P/rho0`. Compact vorticity and the
Biot--Savart expansion give, uniformly outside a ball containing the core,

    q_e(x,t)=[3n(n dot Ivec)-Ivec]/(4pi |x|^3)+O(|x|^-4),              (30)

with `Ivec -> pi kappa r_*^2 e_z`. The pressure of the lab field has the
corresponding finite-energy quadrupolar tail

    p_e(x,t)=[(3n_i n_j-delta_ij) M_ij]/(4pi |x|^3)+O(|x|^-4),
    M_ij=integral q_{e,i}q_{e,j} dx,                                  (31)

provided the displayed first moment is finite. These are actual Euler tails:
there is no radial-swirl pressure patch and no relabeling of the nonsteady
0001 construction.

At the centered impulse frame, the finite axial rotation moment is `j_e n`.
Its connected rigid-rotation stabilizer is `SO(2)` (absent an extra discrete symmetry), so the
rigid orientation orbit is `SO(3)/SO(2)=S^2`. With the 0005 physical KKS
convention, `r_a=a cross x`, `[r_a,r_b]=-r_{a cross b}` and
`Omega_m(X_a,X_b)=-<m,[r_a,r_b]>=(a cross b) dot (j_e n)`. Therefore

    Omega_rot = j_e sin(theta) dtheta wedge dphi,
    integral_S2 Omega_rot = 4 pi j_e.                                  (32)

There is no single global symplectic potential on this sphere. A regular
two-chart representation is

    Theta_N=j_e(1-cos theta)dphi,
    Theta_S=-j_e(1+cos theta)dphi,
    dTheta_N=dTheta_S=Omega_rot,
    Theta_N-Theta_S=2j_e dphi.                                        (32a)

The equatorial transition integral is `4pi j_e`, exactly the sphere period;
this is the action measure before any conditional prequantization rule.
Only the finite stabilizer-axis pairing enters this pulled-back two-form; no
absolute full-space angular-momentum integral for the poloidal tail is claimed.
The statement is restricted to the rigid-rotation suborbit; the full
volume-preserving-diffeomorphism stabilizer is not classified. Also `pi_1(S^2)=0`,
so orientation alone supplies no fermionic rotation loop.

## 5. Gavrilov: a compact steady twisted carrier with nonzero `j`

Gavrilov's exact local solution uses cylindrical radius `s`:

    a=alpha(s/R,z/R),       p=R^4 a/4,
    b=(R^3/4) sqrt(H(a)),
    u=(p_z e_s-p_s e_z+b e_phi)/s.                                    (33)

The Taylor data in the paper give `a=2[(s/R-1)^2+(z/R)^2]+O(3)` and
`H(a)=4a+O(a^2)`, so `b>0` on a sufficiently small positive-pressure shell.
For a smooth nonnegative nonzero cutoff `omega_c(p)` supported in
`[delta,2delta]`,

    u_tilde=omega_c(p)u,       d p_tilde=omega_c(p)^2 dp               (34)

is an exact smooth steady Euler field with compact toroidal support. Both
meridional terms in (33) are retained; they supply the axial pressure balance
that the discarded radial compact-swirl ansatz could not satisfy.

The same-field observables are

    j_G = 2pi rho0 integral omega_c(p)b(p) dnu >0,                     (35)

    H_G = pi rho0 integral omega_c(p)^2
             [|grad p|^2+b^2]/s^2 dnu
        = 3pi rho0 integral p omega_c(p)^2 dnu < infinity,             (36)

where the last equality uses Gavrilov's `|u|^2=3p`. Furthermore,

    Hel(u_tilde)=integral omega_c(p)^2 u dot curl(u) dx
                =4pi integral xi_tilde zeta_tilde dnu,                (37)

with `xi_tilde=omega_c(p)b(p)` and the streamfunction
`Psi_tilde=-integral^p omega_c(s)ds` in the convention of section 1.
Equation (37) is exact; its sign and nonzero lower bound depend on the chosen
cutoff/profile and are not asserted by the existence theorem.

Here both velocity and vorticity vanish outside the toroidal shell, and pressure
is constant there, so all velocity and pressure-gradient tails vanish exactly.
The nonzero positive swirl gives `L=j_G e_z`. The connected rigid stabilizer is
again `SO(2)`, and therefore

    Omega_G = j_G sin(theta)dtheta wedge dphi,
    integral_S2 Omega_G=4pi j_G.                                      (38)

This is the strongest positive carrier in 0010: an actual smooth compact steady
Euler field, with real meridional circulation, intrinsic twist on pressure tori,
finite energy, and finite nonzero rotational KKS action.

The cutoff freedom is not a stability proof. It produces a family of equilibria
on different conserved-label leaves. Applying (8)--(12) requires inverting the
profile `xi=F(Psi)`; at the flat cutoff boundary `F'` degenerates, and no source
proves a globally regular Casimir or a coercive Hessian. Establishing such a
leafwise coercive extension (or an instability mode) is the concrete remaining
LP2 construction for this family.

## 6. Prescribed helical/knotted Beltrami comparator

Enciso--Peralta-Salas construct `curl u=lambda u` with prescribed knotted vortex
tubes and robust KAM line structure. Their global field is explicitly not
finite-energy; the theorem gives only `O(1/|x|)` decay. Consequently the physical
Hamiltonian is not finite, `Hel=lambda integral |u|^2` is not finite, and no
finite angular-momentum/KKS period follows. The stability in their theorem is
structural persistence of vortex trajectories and Lyapunov stability of a tube
boundary under the fixed vector field, not orbital stability under the Euler
PDE.

This refutes only that particular global Beltrami realization as the requested
finite-action particle carrier. It does not refute compact helical Euler fields:
Gavrilov supplies a compact integrable twisted carrier, while prescribed-knot
topology plus compact support remains a separate construction.

## 7. Classical scaling and what the action result does not select

For any fixed compact steady template `u_0`, Euler scaling gives another exact
steady solution

    u_{U,L}(x)=U u_0(x/L),       p_{U,L}(x)=U^2 p_0(x/L).               (39)

Writing shape constants `E_0=(1/2)integral|u_0|^2`,
`J_0=integral(y cross u_0)_z`, one has

    H_{U,L}=rho0 U^2 L^3 E_0,
    j_{U,L}=rho0 U L^4 J_0,
    integral_S2 Omega=4pi j_{U,L}.                                    (40)

If a dimensionless internal streamline has period `tau_gamma`, its physical
turnover period is `T_gamma=(L/U)tau_gamma`, and therefore

    H_{U,L} T_gamma=(E_0 tau_gamma/J_0) j_{U,L}.                        (41)

Equations (40)--(41) are a genuine classical action/period bridge. They also
show the missing scale-selection mechanism exactly: Euler leaves `U` and `L`
free, so `j` can be varied continuously. The rotational Hamiltonian is constant
on the isotropic orientation sphere, hence (38) alone gives no nontrivial
orientation dynamics. A conditional prequantization equation may restrict an
already chosen `j`; it neither chooses `N=1` nor derives a Hilbert/Born sector,
exchange statistics, relativity, or electric/weak currents.

## 8. Route-scoped verdicts and next executable dependency

- **Cao--Zhan/Turkington-generator route — established as a classical route.**
  It supplies an exact time-dependent translating finite-energy Euler family,
  finite nonzero `j_e` for every fixed thickness, finite impulse, explicit
  dipole/quadrupole tails, nonzero helicity, and the rigid KKS sphere. The new
  exact result is the conserved functional (15) and its Schur-complement identity
  (17)--(21). Its thin limit has `j_e -> 0`, and its PDE orbital stability is not
  established.

- **Gavrilov compact-twist route — established as the strongest positive
  carrier.** It supplies a smooth, compactly supported, steady Euler field with
  real meridional flow, finite energy, strictly nonzero finite axial angular
  momentum, and finite KKS period. Its helicity sign and nonlinear orbital
  stability remain profile-dependent/open.

- **Global prescribed-knotted Beltrami route — refuted for finite-action use.**
  The named construction has infinite energy and only vortex-line structural
  stability. This is not a no-go for other compact helical constructions.

The next executable physical dependency is not “derive quantum mechanics from
prequantization.” It is LP2 on one of the two actual carriers: characterize the
dynamically accessible tangent space and prove a uniform coercive modulated
energy estimate (or exhibit the exact growing mode) for (17) or the Gavrilov
analogue, including non-axisymmetric perturbations. Success would leave the
continuous Euler scaling (40) as the next distinct scale-selection mechanism,
followed by a two-carrier exchange configuration and quantum/relativistic
current map.
