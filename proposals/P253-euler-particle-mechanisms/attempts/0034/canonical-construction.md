# Exact axisymmetric canonical Euler representation and its axial end data

Use dnu=r dr dz, bracket{f,g}=(f_r g_z-f_z g_r)/r, swirl labelxi=r u_theta,
poloidal vorticityzeta=omega_theta/r andpsi=K zeta. Material density isrho_m,
so all reduced actions and energies carry the common multiplier2pi rho_m.
The local calculations below are valid on a regular label stripxi_r>0;
flat-core/exterior limits require the indicated physical continuation, not
inversion of a constant label.

## Canonical equations recover actual Euler, including its swirl source

Introduce a scalar beta with

    zeta={xi,beta},
    H[xi,beta]=1/2 int[zeta K zeta+xi^2/r^2]dnu.        (1)

For compactly supported variations the exact cyclic integration identity
int a{b,c}=int c{a,b} yields

    H_beta={psi,xi}, H_xi=xi/r^2+{beta,psi}.           (2)

The physical canonical action is

    A=2pi rho_m int [int xi beta_t dnu-H]dt,           (3)

and its symplectic form is2pi rho_m int d beta wedge d xi, consistent with
i_X Omega=dH andOmega=-d(int xi d beta). Hamilton's equations are

    xi_t=-{psi,xi}, beta_t=xi/r^2+{beta,psi}.           (4)

Jacobi and the derivation property of the bracket give

    zeta_t+{psi,zeta}={xi,xi/r^2}=2xi xi_z/r^4.         (5)

Equations(4)-(5), the whole-space Hodge relation and the physical swirl
velocityxi/r are exactly axisymmetric incompressible Euler. No approximate
amplitude Hamiltonian has been substituted. Pressure is reconstructed by
Euler/Leray from that field. The canonical coordinates parameterize physical
fields through(1); beta->beta+f(xi) is a gauge, so a coordinate phase is not
a new observable degree of freedom.

## The end jump contains a real leaf invariant

At fixed time use(xi,z) as coordinates on a regular strip, writingr=r(xi,z).
Then

    {xi,beta}=(xi_r/r) partial_z beta|xi,
    dnu=(r/xi_r)d xi dz.                              (6)

Thus any given smooth physicalzeta admits the exact local lift

    beta(xi,z)=beta_-(xi)+int_-infinity^z
                   zeta(xi,s) r(xi,s)/xi_r(xi,s) ds,  (7)

whenever the indicated integral converges. Its end jump is

    J(xi)=beta_+(xi)-beta_-(xi)
         =int_-infinity^infinity zeta r/xi_r dz.       (8)

For any smoothD supported in the regular label interval,

    int zeta D(xi)dnu=int D(xi)J(xi)d xi.              (9)

These are the actual axisymmetric Euler mixed Casimirs. They identify the
leaf's distribution, not just one helicity number. A lift with the same
zero beta gauge at both axial ends exists only whenJ=0. NonzeroJ is allowed
in(7) and cannot be erased by a single beta->beta+f(xi) gauge.

This is why a localized velocity does not automatically provide a globally
decaying canonical potential. The finite-action domain must fix the end-jump
distribution or retain its conjugate boundary variables. Equations(1)-(5)
are exact for compact variations irrespective of this issue; a global action
with moving end data requires the surface terms below.

## Physical translation momentum and the boundary term

On a finite meridional domainD, the normalized axial impulse is
I_D=int_D(r^2/2)zeta dnu. Integrating by parts withzeta={xi,beta} gives

    I_D=int_D beta xi_z dnu
        -int_boundaryD (r^2/2) beta d xi,             (10)

with the oriented meridional boundary taken counterclockwise in(r,z).
The sign follows from f dxi wedge dbeta=beta df wedge dxi-d(f beta dxi);
neither a discarded end jump nor a periodic assumption is implicit.
If the boundary term vanishes or is constant on the fixed leaf, the canonical
bulk translation generator isP=int beta xi_z dnu. Its Hamiltonian equations
arebeta_t=deltaP/dxi=-beta_z andxi_t=-deltaP/dbeta=-xi_z; it translates the
physical field in the negative coordinate direction. With this convention,
H-cP hasbeta_t=H_xi+c beta_z,xi_t=-H_beta+c xi_z and represents the frame
moving at speedc. The physical multiplier2pi rho_m occurs in both momentum
and symplectic form and cancels in the field vector.

Before applying this to the exact solitary wave, evaluate(8),(10) on its
actualF/B labels and prove convergence. The gauge/end distribution is a
physical leaf condition; it is not fixed by demanding a favored scalar
amplitude equation. The full solitary linearization will retain perturbations
of the poloidal field and advected label with the appropriate fixedJ row.

## The actual solitary wave has a nonzero end-jump distribution

Let a denote the background radius label, xi=L(a), and write the exact
solitary streamsurface asr=R(a,z). Since its frame streamfunction is
psi_frame=-c a^2/2, the0027 perturbation satisfies

    f(R(a,z),z)=c[R(a,z)^2-a^2]/2.                    (11)

On every compact regular label annulus, smallness and the radial implicit
function theorem giveR_a>0. Its exact vorticity and jump integrand are

    zeta=LL'/(c a)[1/a^2-1/R^2],
    zeta R/xi_r=L/(c a)[R R_a/a^2-R_a/R].            (12)

There is no assumption that the wave belongs to the zero-jump column leaf.
Indeed, withf=mu F_mu(r,z/L_mu), F_mu->f_0(r)A_*(z/L_mu),

    J(L(a))
      =12 L(a)f_0(a)/(beta c_0^2 a^4) mu L_mu
          +o(mu L_mu).                              (13)

The error is uniform on every fixed compact regular annulus. To see this,
first the implicit relation(11) givesR-a=f(a,z)/(c a)+O(f^2+|f f_r|)
there. Expand the bracket in(12): its linear term is2(R-a)/a^2, while
terms involvingR_a-1 begin at quadratic order. The weighted scaled
convergence in0027 controls the z integral. Its interior radial C1 version
follows by solving the scaled radial Euler equation forF_rr: the axial
term isL_mu^-2 F_XX, coefficients are smooth on the annulus, and the
nonlinearity ismu times a smooth quadratic term. Interior one-dimensional
elliptic estimates, iterated with the weighted X Sobolev bounds, upgrade
F_mu convergence to radial C1 on a smaller annulus. Finally
int A_* dX=6/beta gives(13). Quadratic errors areO(mu^2 L_mu), and the
remaining linear error iso(mu L_mu).

All coefficients in(13) are strictly positive. Hence the exact nonzero
solitary family hasJ>0 on every chosen regular annulus for sufficiently
smallmu. A globally equal-end decaying beta chart is therefore refuted for
this family. The two-end chart(7) is the exact repair. This statement does
not refute a physical Euler action: it identifies the end/Casimir data that
that action must retain.

## Full solitary linearization and the real modulation companions

Let(psi_c,zeta_c,xi_c) be the exact wave in laboratory variables, and put
Psi_c=psi_c-c r^2/2. Direct differentiation of Euler gives

    eta_t=-{Psi_c,eta}-{K eta,zeta_c}
                     +(2/r^4) partial_z(xi_c chi),
    chi_t=-{Psi_c,chi}-{K eta,xi_c}.                  (14)

This retains the full Hodge kernel, vorticity-gradient and advected-label
terms. Translation invariance gives the exact zero mode

    L_c partial_z(zeta_c,xi_c)=0.                    (15)

Where the exact branch is differentiable in c, parameter differentiation
also gives

    L_c partial_c(zeta_c,xi_c)
                    =-partial_z(zeta_c,xi_c).         (16)

Equation(16) is on the full phase space. Its speed derivative changes the
mixed-Casimir distribution(13), and need not lie in a fixed coadjoint leaf.
Therefore a phase/speed Jordan companion cannot be inserted into or deleted
from the fixed-leaf operator by notation. The modulation construction must
state whether it fixesJ(xi), retains its variation, or varies the underlying
F/B label functions along with the wave. The same distinction applies to
other relative invariants carried by the branch.

The next propagation estimate applies to(14) on that declared space. Its
column limit is0030's oscillator, but the localized terms couple axial
Fourier modes and the low-frequency inverse loses|k|. Those terms and the
end-jump row remain the named dynamical remainder; no scalar stability
conclusion has been substituted.

## Exact material Routh reduction identifies what the canonical action means

Use material meridional labels(a,b), azimuththeta_0 and an axisymmetric
volume-preserving map

    X=(r(a,b,t), theta_0+vartheta(a,b,t), z(a,b,t)),
    r det[d(r,z)/d(a,b)]=a.                          (17)

On a finite material domain, or for compact material variations, the full
physical material action, including its incompressibility constraint, has
kinetic part

    S=2pi rho_m int dt int a da db
                     [r_t^2+z_t^2+r^2 vartheta_t^2]/2. (18)

Every material angle is cyclic. Its conserved momentum is exactly
xi(a,b)=r^2 vartheta_t. Holding that distribution fixed and subtracting its
momentum-times-angle-rate gives the exact Routhian

    R=S-2pi rho_m int xi vartheta_t a da db dt
     =2pi rho_m int [r_t^2+z_t^2-xi^2/r^2]/2 a da db dt. (19)

The centrifugal acceleration from(19) is+xi^2/r^3; the volume constraint
supplies the full pressure force. This is an exact symmetry reduction of
Euler, with no elastic constitutive approximation. Its Hamiltonian restores
the positive swirl energyxi^2/(2r^2), exactly as in(1).

For the whole-space reference column, take the reference material map
`r=a,z=b` and hold the same label momentum `xi(a,b)` fixed. The finite object
is the relative Routhian

    R_rel=2pi rho_m int dt int a da db
      [(r_t^2+z_t^2)/2-xi^2/(2r^2)+xi^2/(2a^2)],     (19a)

together with the relative volume-constraint term. The last summand is fixed
label data. It changes neither the material equations nor the centrifugal
sign, but it subtracts the logarithmically and axially divergent reference
column contribution. Equation(19a) is used on maps for which the displayed
relative integral and the relative constraint pairing converge. The
first-order action(3) is likewise a fixed-end relative action, or an action
for compact variations, rather than an unrenormalized absolute whole-space
integral.

The first-order canonical action(3) describes this reduced meridional/swirl
system. It is therefore not automatically equal to the unreduced material
action(18): fixing a cyclic momentum changed the endpoint term in(19).
Likewise beta is a coordinate for the advected-label deformation, not the
material azimuthvartheta. An arbitrary beta gauge does not create a physical
internal rotation period. If one later reconstructs material angles, their
exact equation isvartheta_t=xi/r^2 along the actual material trajectory, and
the momentum/endpoint action in(19) must be restored.

This distinction is relevant to the proposed quantum bridge: a phase in a
canonical chart, a material rotation and a nonconstant compact orbit of Euler
fields are different objects. Their action integrals agree only through the
explicit reduction/reconstruction above. No action quantum or compact
state-space spin orbit follows by naming one of these coordinates a phase.

## Boundary of the canonical statement and the next restoring operator

Equations (1)--(10) are an exact regular-strip canonical chart with explicit
end data.  The physical Euler operator (14) is global through the Hodge
recovery.  On a flat exterior swirl label, `beta` has additional gauge freedom
because `{xi,beta}=0` for every `beta`; a single nondegenerate global canonical
chart is not asserted.  The material Routh reduction is global on the declared
finite-excess domain only in its relative form(19a), with the relative volume
constraint retained. A global first-order chart action for arbitrary nearby
states still requires compatible strip charts and their end/exterior gauge
reduction.

The immediate restoring problem is now precisely typed.  On a fixed mixed-
Casimir leaf, evolve (14) with `delta J=0`, quotient the translation zero mode
(15), and retain the nonlocal Bessel exterior.  The branch-speed derivative
(16) lies outside that fixed leaf at leading order because (13) varies with
the branch.  It is therefore a modulation parameter with a Casimir companion,
rather than an unqualified generalized eigenvector.  The column limit has the
positive propagating metric and strict phase/group-speed ceiling established
in 0030; the localized coefficients are small but couple axial frequencies at
the zero-frequency threshold.  Controlling that threshold coupling is the
named next construction.  Neither the exact canonical chart nor its six
symbolic tests establishes nonlinear persistence by itself.
