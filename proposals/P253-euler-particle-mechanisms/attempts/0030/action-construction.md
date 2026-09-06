# Actual Euler action for the solitary branch and the stable-column comparison

This calculation is conditional on0027's actual existence hypotheses, under
independent0028 review. It does not transfer the geometric NLS Hamiltonian.
Use axisymmetric variables xi=r u_theta, zeta=omega_theta/r, dnu=r dr dz,
psi_lab=K zeta, and the physical common multiplier2pi rho_m. Normalized energy
and axial impulse are

    E=1/2 int(zeta K zeta+xi^2/r^2)dnu,
    I=1/2 int r^2 zeta dnu.                           (1)

The convergent relative integrals are understood about the fixed column.
The exact axisymmetric Casimirs have forms int C(xi)dnu and int zeta D(xi)dnu.
On a regular label interval where the traveling branch xi=F(psi_frame) can be
inverted, write Psi=F^-1. Its steady relation is

    psi_frame=psi_lab-c r^2/2=Psi(xi),
    zeta=F F'/r^2-B'.                                (2)

## Exact translating critical functional

The conserved Euler energy-Casimir functional is

    A=E-c I-int zeta Psi(xi)dnu-int B(Psi(xi))dnu.      (3)

Its first variations are

    A_zeta=psi_frame-Psi(xi),
    A_xi=xi/r^2-zeta Psi'(xi)-(B composed Psi)'(xi).    (4)

Using(2), both vanish on the actual wave in the regular interior. This is a
physical full-field statement; outside the invertible label interval one must
retain the corresponding boundary/variational inequality rather than assign
an inverse to a constant exterior label.

The quadratic form on q=(eta,chi) is

    Q=int[eta K eta-2Psi' eta chi
       +(1/r^2-zeta Psi''-(B composed Psi)'')chi^2]dnu. (5)

Physical Q is2pi rho_m times(5). The exact accessible axisymmetric tangent is

    chi={xi,g}, eta={zeta,g}+{xi,a},
    {f,g}=(f_r g_z-f_z g_r)/r.                       (6)

Arbitrary eta,chi are not automatically independent physical directions.
However where xi_r!=0, high-z-frequency g and a vary their leading ratio
independently. Since K has order minus two, the off-diagonal term in(5) has
both unbounded signs in the kinetic norm eta K eta+chi^2/r^2. This rules out
using this particular full-leaf traveling functional as a uniformly positive
kinetic metric. It does not prove instability of the wave.

## A positive background functional is available in a different frame

The reference column is stationary in the laboratory frame. Let L(r)=xi_bg
be strictly increasing where w>0 and let Omega=L/r^2 be nonincreasing. This
last condition is satisfied by a smooth nonincreasing w(r^2), constant near0
and flat at its support edge R_0. Define on0<xi<L_infinity

    C'(xi)=-Omega(r(xi)),
    C''(xi)=-Omega'(r(xi))/L'(r(xi))>=0.              (7)

C has its continuous convex endpoint extension at L_infinity. On the core,
E+int C(xi) is critical at the column; its quadratic coefficient is

    1/r^2+C''(L(r))=2L/(r^3 L')>0.                  (8)

The linearized exact column equations are

    eta_t=(2L/r^4) partial_z chi,
    chi_t=(L'/r) partial_z K eta.                    (9)

Hence the positive quadratic form

    1/2 int[eta K eta+2L chi^2/(r^3 L')]dnu          (10)

is exactly conserved by integration by parts in z. It controls the kinetic
norm when w is nonincreasing because2L/(r L')>=1. Near the flat edge the
weight can diverge; exact accessible chi contains L', and the domain of(10)
must retain that weight. It is not equivalent to an unrestricted unweighted
kinetic topology at the edge.

There is also a nonlinear conditional energy inequality. For axisymmetric
states with0<=xi<=L_infinity and convergent relative energy/Casimir integrals,
let chi=xi-L and psi=K zeta. On r<R_0 the linear term cancels exactly; convexity
of C gives its nonnegative Taylor remainder. On r>=R_0, L=L_infinity and

    [L/r^2+C'(L_infinity)]chi
       =[L/r^2-L/R_0^2]chi>=0,

because chi<=0. Thus the conserved relative functional, when its transport
and energy fluxes vanish at infinity, bounds

    Delta(E+int C)>=1/2 int[zeta K zeta+chi^2/r^2]dnu. (11)

This gives a physical kinetic Lyapunov bound about the background column for
as long as the classical axisymmetric evolution and the declared integral
conditions hold. It neither proves global regularity with swirl nor stability
around the nonzero solitary excitation. Those are different positive tasks.

## The traveling impulse is an actual momentum plus a Casimir

For the pure-swirl background of0027, Psi(xi)=-c r(xi)^2/2. Moreover
(B composed Psi)'=Omega(r(xi)), independent of c up to a constant.
Therefore(3) can be written

    A=H_C-c P_C,
    H_C=E+int C(xi)dnu,
    P_C=I-(1/2)int zeta r(xi)^2 dnu.                 (12)

P_C differs from the translation momentum only by a Casimir. It generates the
same spatial translation on the physical leaf. At the column its leading
nonzero term is

    P_C=-int [r/L'] eta chi dnu+higher terms.          (13)

The positive background energy and the indefinite traveling quadratic form
are therefore compatible: subtracting c times the actual momentum changes
its second variation. The conserved scalar steady functional in0027 is not
substituted for this physical Hamiltonian by notation.

## Next exact construction

To resolve the solitary wave, derive its complete Hamiltonian generator and
continuum scattering or positive propagation norm using(5)-(6),(12). The
supercritical axial speed and actual nonlocal exterior are possible tools;
the absence of a positive full-leaf Hessian is not a terminal result. First
retain all axisymmetric continuous-spectrum modes and the phase companion;
then inspect the nonaxisymmetric azimuthal sectors. A scalar homoclinic
stability theorem only becomes relevant after its actual Euler dynamical and
physical energy/current map is proved.
