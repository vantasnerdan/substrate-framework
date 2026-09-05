# Actual stationary-tag fourth harmonic in the elliptic affine core

This is a constructive affine-core result and its explicit nonlinear-cell
transfer interface, not a completed fixed-cell second spatial jet.
Its physical observation differs materially from0155's finite-X marker.

## 1. Actual material-tag equations before a mode is selected

In the physical normal plane (a,b) use the exact cell coordinates from
0205, with axial X. Put u_perp=J grad psi, J=[[0,-1],[1,0]], T=u.grad,
H=-Delta_perp. For an X-independent Euler velocity
w_perp=J grad phi, the full horizontal vorticity equation on the
NONLINEAR cell is

    H phi_t=-T(H-1)phi.                              (1)

This includes the actual pressure through velocity recovery; it is not
an arbitrary transport oscillator. The longitudinal velocity is retained
as a separately forced Euler row and does not enter this closed planar
vorticity equation.

Take a literal stationary nonnegative material fraction chi(psi), with
compact transverse support inside the elliptic core and the FULL axial
circle. Its physical unperturbed centroid and covariance are stationary.
Set Q=integral chi(a²-b²)dA, and assume the actual Q is nonzero. The
linear centered covariance angle is

    theta=integral delta_chi a b dA / Q.              (2)

Centroid variations do not contribute linearly since both reference
first normal moments vanish. For a solenoidal horizontal displacement
xi=J grad s, delta_chi=chi_psi T s. Its exact transported density obeys

    delta_chi_t+T delta_chi=chi_psi T phi.            (3)

The actual mechanical axial spin is

    S=rho integral [chi(a w_b-b w_a)
                         +delta_chi(a u_b-b u_a)]dA. (4)

No angular inertia is inserted in(2)–(4). All quantities are per actual
axial length; divide by the physical transverse cell area for densities.

For this stationary tag, an especially useful connection vanishes by
an exact integral, not by renaming momentum:

    integral chi (xi cross u)_X = -integral chi T s=0,
    S=G_t, G=rho integral chi(a xi_b-b xi_a)dA.        (5)

The transport is divergence free and T chi=0. This is why the stationary
full-X tag has a different physical current license from painted finite-X
lobes. G still has its literal initial value.

## 2. Exact affine Euler field and full anisotropic Poisson inversion

The leading local physical field is

    u_aff=(W0,b,-Omega² a), p_aff=Omega²(a²+b²)/2,
    Omega=1/10, W0=101/100.

It is actual stationary Euler on the covering affine domain. It is not
claimed to be a finite-energy periodic background. Define area-preserving
normal coordinates (R1,R2)=(sqrt(Omega)a,b/sqrt(Omega))/ell,
r²=R1²+R2², Gmetric=diag(Omega,1/Omega),

    T=-Omega partial_theta,
    Delta_phys=ell^(-2) Delta_Gmetric.

The affine normal vorticity is constant, so the ACTUAL linear vorticity
perturbation obeys zeta_t+T zeta=0. Select

    zeta=zeta0 r^4 exp(-r²/2) exp(i4theta+i4Omega t).  (6)

Its velocity is the full decaying plane Biot-Savart field, recovered from
Delta_phys phi=zeta. The physical pressure then follows from the full
linear Euler equations. Neither a radial wall nor a local inverse is used.

Let tau=(1-Omega)/(1+Omega), so0<tau<1. In Fourier polar coordinates
the exact inverse metric symbol has the convergent Poisson expansion

    1/[Omega cos²vartheta+Omega^(-1)sin²vartheta]
       =sum_(j in Z) tau^|j| exp(i2j vartheta).       (7)

The unitary plane Fourier transform of the Gaussian fourth harmonic in
(6) is k^4 exp(-k²/2) exp(i4vartheta). Inverting the FULL symbol gives
these three exact angular coefficients of phi/(zeta0 ell²):

    phi_2=tau r² exp(-r²/2),
    phi_-2=tau³ r² exp(-r²/2),
    phi_0=-tau²(2-r²)exp(-r²/2).                    (8)

Other even angular sectors remain in the actual velocity and action.
Equations(8) are projections of the complete Poisson solution, not a
three-mode approximation to it. The defining Hankel integrals are
integral k³ exp(-k²/2)J2(kr)dk=r² exp(-r²/2) and
the corresponding J0 integral=(2-r²)exp(-r²/2).

## 3. Select actual Lin data, retaining the unobserved resonance

For every l!=4 choose the initial stream displacement

    s_l=phi_l/[i(4-l)Omega].                          (9)

Then its actual Lin history is that same coefficient times exp(i4Omega t),
since s_t+T s=phi. The l=4 sector is resonant. Retain it explicitly:

    s_4(t)=exp(i4Omega t)[s_4(0)+t phi_4].           (10)

There is no assumed full-field monochromatic mode or discarded secular
state. The stationary covariance observation sees only l=±2, and its
spin/G rows see l=0. Thus(10) is unobserved by these particular physical
rows in this exact affine problem. This distinction, not isolated
spectral autonomy, licenses the observed harmonic.

Put z=T s, so delta_chi=chi_psi z. Equations(9) give

    z_2=-phi_2, z_-2=phi_-2/3, z_0=0.               (11)

All full Fourier sectors stay present in the actual preparation.

## 4. Positive measured spin, physical angle and initial current

Choose a smooth nonnegative chi(x), x=r², supported in0<=x<4, with
nonzero mass; more general fixed profiles with the same nonzero response
integral are allowed. This is chi of the affine streamfunction. Define

    I_chi=integral chi(x)x(4-x)exp(-x/2)dx>0,
    Q=pi ell^4(Omega^(-1)-Omega)/2 integral chi(x)x dx>0.

Use the exact physical Jacobian and (2),(4), not a calibrated tensor
angle. The complex angle and spin amplitudes are

    theta=-i C exp(i4Omega t),
    C=pi zeta0 ell^4(tau+tau³/3)I_chi/(4Omega Q),
    S=pi rho zeta0 ell^4 tau² I_chi exp(i4Omega t).   (12)

For example, the velocity contribution to spin is obtained by integrating
chi r.grad phi by parts; the resulting radial functional is
-2pi rho integral(chi+x chi_x)phi_0 dx. The moving-tag contribution
selects z_0 and is zero by(11), not omitted. The angle selects z_2-z_-2.
These distinct integrals yield

    S=Delta theta_t, G=Delta theta,
    Delta=rho Q tau/(1+tau²/3)>0.                   (13)

This is the SAME stationary physical tag and actual total mechanical
spin. It is not a chosen inertia or a phase winding. Equation(5) confirms
the complete current relation, including its initial G. The two real
phases of(12) have nonzero angle/rate determinant4Omega C².

The tag's mass is O(rho ell²) per axial length and Delta is O(rho ell^4).
On the fixed whole cell j=Delta/(transverse cell area). Hence every
finite ell gives nonzero actual j, but j tends to zero as ell shrinks.
No nonzero-density limiting inertia follows from a small core alone.

## 5. Actual phase normalization without changing those observations

Use the full initial Lin phase, with projected cotangent
pi/rho=w+Du_aff xi. For its two real phases,

    beta=rho Im integral xi*. (w+Du_aff xi)dA.

The constant-strain term has zero imaginary integral: its antisymmetric
part is proportional to the integrated Jacobian of the two decaying
stream displacements, a boundary integral. The other term is

    integral grad s*.grad phi=-integral s* zeta.

Only angular l=4 contributes. Choose the ACTUAL initial unobserved
resonant row s_4(0)=i c zeta_4, c real. It follows that

    beta=rho c ||zeta_4||²,
    ||zeta_4||²=24pi zeta0² ell².                    (14)

Thus c=Delta(4Omega)C²/[rho||zeta_4||²]>0 makes the actual physical
angle/rate chart mass precisely Delta. In the Re/Im initial-column
convention its Wronskian is+4Omega C²; the sign is derived in that
convention rather than imported from0155's different encoding. This
extra s_4 changes none of theta,S,G in(12),(13), but its complete history
(10), phase and energy remain. It scales as ell²/Omega, the same physical
order as the nonresonant displacement; no inverse shrinking-core norm
is hidden in this normalization.

The actual observed scalar action is consequently

    L_theta=Delta[theta_t²-(4Omega)²theta²]/2,         (15)

on this exact two-phase prepared affine Euler history. The full conserved
Jacobi energy is a separate complete quadratic form;0205's actual
energy returns are an available repair if it differs from(15)'s mechanical
energy. Neither positivity of(15) nor phase matching silently changes it.

## 6. Fixed-cell continuation still to derive

The affine source is a concrete actual finite-action Euler perturbation,
not yet a theorem about the nonlinear periodic cell. Its full velocity
and Lin tails, the nonlinear transverse vorticity gradient, and the
longitudinal Euler response must be controlled on the same cell before
transferring(12)–(15). A finite-core approximation must compare angle,
spin and phase errors to C and Delta, not only to a bulk velocity norm.

The common-K optical second jet, its sign, the nonlinear-cell clock and
the positive-density/compact-EPS joining remain actual constructions.
In particular the small-core leading equation alone cannot replace a
persistent clock error by an error smaller than j times a small spatial
coefficient. The registered fixed-cell observed-response and whole-family
phase-correlation competitors remain active routes for that step.
