# Exact physical column waves and the solitary speed boundary

This is a direct continuation of0030's Euler action calculation on the same
column as0027. It is a full-pressure axisymmetric linear result about that
background, not a stability theorem about the nonzero solitary wave. Use
L=r U_theta, Phi=2LL'/r^3, K=(-Delta_*/r^2)^-1 and dnu=r dr dz. The physical
energy multiplier remains2pi rho_m.

## Full-pressure restoring operator

Writing eta=delta(omega_theta/r), chi=delta(r u_theta), the actual Euler
linearization is

    eta_t=2L/r^4 chi_z,  chi_t=L'/r (K eta)_z.             (1)

No pressure term is discarded: K reconstructs the complete poloidal velocity
on the whole radial half-line. The positive metric in action-construction.md
makes(1) skew with respect to

    ||(eta,chi)||_C^2=int[eta K eta+2L/(r^3 L') chi^2]dnu. (2)

Its closure is taken from smooth accessible perturbations with finite(2).
The flat exterior carries no independent accessible swirl perturbation because
L'=0 there; its poloidal kinetic energy remains in K. Multiplying(1) by the
metric, integrating in z, and then completing in(2) gives a unitary group on
this column energy space. Equivalently each nonzero Fourier sector is the
skew-adjoint oscillator constructed from the positive self-adjoint radial
operator below. This construction provides existence and conservation, rather
than assuming existence merely from the formal energy identity.

For a mode exp(i k z-i omega t), psi=K eta satisfies

    c_ph^2[-psi''+psi'/r+k^2 psi]=Phi psi,
    c_ph=omega/k.                                      (3)

On r>R, beyond the compact vorticity support, its exact finite-energy solution
is the0027 Bessel exterior. Eliminating it gives the generalized radial form
on V={f(0)=0, int_0^R f'^2/r finite}:

    A_k[f]=int_0^R(f'^2+k^2 f^2)/r dr+T_R(k)|f(R)|^2,
    B[f]=int_0^R Phi f^2/r dr,
    T_R(k)=|k| K0(|k|R)/(R K1(|k|R)).                   (4)

The compact positive operator A_k^-1/2 B A_k^-1/2 supplies the nonzero
oscillator frequencies. Its largest eigenvalue is c_1(k)^2. The radial
Sturm equation makes each positive radial eigenvalue simple; modes are
indexed in decreasing speed. A kernel, if present in the completed ambient
space, gives stationary rather than oscillating directions. None is removed
by assigning it a positive frequency.

## The solitary branch is strictly faster than every linear axisymmetric wave

The variational threshold in0027 is exactly

    c_0^2=sup B[f]/A_0[f].                              (5)

For every k!=0, A_k>A_0 on a nonzero mode. Thus

    0<c_n(k)<=c_1(k)<c_0.                               (6)

This is not a shallow-water replacement: the complete exterior symbol occurs
in(4). The limit c_1(k)->c_0 follows from form convergence as k->0. The
threshold is a resonance of the whole exterior problem, not an L2 radial
eigenfunction at k=0.

There is also an exact group-velocity bound. The exterior variational problem
is the minimum of int_R^infinity(f'^2+k^2 f^2)/r at fixed unit trace. Its
envelope derivative gives

    0<=k T_R'(k)<=2T_R(k), k>0.                         (7)

Consequently 0<=k A_k'[f]<=2 A_k[f]. Differentiating a simple generalized
eigenpair A_k f_n=c_n^-2 B f_n gives

    d(k c_n)/dk
       =c_n[1-k A_k'[f_n]/(2 A_k[f_n])].               (8)

Hence every positive-frequency branch has group velocity between0 andc_n,
and every negative branch has the opposite velocity. In particular all
axisymmetric group speeds have magnitude at mostc_0. The actual solitary
speed c>c_0 outruns them. In its translating frame,

    |c k +/- omega_n(k)| >= (c-c_0)|k|.                (9)

Equation(9) is a frequency-dependent nonresonance estimate. It loses its gap
at k=0 and therefore does not furnish a uniform inverse on the continuous
spectrum, nor does it settle the nonlinear solitary stability problem.

## Explicit next use

The exact positive background metric and the supercritical speed separation
supply the natural incoming/outgoing variables for a full solitary-wave
modulation calculation. The next construction retains the k=0 long-wave
channel, its translation mode and speed companion, the frequency-dependent
Bessel pressure kernel, and the other radial channels. An actual Euler
wave-to-amplitude map and an estimate for the retained complement are needed
before the scalar homoclinic action can control disturbances of the wave.

## Oscillator domain rather than formal conservation alone

For fixed k!=0, complete core pairs using(2). Put p=K_k^(1/2) eta and
q=a^(1/2) chi, wherea=2L/(r^3 L') and K_k is the full exterior inverse
restricted back to the core. Equations(1) become

    p_t=i k D q, q_t=i k D^* p,
    D=K_k^(1/2) (2L/r^4) a^(-1/2).                  (10)

Indeed a^(1/2)L'/r=(2L/r^4)a^(-1/2); all multiplications are real.
The squared singular values ofD are the generalized eigenvalues in(3)-(4),
so||D||=c_1(k)<=c_0. The real symmetric block[[0,D],[D^*,0]] is bounded
self-adjoint for this fixedk. Its exponential defines a unitary group in(p,q),
and the measurable direct integral overk defines the full z-evolution on
the weighted energy completion. Its generator domain is the natural direct
integral domain, not an asserted bounded generator across allk. Smooth
accessible data form the physical core; density extends the exact group.
