# A finite-Fourier stationary Euler corrector on unequal circular waves

This is the failure-derived method repair inside the registered unequal-
wave family. It does not extrapolate the equal-wave0184 obstruction.

Let A,B>0, A!=B, and use the actual stationary Euler field

    u=(psi,A sin Z,-B sin Y),
    psi=B cos Y+A cos Z, alpha=B cos Y-A cos Z,
    T=A sin Z partial_Y-B sin Y partial_Z,
    H=-Delta, mathcalB=H-1.

It is constant-curl, curl u=-u, with p=-|u|²/2. The actual planar
first-cell velocity equation, for a translation strain coefficient d,
is H phi_t=-T mathcalB phi+d Talpha. No Euler mode or pressure is
truncated to obtain this equation.

## 1. The first-integral moment matrix is nonsingular when A!=B

Seek a stationary solution with

    mathcalB phi=d alpha+c1 psi+c3 psi³.                (1)

Because Tpsi=0, this solves the exact transport equation if its first
Fourier shell vanishes. The two nontrivial kernel moments are

    <cosY psi>=B/2, <cosZ psi>=A/2,
    <cosY psi³>=3B(B²+2A²)/8,
    <cosZ psi³>=3A(A²+2B²)/8.

Their determinant is 3AB(B²-A²)/16, not zero for the selected field.
Solving the actual two equations gives

    c1=-3d(A²+B²)/(B²-A²), c3=8d/[3(B²-A²)].           (2)

All constant and sine kernel moments vanish by parity. Therefore(1)
has the explicit finite-Fourier solution

    phi=d{[B³ cos(3Y)+A³ cos(3Z)]/[12(B²-A²)]
          +[B²A cos(2Y)cosZ+BA² cosY cos(2Z)]/(B²-A²)}. (3)

There is no small-divisor inverse or fitted spectrum: the remaining modes
have mathcalB eigenvalues8 and4. The singular A->B limit is real and
explains why the equal-amplitude parity obstruction was not contradicted.

## 2. Complete velocity, not just planar vorticity

Set r=-d alpha/2 and b=r+Hphi. Then

    z=(b,-phi_Z,phi_Y)

is a smooth mean-zero solenoidal field and satisfies the FULL stationary
Euler corrector equation Lz+dF0=0, with F0 the actual cross-wave force.
Indeed Tr=-(d/2)Talpha is its required axial-plus-vorticity equation,
and(1) supplies its planar curl. The pressure is the full Leray pressure;
curl reduction did not omit an axial equation or a harmonic mean.

One may add first-shell phi components while preserving the planar
stationary equation. Choosing b=r+Hphi adds their actual Beltrami velocity
tangents, not a scalar output adjustment. For example phi -> phi+d eta
alpha changes the WHOLE physical static current by

    delta R_D=-eta(A²+B²)/10.                           (4)

The particular(3) has no first-shell projection, so its whole current
correction is zero; the fixed translation lift alone gives
(A²+B²)/6. Thus a first-shell physical preparation can change this static
D stress, but its energy and the common-V column have not been supplied
by the sign of(4).

## 3. Full energy and the pending actual velocity partner

For eta_D(0)=D, w_D=T_D+ik(q_D+z), D transverse, the actual initial
material rate is ik(q_D-aD+z), where a=kappa.u. Its leading divergence
vanishes. The complete inherited Hamiltonian coefficient is

    H_DD^(2)=rho/2{||q_D-aD+z||²-||aD||²}.             (5)

This follows from the constrained Legendre transform, including the
normal pressure projection. Its quadratic dependence on the first-shell
freedom cannot be replaced by the linear observed current(4).

The other required column begins with the genuine mean velocity V.
At zero K, Euler enforces w_V(t)=V-t(Du)V when no extra leading
microscopic velocity is introduced. At first K an actual Jordan/
configuration-return equation must connect this column to the chosen
D preparation. The full initial phase and energy, physical mean and
stationary tag angle/spin rows are the tests. A stationary D corrector
or favorable static stress alone does not imply any of them.
