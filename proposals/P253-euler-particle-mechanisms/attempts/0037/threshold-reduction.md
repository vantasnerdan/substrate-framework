# The Euler threshold forces a KdV leading system

This derivation uses only the exact `0027` branch and the exact `0030/0034`
Euler equations. It identifies the complete leading scalar threshold problem;
it does not replace the full Euler complement by that scalar equation.

## 1. Long-wave dispersion fixes the linear coefficient

Let `mu=lambda(c)>0`, where `lambda(c0)=0` is the simple radial eigenvalue
used in `0027`, and put

    sigma=(lambda'(c0))^-1>0,
    lambda'(c0)=2 int Phi f0^2 dr/r / c0^3.          (1)

For the first right-moving column branch, the exact generalized radial
eigenproblem and exterior Dirichlet-to-Neumann form give

    lambda(c_ph(k))+k^2+T_R(k) f0(R)^2+o(k^2 log(1/|k|))=0,
    T_R(k)=k^2[log(2/(|k|R))-EulerGamma]+o(k^2).     (2)

The complement and eigenfunction corrections in `0027` are lower order in
this scalar row. Hence

    c_ph(k)-c0
      =-sigma f0(R)^2 k^2 log(1/|k|)+O(k^2).         (3)

Choose the exact branch scale `L_mu` by

    mu L_mu^2=f0(R)^2 log L_mu.                      (4)

With

    X=(z-c0 t)/L_mu,  T=mu t/L_mu,  k=kappa/L_mu,

equation (3) gives the scaled frequency

    Omega_mu(kappa)
      =(omega(k)-c0 k)L_mu/mu ->-sigma kappa^3       (5)

uniformly on compact nonzero `kappa` sets, and continuously at zero. The
ratio `log(L_mu/|kappa|)/log L_mu` tends to one. Thus the limiting linear
amplitude evolution is `A_T+sigma A_XXX=0`. The logarithmic Bessel exterior,
rather than a local wall condition, supplies this coefficient.

## 2. Euler quadratic projection fixes the nonlinear coefficient

The exact axisymmetric perturbation equations about the column are

    chi_t=(L'/r) psi_z-{psi,chi},
    eta_t=(2L/r^4)chi_z+(2/r^4)chi chi_z-{psi,eta},
    psi=K eta.                                       (6)

The right-moving critical vector has

    psi=mu f0 A+...,
    chi=-mu L'f0 A/(c0 r)+... .                      (7)

Projection of the quadratic terms in (6) onto the normalized radial kernel
gives the same coefficient already exposed by the exact steady reduction,

    beta=(1/2) int J_c0(r) f0(r)^3 dr/r>0.           (8)

Translation invariance and the conservative first-order Euler system put the
scalar quadratic term in divergence form. Its normalization is fixed without
fitting by the exact branch: a profile moving at scaled speed `sigma` must
satisfy the already proved equation

    A-A_XX=beta A^2.                                 (9)

Combining (5), (8), and (9) fixes the only compatible leading scalar equation
as

    A_T+sigma partial_X(A_XX+beta A^2)=0.            (10)

Indeed a wave `A(X-sigma T)` reduces (10) exactly to (9). This coefficient
match is an Euler solvability identity. The claim that the full Euler flow is
approximated by (10) on `T=O(1)` still requires the complement estimate below.

## 3. Exact scalar restoring structure

The positive solitary profile is

    A_*(X)=3/(2 beta) sech^2(X/2).                   (11)

Equation (10) is Hamiltonian with Poisson operator `sigma partial_X` and

    E_KdV=integral[(A_X^2)/2-(beta A^3)/3]dX,
    P_KdV=(1/2)integral A^2 dX.                      (12)

The profile is critical for `E_KdV+P_KdV`. Its Hessian is

    L_*=-partial_XX+1-3 sech^2(X/2).                 (13)

It has one simple negative eigenvalue `-5/4`, the simple translation zero
mode `A_*'`, one further discrete positive eigenvalue `3/4`, and essential
spectrum `[1,infinity)`. For the speed-`v` family

    A_v=3v/(2 sigma beta)
        sech^2(sqrt(v/sigma)X/2),                    (14)

the exact scalar invariants are

    integral A_v dX=6 sqrt(v/sigma)/beta,
    P_KdV(A_v)=3(v/sigma)^(3/2)/beta^2.              (15)

Both slopes are positive. The standard one-negative-direction constrained
calculus therefore gives scalar orbital coercivity modulo translation. Here
that statement is proved only for (10); it is not imported as Euler
stability.

The two constraints have direct Euler candidates. The first variation of the
`0034` end distribution is proportional to the first integral in (15), while
the exact translation momentum/Casimir `P_C` supplies the second. A full
Euler proof must calculate those pullbacks through the next order and show
that their constraint matrix stays nonsingular.

## 4. The named full-Euler remainder

Write the solution as its critical component plus radial complement and use
the exact `0027` rank-one exterior inverse. The missing theorem is a uniform
estimate, on physical times `|t|<=C L_mu/mu`, for all of:

1. the time-dependent complement Duhamel operator, including the derivative
   of the Bessel symbol at `k=0`;
2. the nonlinear chart/end-distribution rows enforcing `delta J=0`;
3. the translation modulation and exact momentum constraint;
4. quadratic and cubic terms with one complement factor; and
5. reconstruction in the finite-excess velocity norm.

The static multiplier convergence used for existence controls no growing
time interval by itself. A bound `O(1/L_mu)` at fixed time can accumulate to
`O(1/mu)` on the modulation horizon. Route B therefore establishes the
forced leading system and its exact restoring structure, while the Euler
transfer remains active at this single operator-valued long-time remainder.

Route A attacks that remainder through the constrained resolvent. Route C
attacks it through an Evans/Birman--Schwinger construction. Neither route may
drop the `k=0` exterior, replace the mixed-Casimir row by parity, or infer
three-dimensional stability from this axisymmetric calculation.
