# Actual velocity correctors and the observed-current quotient

## Frozen routes and first-cell conventions

This is the registered0179 construction on v0.178.0. It changes actual
microscopic initial data, not Euler's equation, the measured mean, or the
parent objective. The stationary velocity normal form and the weaker
observed-current quotient are distinct candidates. A failure of the former
is not a new acceptance condition for the latter.

For curl u=lambda u and a solenoidal velocity w the exact generator is

    L w=P[u cross(curl-lambda)w].                         (1)

For a translation phase T_D=-D.grad u, use the divergence-compatible lift

    T_D(K)=curl_K(T_D/lambda),
    w_D=T_D+ik[q_D+z], q_D=kappa cross T_D/lambda,
    div z=0, mean z=0.                                   (2)

The actual first-cell equation is

    z_t=Lz+F_D,
    F_D=-P[u cross(kappa.grad T_D)]/lambda.               (3)

It follows by expanding the exact full Bloch generator, not differentiating
an isolated projector bound. The measured mean acceleration is

    m_t=k^2 P_kappa <(kappa.u)(q_D+z)
                          +u[kappa.(q_D+z)]>+O(k^3).     (4)

All terms are the physical Eulerian mean stress. A chosen corrected velocity
has its own exact initial circulation class; the corresponding initial
material displacement remains D, and Lin's formula fixes its initial time
derivative. The mean phase pairing with an independent common-V phase
remains rho at zeroth order. Higher phase/current terms must be calculated,
not identified with that leading mass.

For the specific preparation considered here the exact initial phase can
in fact be kept unchanged, not merely its leading mass. Set eta_D(0)=D,
eta_V(0)=0, velocity_V(0)=V, and complete the microscopic D-velocity as
T_D(K)+ik P_K z0. It has zero mean and is exactly Bloch solenoidal.
Lin's formula supplies eta_D,t(0). Its canonical momentum has only
microscopic nonzero modes besides the common V momentum, since mean u=0.
Consequently the complete initial macro pairing is rho J. In the0180
notation this is the special choice chi0=0: neither its physical-mean slip
nor its internal initial one-form is silently assigned zero for arbitrary
corrected tags. Each changed microscopic velocity has its own actual
initial circulation class.

## 1. Exact stationary-cell obstruction for the two-wave field

Shift y by pi/2 and use planar coordinates Y,Z. The actual field is

    u=(cosY+cosZ, sinZ, -sinY), lambda=-1,
    psi=cosY+cosZ, alpha=cosY-cosZ,
    A=sinZ partial_Y-sinY partial_Z.

For a planar solenoidal correction with streamfunction phi, its x-vorticity
is Delta phi. The exact planar Euler vorticity generator is

    curl_x(L z)=A B phi, B=-Delta-1.                      (5)

Let d=D_Y*kappa_Y-D_Z*kappa_Z. The forcing(3) gives
curl_x F_D=-d A alpha. Consequently a stationary first velocity corrector
would require

    A(B phi-d alpha)=0.                                  (6)

Every regular level psi=c on the square torus, c in(-2,0) or(0,2), is
one connected closed streamline. Thus an L2 first integral is f(psi)
off the measure-zero separatrices. The exchange Y<->Z preserves psi and
reverses alpha, so <f(psi)alpha>=0. But B alpha=0 and B is self-adjoint.
Taking the alpha pairing of B phi=d alpha+f(psi) therefore gives

    0=d <alpha^2>, <alpha^2>=1.                           (7)

For d!=0 this is impossible. This is a full transport/Fredholm obstruction,
not absence of a solution in one Fourier shell. Any three-dimensional
periodic corrector would imply a planar one by averaging in the invariant
x direction; that average commutes with the planar Euler equation. The
obstruction thus includes that larger smooth periodic class. For example
kappa proportional to(0,1,1), D proportional to(0,1,-1) is physically
transverse and has d!=0.

This does not refute constant observed mean stress. It refutes only the
stronger stationary-microscopic-velocity route for these generic strains.

## 2. The solvable elementary-wave stationary route selects the wrong sign

For u=(cos z,sin z,0), lambda=-1 and n=e_z, the forcing(3) is exactly0.
Thus z=0 gives an actual stationary first-cell velocity. Put c=<|u|^2>/2.
Direct evaluation of(4) yields

    R_D=c(D.n) P_kappa n.                                (8)

It has the inverted sign: m_t=+k^2 R_D. Whole-field Haar averaging gives
+c D/3 on the physical transverse bundle. The exact bare Kelvin-D/common-V
preparation from0151 instead gives the positive restoring response
-2<|u|^2>D/15. Hence imposing microscopic stationarity can destroy the
desired physical branch even in a geometry where the equation is solvable.
This is a representation/selection failure, not an Euler no-go.

## 3. The weaker actual-current criterion

For fixed physical kappa,D, denote the linear part of(4) by C z. Actual
constant observed stress requires

    C exp(t L)(L z0+F_D)=0                               (9)

on the specified time interval. It does not require Lz0+F_D=0. Equivalently
the residual belongs to the unobserved invariant subspace, defined by the
actual semigroup and physical current row. This is a useful exact criterion
only when that subspace or an adjoint current closure is constructed; it
is not a licence to assume a favorable solution.

For an exact finite Fourier diagnostic, freeze the zero-mean solenoidal
first-shell initial class and the physically transverse direction
kappa=(0,1,1), D=(0,1,-1), with normalization restored in any coefficient.
Use the full, untruncated Euler generator for each successive derivative.
The exposing rows are

    C L^j(L z0+F_D), j=0,...,8.                           (10)

The linear unknowns are only the initial first-shell coefficients; all
generated modes remain. Compare the coefficient and augmented ranks as
rows are added. A finite row mismatch refutes this INITIAL class for exact
current closure; a consistent finite fit remains only a finite fit, not
an all-time construction. This diagnostic is exact algebra, with no solver
frequency, empirical comparator, tolerance or numerical spectral design.
The next candidate after this test is an actual adjoint-current graph
corrector, with the successful generalized0163 construction as the explicit
different-background comparison. No same-EPS cohomology transfer is assumed.
