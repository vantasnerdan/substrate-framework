# Moving-frame action, complement and physical observation

Use the Euler/KKS sign convention of accepted C-CST-009. The ambient
quadratic action is

    L(x,xdot)=-x^T Omega xdot/2-x^T H x/2,
    A=-Omega^-1 H, xdot=A x.

Omega is constant, skew and nondegenerate; H is symmetric. All finite
matrix identities below extend to the appropriate Euler bilinear forms
only on their actual common domains. The canonical implementation tests
the identities, not domain/existence hypotheses for a PDE.

## Complete time-dependent pullback

For x=E(t)z put Omega_E=E^T Omega E and Q=E^T Omega Edot. Direct
substitution gives

    L_E=-z^T Omega_E zdot/2-z^T H_eff z/2,
    H_eff=E^T H E+sym Q,
    dotOmega_E=Q-Q^T.

The symmetric part of Q alone enters the scalar Hamiltonian. Variation
also differentiates the kinetic one-form; the complete equation is

    Omega_E zdot+(H_eff+dotOmega_E/2)z=0.

Consequently, with Pi=Omega_E^-1 E^T Omega and P=E Pi,

    A_E=Pi(AE-Edot),
    R=Edot+E A_E-AE,
    E^T Omega R=0.

The last identity is symplectic orthogonality, not R=0. In particular
positive H_eff or elliptic A_E supplies no invariance statement by itself.
When a proposed Floquet basis obeys

    Edot=AE-E A_F+r,

the exact identity H=-Omega A cancels all large common advective terms:

    H_eff=-sym(Omega_E A_F)+sym(E^T Omega r).

This is the appropriate starting point for an analytic small-energy bound.
Subtracting two separately large numerical energies would neither prove
the cancellation nor establish a sign below the evaluator floor.

## The complete moving complement

For an unrestricted ambient solution, set z=Pi x and r=(I-P)x. Then

    zdot=A_E z+(dotPi+Pi A)r,
    rdot=-R z+[(I-P)A-dotP]r,
    Pi r=0.

The moving-fiber term -dotP is essential: its equation preserves
dotPi*r+Pi*rdot=0, even though rdot need not lie in the instantaneous
complement. The propagator of this full time-dependent complement gives
its Volterra memory, including the homogeneous evolution of initial r.
An actual observation is

    observable=O(t)E(t)z+O(t)r.

Deleting the last term or setting only r(0)=0 does not close the full
Euler realization. For a finite-time semigroup bound and an actually
measured residual norm, Duhamel gives the corresponding approximation
estimate. Those bounds are application data, not supplied by this algebra.

## Floquet gauge and the physical clock

Under a differentiable invertible coordinate frame F(t), replace E by EF.
The exact transformations are

    Omega_E'=F^T Omega_E F,
    H_eff'=F^T H_eff F+sym(F^T Omega_E Fdot),
    A_E'=F^-1 A_E F-F^-1 Fdot,
    observation'=O E F.

For periodic F the return maps are conjugate, but Floquet logarithms may
differ by integer winding. Thus a coordinate-Hamiltonian sign is not an
invariant of free Euler motion under arbitrary periodic frame choices.

The exposing exact example uses Omega=[[0,1],[-1,0]], H=I and the periodic
physical frame E(t)=exp(2 Omega t). The pulled-back H_eff is -I, while

    Edot+E A_E=A E

holds exactly and its physical observed motion is unchanged. Both energy
signs describe the same underlying oscillator in different frames. A
physical positive comoving rotor therefore specifies its frame, torsion
lift, carrier phase and observation before selecting a sign. Choosing a
different actual lambda/polarization/geometry is a physical candidate;
adding a frame winding just to make H_eff positive is not such a candidate.

## Evidence and next application

The new euler_phase.py derives the matrices from supplied Omega,H,E,Edot.
Its tests independently vary the original time-dependent action with a
nonconstant Omega_E, reconstruct the full moving complement, and evaluate
the periodic winding example. The missing dotOmega_E/2 mutation changes
the defining Euler-Lagrange equation.0114 applies this bookkeeping to the
actual EPS packet forms; its residual, sign and physical moments remain
its own constructive obligations rather than being assumed by this API.
