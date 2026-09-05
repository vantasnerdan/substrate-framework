# Velocity-interface representation repair before freeze

The D construction inherited0179's divergence-compatible lift
curl_K(T/lambda), whose first correction is q=kappa cross T/lambda.
The initial velocity-interface derivation instead used the first
correction of the DIFFERENT Kelvin lift P_K(D cross omega). This made
its b_V row incorrect even though the standalone moment algebra passed.

The actual Leray derivative is
L1 V=-grad Delta^(-1)(kappa.T)-P[(kappa.u)V]. Since div q=-kappa.T,
the corrected row is b_V=-P[(kappa.u)V+q]. The new verifier directly
assembles the Fourier pressure derivative and checks this identity.

For the exposing planar input the corrected stream row is d alpha,
so every Jordan moment applies to phi_D-d alpha rather than phi_D.
This shifts the two solved first-shell controls by eta -> eta+1.
The psi^5 and psi^7 residual polynomials and the no-Jordan verdict for
the cubic family remain unchanged. Earlier stdout files retain the
original row; the final execution checks the corrected physical row.

The full stationary Euler D fields, direct initial current, complete
D energy, correlated helicity projection, and positive degree-seven
construction never used b_V and are unchanged. This repair does not
license a common-V or physical tag closure; those remain actual equations.
