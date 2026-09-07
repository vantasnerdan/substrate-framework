# Validation and claim boundary

## Exact checks

The attempt verifier derives, rather than numerically samples:

1. the full theta-curl identity (14), including the non-gradient electric
   row;
2. the source-side chain rule `partial_z[g K(P)/r]=g chi(P)P_z/r` when
   `K'=chi`, which is the right-hand row used in the radial Ampere primitive;
3. Gauss plus axial Ampere gives the scalar elliptic equation (17);
4. substitution of the vorticity first integral gives Bernoulli (18);
5. the toroidal-vorticity regular-action stabilizer equations imply tag
   locking only when `zeta(I)` and `zeta'(I)` are nonzero;
6. exact period integration cancels `partial_t omega`, after which the
   degree-minus-six fixed-tail compatibility row contains no faster remainder
   term;
7. the leading circulation/impulse matrix sees radial translation with
   nonzero coefficient for positive circulation and radius;
8. an `O(g)` Maxwell field gives `O(g^2)` Lorentz backreaction.

These checks expose sign, missing-curl, charge-scaling, and false-stabilizer
errors.  They do not prove the exact complement isomorphism HSE or the Schur
determinant (31).

## Analytic evidence

The axisymmetric reduction is a direct consequence of the reviewed joint
action and the exact cylindrical Maxwell and Euler equations.  The affine
potential domain retains the charged `1/r` tail and the scalar Gauss row.
The traveling-potential identity
`E+c e_z cross B=-grad(phi-c A_z)` and the full left-hand radial Ampere
primitive (16) are analytic calculations in the derivation, independently
reviewed as equations; the verifier does not assert either complete identity.
Cao's exact equation supplies the uncharged base and the limiting kernel.  The
source's auxiliary parameter system is explicitly excluded from the branch
proof.

## Verdict boundary

Unit A and the exact steady reduction are established.  Route B1 is blocked on
the exact complement theorem HSE and one named analytic determinant.  The
conditional IFT implication is exact,
but no charged branch, persistence, P2, P4, P5, electron, or neutrino result is
claimed.  No production numerics were run or licensed.
