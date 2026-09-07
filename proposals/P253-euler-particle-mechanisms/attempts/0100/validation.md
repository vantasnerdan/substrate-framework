# Validation ledger

The exact verifier checks eleven derived predicates:

1. the arbitrary-matrix Euler stretching and transported-gradient terms in
   `D_t(omega dot grad chi)` cancel;
2. the product rule gives `div(chi omega)=q_E+chi div omega`, and zero boundary
   flux gives zero total charge in the regular divergence-free class;
3. `D_t q_E=(curl f) dot grad chi` is a divergence and fixes the conserved
   flux sign as `q_E u-f cross grad chi`;
4. a true scalar tag gives pseudoscalar density/axial spatial current, while a
   pseudoscalar tag reverses that parity;
5. the preserved constitutive lock has residual `-chi D_t lambda`;
6. the forced lock has the distinct source condition
   `chi D_t lambda=(curl f) dot grad chi`;
7. the closed-vorticity-line multiplier is the exponential of the integrated
   lock coefficient; and
8. `omega=r zeta e_theta=zeta partial_theta` annihilates every axisymmetric
   tag;
9. the charged Cao correction has the exact sign and curl-superpotential form
   and therefore zero divergence; and
10. the circle-valued azimuthal phase on the punctured domain gives
   `dtheta(omega)=zeta` and equilibrium integral `2*pi*kappa`; and
11. the fixed geometric azimuth has `D_t theta=u_theta/r`, so it is a
    material label only on the no-swirl subspace.

The focused API tests repeat these identities at the import boundary and check
that an untyped parity is rejected.  They are exact symbolic/unit checks, not
evidence for a persistent nonaxisymmetric carrier, a dynamically advected
moving defect, a Lorentz
chiral representation, or P6.

No production numerics are used.  Repository schema validation checks artifact
agreement and workflow structure; it does not prove the scientific identities.
