# Validation and claim boundary

The exact repository-interpreter oracle derives the Newton Hessian, contracts
the antisymmetric vorticity moments, and checks the signs and powers of the
impulse velocity, two-carrier cross energy/force, physical-pressure
quadrupole, radial flux, the translation Ward identity for a compact internal stress source, and annihilation of the general local isotropic scalar source by the transverse projector. All nine assertions pass with exit `0`.

The importable API has six focused tests covering axial and transverse dipole
signs, cubic separation scaling, physical density and isotropic pressure
cancellation, monopole flux, and invalid domains. The focused run passes six
tests with exit `0`.

These algebraic checks expose a false Coulomb power, a sign-reversed dipole,
omitted physical density, a hidden flux source, or an isotropic pressure tail.
The analytic derivation supplies the integration-by-parts hypotheses and the
collective-source distinction. Neither the oracle nor the API constructs an
autonomous defect, a massless scalar, charge, magnetic coupling, or a particle.
Independent scientific review is required before consuming the established
route-scoped multipole conclusions.
