# 0270 materialization repair

The original implementation is preserved as failed evidence by the prior
`materialization-review.md`: its `j` formula incorrectly used kinetic energy
`E` and its test copied that false result.  The repaired API now returns

`j = 2*rho*mean_tag_fraction*Q/(3*cell_volume)`,
`a = E/(3*cell_volume)`, and `ell_tag_sq = Q/tag_volume`,

with explicit positivity checks for known nonpositive density, volumes, and
fraction.  Unknown symbolic signs remain caller hypotheses.  The docstring
keeps these as literal geometric/energy normalizations, not a constitutive
law or PDE-existence result.

The tests add an independent sensitivity check: changing `E` leaves `j`
unchanged, while changing the supplied pressure in the analytic solid
rotation exposes the radial centrifugal residual.  A nonzero-poloidal test
(`psi=r*z`, `I=r**2`) gives `u_r=-1`, `u_phi=r` and residual
`R_phi=-2`, explicitly detecting the `u_r*u_phi/r` cylindrical metric term.
Known invalid volumes/density are rejected.

Repaired hashes: module
`1a6e4795ae77fd8e946277377902147fc90b5035fe26edf4d126c2aff61e6d73`, tests
`fb39d6244312967606d065ce53a32b7ed5d1086c4ac7a086b6eabe07734bdd7c`.
Focused pytest was run locally because Herdr pane `w3:p2` was unavailable;
stdout hash `39aee92ac12098bf99d6b37b35c290a45bbeb2b70850c6c63fa392a791e7402b`,
exit hash `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

`route_verdict: repaired conditional additive materialization`.
