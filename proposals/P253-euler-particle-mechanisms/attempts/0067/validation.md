# Validation and scope

The final repository-interpreter oracle checks the exact longitudinal/transverse
Oseen contractions, the isotropic trace produced by a declared orthonormal
three-field multiplet, sign reversal for opposite vector sources, finite
smooth-form-factor self energy, and cancellation of carrier/field translation
momentum exchange. These five checks pass with exit `0`. The first run is
preserved with exit `1`: its sign-reversal assertion omitted the common leading
minus sign in the on-shell energy and was repaired by comparing the same-source
and opposite-source energies explicitly. These algebraic checks do not construct the carrier
action, source constitutive law, field triplet, or frame lock.

The fifth check is only the local algebra `exchange+(-exchange)=0`; it does not
prove the field momentum is finite or the translation Noether theorem. The
action variation and Noether identity are derived analytically under
translation invariance, sufficient decay/finite field momentum, and the full
Euler--Lagrange equations. No
production numerics or empirical values are used. Independent review is
required before consuming the route-scoped conditional recoil theorem.
