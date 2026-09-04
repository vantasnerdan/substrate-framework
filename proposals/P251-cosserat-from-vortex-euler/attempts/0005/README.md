# Attempt 0005 — N2 block B2 (twist stiffness), gauge-identity result and route pin

## Route probed

Bundle line-stretch heuristic: each vortex line at radius r winds at rate
chi' = k; its arc length per z grows by (1 + chi'^2 r^2 / 2); multiplying the
local core energy by that stretch suggested Delta E_local/L = rho G^2 chi'^2
/ (8 pi). DIMENSIONAL CHECK refutes this route: rho G^2 chi'^2 has units
force/length^2, not force (energy/length) -- the stretch heuristic
double-counts the shared-core local field, which to O(k^2 a^2) is the field
of a straight line (long-wave limit).

## Gauge-identity mechanism (the load-bearing finding)

For an AXISYMMETRIC continuum bundle, the rigid twist at rate chi' maps the
configuration at height z to the initial configuration rotated by chi' z.
The axisymmetric bundle is rotation-invariant, so the twisted state is the
SAME state: the twist energy of the axisymmetric bundle is identically
zero. Combined with attempt 0004 (axisymmetric redistribution gives
negative Delta E -- a flux-redistribution artifact), this establishes:

  The twist degree of freedom is NOT representable in any axisymmetric
  sector. It is carried by the core POLARIZATION -- the m=2 (elliptical)
  surface mode of the tube -- whose rotation rate is the microrotation
  variable Phi. This is consistent with the Comparsi intake decomposition
  v = V_macro + Phi x r + v' (spin of a LOCALIZED, polarized core).

## Consequence for the preregistered candidates

Both candidates converge on the polarization carrier and differ only in
what carries it: candidate A = triad/frame polarization of the resolved
line structure; candidate B = m=2 core-surface polarization (Rankine
vortex). The frozen selection criteria and N3/N4 identification tests are
unchanged. This eliminates the axisymmetric sub-family from the candidate
universe (append-only expansion: "axisymmetric twist sub-family, refuted
by gauge identity, 0005").

## Next route (0006, well-posed closed form)

Exact m=2 surface mode of a Rankine vortex tube (vorticity 2Omega = G/(pi
a^2) inside, potential flow outside, vortex-sheet boundary at r = a):
linearized Euler inside, Laplace outside, kinematic + dynamic boundary
conditions; the dispersion relation is the classic Kelvin 1880 problem and
is ELEMENTARY (algebraic in omega for each m, k). The twist modulus:
second variation of the mode energy with respect to the polarization
rotation rate; the microrotation inertia J_i is the matching moment of the
same core. SymPy-derivable end to end; two-route check via the energy of a
statically elliptical core (quadrupole line tension).

## Status

Attempt 0005 closes with a structural result (gauge identity + route pin),
not a modulus. N2 remains active.
