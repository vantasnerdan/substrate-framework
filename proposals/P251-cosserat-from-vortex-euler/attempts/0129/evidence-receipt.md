# Bounded result and exact receipts

`inertial-centroid.md` SHA256:
`b689cab53d0f95fef06679d4d7753df51e0e01ac3265a288ed4ee88a6f65c1cd`.

`periodic-cell-continuation.md` SHA256:
`bb7dc305660899bc5318fa3fe7ce3b34d81a15df7b7707743da9695ae000a80f`.

`verify.py` SHA256:
`fcb21f60bb7fe84d5e25ea8a282e661ec28f8c81e5b2e1c2f651f4c7078f2843`.

The initial Bessel/Euler/centroid/action oracle passed 25 checks. After
the failure-derived cell/Rayleigh construction, the expanded oracle
passed all 39 checks on its first execution. Both runs are preserved.
Canonical `CheckLedger` supplies checks and status-zero completion;
`moving_phase_pullback` supplies the exact physical frame/action map.
Ruff and all new-file no-index whitespace checks pass (exit 1 for a
new-file diff, with no whitespace diagnostics, is not a failed check).

The analytic oracle establishes actual equations, pressure and material
maps before dispersion is interpreted. The periodic ODE inverse and
analytic implicit-function proof supply existence of the helical
growing/decaying Euler branch; the symbolic tally is not used as a
substitute for that argument. No numerical roots, energy-floor signs,
or stability-window searches were performed.

Strongest results:

- Exact low-frequency m=1 Euler/Lin Bessel family, a nonzero observed
  subparcel centroid, and its correct slip-boundary sideband.
- Exact mechanism preventing that fixed family from being a positive
  finite-mass acoustic continuum: nonzero KKS at vanishing physical
  momentum, laboratory Hamiltonian sign change and complete-carrier
  mean cancellation.
- Actual periodic translation Jordan chain and the first Euler cell
  equation with its physical mean-momentum stress coefficient.
- An explicit single-helical-field cell solution and exact Euler/Kelvin
  continuation with omega=+/-i V K/sqrt(2)+O(K^3), so this DIFFERENT
  physical acoustic candidate has negative rather than positive stiffness.

These are two route-specific scientific mechanisms and reusable exact
calculations, not an exhaustion result for stationary Euler ensembles.
The positive acoustic response of other declared cellular, multiwave or
vortex-array backgrounds remains an active parent construction.
