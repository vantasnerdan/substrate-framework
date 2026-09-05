# Exposing-tangent repair after the finite-space extension

The first thirteen action/current predicates passed in first.stdout.
The finite-space extension then passed fifteen predicates and failed the
conjunction claiming its chosen preparation slip was nonzero AND had zero
observed stress. For kappa=(3/5,0,4/5), D=(4/5,0,-3/5), the actual
P_n[aD-kappa(u.D)] is identically zero: D and kappa have parallel
transverse projections. The stress cancellation was not false; the probe
did not expose its nontrivial case. Use D=(0,1,0) for that check, giving
the actual nonzero slip (-3 sin z/5,3 cos z/5,0). Its stress row still
vanishes exactly. The nilpotent-subspace checks retain their original data.

finite-space.stdout preserves the fifteen passes. The next energy-extension
run was launched before the earlier process's failure was polled and hit
the same unchanged assertion; isotropic-energy.stdout preserves it too.
Both exited1 with CheckFailure on the nonzero-slip conjunction. The repair
changes only this exposing test vector, not an equation or claim scope.
