# Validation and claim boundary

The analytic derivation is load-bearing for the route verdicts.  The focused
API tests and `verify_action_selection.py` expose the finite identities that
would otherwise permit false charge/action selection claims.

The verifier independently checks:

1. the physical dimensions of `Q^2/(epsilon_EM*c_EM)`;
2. invariance of `c_EM^2` and `g_tag^2/epsilon_EM` under gauge-field rescaling;
3. the continuous transported-tag charge family;
4. exact compact-phase gauge cancellation and signed integer-character row in
   the typed convention `Q=g_0*m*N_phase=kappa_m*J_phase`,
   `J_phase=S_0*N_phase`;
5. `H^2(S^3)=0`, `H^2(S^2)=Z`, and the smooth whole-space flux row;
6. the fixed-charge independent action column of the Jacobian;
7. the same nonzero imposed action root for positive and negative
   charge-per-action coefficient `kappa_m`;
8. failure of ordinary Euler vorticity to supply compact BF data.

The public tests additionally reject a zero charge-per-action coefficient in
the imposed-root helper while accepting both signs, and state that the field
rescaling input is a positive normalization magnitude rather than a
charge-sign convention.

These checks do not quantize a classical phase momentum, build a compact
two-form gauge field, select a BF level, or derive a nonzero carrier action.
They verify the algebra and the exact counterdirections used by the analytic
argument.  No production numerics or small-ratio claim is present.

The parent campaign stays active.  In particular this attempt neither
supplies nor replaces the 0088 doublet controls, the source-specific outgoing
Feshbach matrix, a nonlinear persistent invariant/center-stable manifold, a
different coercive carrier, Born/reset dynamics, exchange statistics, or the
neutrino sector.
