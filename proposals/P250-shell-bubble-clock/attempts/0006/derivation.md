# Attempt 0006 — exact stitching and repaired spherical bag family

This review attempt repairs the scientific joints exposed by independent
review of PR #196. It preserves attempts 0001–0005 as immutable route evidence
and replaces only the promoted interpretation of their affected outputs.

## Exact wall-sector corrections

`exact_stitching.py` derives every statement below from the canonical action.

- The aligned real-`psi` slice is an exact invariant reduction because all
  four transverse first derivatives of the rotating-frame potential vanish.
  This does not prove that every full-theory planar wall lies on the slice.
- The two off-axis shear Hessians differ and both include a
  `-omega_sq` term. The single positive expression archived in attempt 0001
  is the first static shear channel only, so it cannot support a universal
  transverse-positivity claim along the rotating wall.
- A spatial orbit angle has exact gradient excess
  `iota(q')^2/2`, with `iota=f^2+4b^2`. Thus the phase-slip estimate must bound
  the full clock-generator amplitude, not the scalar amplitude alone.
- Along a stationary bulk branch,
  `p=-V_omega(B(omega))` and `dp/d(omega_sq)=iota_B/2`. Therefore pressure is
  the integral of the changing inertia from the crossing; the endpoint
  product `delta*iota_B(delta)/2` is not exact.
- If `F_omega=E0-omega^2 I/2` and `Q=omega I`, then the physical energy is
  `E=F_omega+omega Q`, whereas `F_omega+omega^2 I/2` is only `E0`. This repairs
  the factor-of-two energy identification in the original C-M5W-005 wording.
- The radial identity is
  `d[r^3(T-V)]/dr=-r^2(T+3V)`, so the Derrick integral is
  `int r^2(T+3V)dr=0` rather than the archived `int r(T+V)dr` diagnostic.

The exact script also reproduces the canonical Maxwell equations with their
normalizations. Its first execution is retained in
`exact_stitching_failure.txt`: a visually identical but assumption-distinct
SymPy `f` symbol failed the Maxwell comparison. Reusing the symbols returned
by the canonical API repairs the implementation, and all checks then pass.

## Numerical representation repair

The old fixed-window Dirichlet BVP left the wall translation mode nearly
singular, hit the mesh ceiling on six of seven rungs, and introduced a clamp
derivative discontinuity. `bag_bvp_repair.py` instead solves in the moving
coordinate `z=r-R`, treats `R` as a BVP parameter, imposes an auxiliary phase
condition, and uses the positive matrix square roots of the exact endpoint
linearizations for radial Robin data.

All seven production rungs converge with solver status 0 and maximum
collocation residual at most `1.00e-8`. The largest relative radial Derrick
residual is `3.25e-10`. At the middle rung, changing the collocation tolerance
from `1e-6` through `1e-10` moves `R_m` by `1.43e-7` relative; changing the
matching exponent from 12 through 16 moves it by `7.31e-10` relative.

Using exact canonical pressure, the radius falls from `1218.761785` to
`174.138125`, and `chi=R_m p/(2 sigma_0)` runs from `1.000328` to `1.002286`.
The endpoint log-log slope of `|chi-1|` is `0.99768`, the expected leading
curvature order. Physical charge `Q=omega I` falls strictly from
`2.343004246e10` to `6.862806275e7`. The five centered finite-difference
checks of `dE/dQ=omega`, using `E=F_omega+omega Q`, have maximum relative
error `1.847e-4`. The first-rung rotating functional agrees with
`16 pi sigma_0^3/(3p^2)` to `0.15494%`.

## Route verdicts

- `WALL_EXACT_SLICE`: established as an invariant reduction with the exact
  first integral; the original universal full-wall quantifier is withdrawn.
- `MAXWELL_LOCAL_CROSSING`: established exactly and locally by the attempt
  0001 interval certificate. Exact global equality of the branch crossing
  with the global ratio infimum is outside this campaign boundary by explicit
  owner direction and remains a separate successor campaign.
- `LOCUS_PHASE_SLIP`: established as exact H1 energy bookkeeping when the
  interface trace lies on the orbit-fixed locus; no global weak-solution or
  stationary two-frequency-wall claim is made.
- `SPHERICAL_BAG_WINDOW`: established numerically on the diagonal slice for
  the seven displayed rungs by the repaired phase-fixed Robin BVP.
- `RADIAL_STABILITY_SPLIT`: the fixed-frequency thin-wall radius is an exact
  radial saddle, and the repaired family has six negative charge secants.
  Neither fact is promoted to a fixed-charge constrained-minimum theorem.
