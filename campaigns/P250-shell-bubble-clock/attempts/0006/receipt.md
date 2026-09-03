# Attempt 0006 — review repair receipt for the spherical bag family

Parent: PR #196 review of C-M5W-007 and its C-M5W-008 consumer.

## Positive object

A seven-rung stationary spherical bag family on the C-M5W-001 diagonal slice,
at `omega^2 = omega_*^2 + delta`, `delta = 0.001,...,0.007`, with every rung
anchored to `solve_bvp` status 0 and pressure evaluated directly from the
canonical rotating-frame potential.

## Failure mechanisms repaired

1. Attempt 0003 used `delta*iota_B(omega)/2` as an exact pressure.  Along a
   stationary branch the envelope theorem instead gives
   `dp/d(omega^2)=iota_B/2`; the pressure is the integral of this derivative,
   or equivalently `-V_omega(B(omega))`.  The endpoint product is only a
   first-order approximation.
2. A fixed, Dirichlet-clamped wall window left the translational wall mode
   nearly singular.  Mesh growth reached the node ceiling on six rungs and
   the artificial clamps created a derivative jump against the analytic flat
   interior.
3. The archived producer formed `Q_trap = omega^2 I` while the accepted charge
   is `Q = omega I`.
4. The archived virial diagnostic integrated `r(T+V)`.  The radial Derrick
   identity for `E_omega = 4 pi int r^2(T+V)dr` is
   `int r^2(T+3V)dr = 0`, equivalently
   `d[r^3(T-V)]/dr = -r^2(T+3V)` on solutions.

## Frozen numerical representation

Use the moving wall coordinate `z=r-R` and solve for the translation parameter
`R` together with the eight first-order wall variables.  An auxiliary phase
condition makes the translation mode nonsingular.  At both window ends use
the exact linear radial asymptotic logarithmic derivatives: the regular
interior mode has `u'=(sqrt(D_B)-1/r)(u-B)` up to exponentially small
origin reflection, while the exterior mode has
`u'=-(sqrt(D_A)+1/r)(u-A)`, with `D=(2K)^(-1) Hess(V_omega)`.

The production matching amplitude is `exp(-14)`.  Boundary sensitivity is
measured at matching exponents 12, 14, and 16; collocation sensitivity is
measured at tolerances `1e-6`, `1e-8`, and `1e-10`.  Every quoted production
rung must have solver status 0, finite output, maximum collocation residual
below `1.1e-8`, boundary residual below `1e-10`, and consistent branch identity.

Pressure is `p=-V_omega(B(omega))` from the canonical API.  The fixed-frequency
functional, inertia, physical charge, and physical energy are

`E_omega=4 pi int r^2(T+V_omega)dr`,
`I=4 pi int r^2 iota dr`,
`Q=omega I`, and `E=E_omega+omega Q`.

Flat-bulk contributions from `0` to the inner matching radius are integrated
analytically.  The omitted nonlinear tails are bounded by the measured
matching-exponent sensitivity; no large self-energy subtraction is used.

## Maximum verdict

`BAG_EXISTS_NUMERICALLY` on the spherical diagonal slice and the displayed
seven-rung window.  This does not license a full-space minimum, a constrained
stability theorem, non-spherical modes, dynamics, or the thick-wall regime.
