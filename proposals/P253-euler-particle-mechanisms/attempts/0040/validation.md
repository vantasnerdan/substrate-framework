# Validation receipt

## Exact oracle

The frozen receipt execution used the repository interpreter:

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P253-euler-particle-mechanisms/attempts/0040/verify_kelvin_bridge.py

Its command, stdout, empty stderr, and exit `0` are preserved in the
`kelvin-bridge.*` quartet.  The verifier checks twenty-five exact predicates:

1. both source filament coefficients have their stated logarithmic slopes;
2. both fixed-mode large-log limits are correct;
3. their product gives the stated Kelvin frequency;
4. the source canonical Hamiltonian reproduces both oscillator equations;
5. direct angular integration gives the nonzero quadratic `J_z` coefficient;
6. the rotating-frame sign annihilates the mode at `Omega=-nu/l`;
7. the CR parameter derivative is the nonzero multiplier `-i*l`;
8. the physical push-forward gives `X_J=-[R,omega]`, hence
   `dJ_z/da=+l*sigma_l*a`, and the explicit filament KKS coefficient gives
   exactly the same positive `J_z` coefficient as direct angular integration;
9. the vorticity-side moment-map potential curls to physical axial rotation;
10. the Cao Lane--Emden equation implies the profile identities used to derive
   the Richardson function and its dimensionless H2 sign criterion;
11. the radial Lane--Emden coefficient recurrence through the cubic spatial
    jet generates the quoted `y`, `t`, and positive H2-comparison center
    coefficient; any hypothetical first zero then has positive outward
    derivative, and the resulting bound makes the Richardson logarithmic
    derivative strictly negative; and
12. the exterior `m=1` modified-Bessel jet satisfies its differential equation
    to the retained order and its exact Dirichlet-to-Neumann ratio has the
    universal `k^2[log(k a/2)+gamma_E]` threshold coefficient.

The oracle is exact SymPy algebra, not production numerics.  It deliberately
does not claim the missing spectral projection, strict-tail regularization,
curved-core inverse, nonlinear branch, or quantum interpretation.
The author-stage sign repair, the exposing order-truncation failure, and the
new exterior-threshold result are recorded append-only in
`author-sign-and-exterior-receipt.md`.

## Exposing analytic checks

- Setting `l=1` removes the leading logarithm from `B_l`, exposing the
  translation sector; the branch target therefore retains `l>=2`.
- Dropping the minus sign in `Omega_0=-nu/l` leaves the kernel scalar nonzero.
- Omitting the phase-quadrature factor `nu/B_l` makes (20) fail Hamilton's
  equations and changes the angular moment.
- A full-space filament/operator norm comparison is rejected analytically
  because the exact column has an essential interval and fast internal core
  modes absent from the two-coordinate model.
- Deleting the logarithmic term in the exterior Bessel jet makes both its
  differential-equation residual and Dirichlet-to-Neumann limit fail.
- The helical contour construction is rejected as a compact carrier by the
  exact identity `H_(2*pi)(x)=x+2*pi*h*e_z`, not by a numerical observation.

`git diff --check` is the formatting oracle.  No central file, source module,
accepted API, test, or artifact outside attempts 0032 and 0040 was edited.
