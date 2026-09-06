"""Exact checks for the P253/0037 forced threshold normal form."""

import sympy as sp


def main() -> None:
    x, beta, sigma, v = sp.symbols("x beta sigma v", positive=True)
    profile = 3 * sp.sech(x / 2) ** 2 / (2 * beta)

    # The exact profile equation forced by the 0027 branch.
    assert sp.simplify(profile - sp.diff(profile, x, 2) - beta * profile**2) == 0

    # A speed-sigma traveling wave of the derived KdV equation is that profile.
    traveling = -sigma * sp.diff(profile, x) + sigma * sp.diff(
        sp.diff(profile, x, 2) + beta * profile**2, x
    )
    assert sp.simplify(traveling) == 0

    # Exact mass and translation-momentum normalizations.  Use
    # t=tanh(x/2), dx=2 dt/(1-t^2), so no improper hyperbolic integral is
    # delegated to a heuristic integration branch.
    t = sp.symbols("t", real=True)
    profile_t = 3 * (1 - t**2) / (2 * beta)
    dx_dt = 2 / (1 - t**2)
    mass = sp.integrate(profile_t * dx_dt, (t, -1, 1))
    momentum = sp.integrate(profile_t**2 * dx_dt / 2, (t, -1, 1))
    assert sp.simplify(mass - 6 / beta) == 0
    assert sp.simplify(momentum - 3 / beta**2) == 0

    # General-speed invariant slopes are positive.
    mass_v = 6 * sp.sqrt(v / sigma) / beta
    momentum_v = 3 * (v / sigma) ** sp.Rational(3, 2) / beta**2
    assert sp.simplify(sp.diff(mass_v, v)) > 0
    assert sp.simplify(sp.diff(momentum_v, v)) > 0

    # Pöschl--Teller rescaling: 4 L=-d_yy+4-12 sech^2(y).
    # Its n=3 discrete levels are 4-(3-j)^2, j=0,1,2,3.
    levels = [sp.Rational(4 - (3 - j) ** 2, 4) for j in range(4)]
    assert levels == [sp.Rational(-5, 4), 0, sp.Rational(3, 4), 1]

    print("6 exact threshold-reduction checks passed")


if __name__ == "__main__":
    main()
