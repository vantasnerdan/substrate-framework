"""Derive correlated observed rows, positive phase and actual energy equality."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0187-correlated-family-action")
    nu, t, q, k, channel = s.symbols("nu t q k c", positive=True)
    amp, slope, curvature, speed, bend = s.symbols("a d e v b", real=True)
    profile = amp+slope*q+curvature*q*q/2
    frequency = nu+speed*q+bend*q*q/2
    for name, row in (("angle", s.cos(frequency*t)),
                      ("rate", s.sin(frequency*t)/frequency)):
        derived = s.diff(profile*row, q, 2).subs(q, 0)
        base = row.subs(q, 0)
        target = curvature*base+(amp*bend+2*slope*speed)*s.diff(base, nu)
        target += amp*speed**2*s.diff(base, nu, 2)
        checks.check(f"actual carrier derivative of {name} observation",
                     s.simplify(derived-target) == 0)
    b_eff = s.symbols("B", real=True)
    f = s.cos(nu*t)+channel*k*k*b_eff*s.diff(s.cos(nu*t), nu)/2
    g = s.sin(nu*t)/nu+channel*k*k*b_eff*s.diff(s.sin(nu*t)/nu, nu)/2
    wronskian = s.expand(f*s.diff(g, t)-g*s.diff(f, t))
    checks.check("signed variance cancellation gives a stationary full-time Wronskian",
                 s.simplify(wronskian.coeff(k, 0)-1) == 0
                 and s.simplify(wronskian.coeff(k, 2)) == 0)
    effective_gap = nu*nu+channel*k*k*nu*b_eff
    for name, row in (("angle", f), ("rate", g)):
        residual = s.expand(s.diff(row, t, 2)+effective_gap*row)
        checks.check(f"observed {name} obeys the stationary equation through spatial two",
                     s.simplify(residual.coeff(k, 0)) == 0
                     and s.simplify(residual.coeff(k, 2)) == 0)
    j, jp, jpp = s.symbols("j jprime jsecond", real=True)
    mass = (j+jp*q+jpp*q*q/2)*profile**2
    energy_difference = s.diff(mass*frequency**2, q, 2).subs(q, 0)
    energy_difference -= nu**2*s.diff(mass, q, 2).subs(q, 0)
    target = 2*j*amp**2*(speed**2+nu*bend)
    target += 4*nu*(jp*amp**2+2*j*amp*slope)*speed
    checks.check("complete physical energy keeps mass and amplitude derivatives",
                 s.expand(energy_difference-target) == 0)
    w1, w2, a1, a2, j1, j2, v1, v2, total = s.symbols(
        "w1 w2 a1 a2 j1 j2 v1 v2 J", nonzero=True)
    coefficient = s.Matrix([[2*w1*v1, 2*w2*v2],
                           [2*nu*w1*v1*(2*j1*a1-total),
                            2*nu*w2*v2*(2*j2*a2-total)]])
    checks.check("curvature and inherited energy can be matched independently",
                 s.factor(coefficient.det()-8*nu*w1*w2*v1*v2*(j2*a2-j1*a1)) == 0)
    amplitudes = [s.Rational(-2, 3), s.Rational(8, 3)]
    derivatives = [s.Rational(9, 20), s.Rational(1, 10)]
    speeds = [s.Integer(2), s.Integer(1)]
    observed = sum(amplitudes)/2
    variance = sum(ai*vi**2 for ai, vi in zip(amplitudes, speeds, strict=True))/2
    response_b = sum(di*vi for di, vi in zip(derivatives, speeds, strict=True))
    phase = sum(ai**2 for ai in amplitudes)/2
    energy = sum(ai**2*vi**2+4*ai*di*vi for ai, di, vi in
                 zip(amplitudes, derivatives, speeds, strict=True))/2
    checks.check("positive-probability example cancels observed variance",
                 observed == 1 and variance == 0 and phase > 0)
    checks.check("same example retains positive curvature and the actual Euler-energy interface",
                 response_b == 1 and energy == phase*response_b)
    checks.check("freezing amplitude slopes would expose false energy and curvature claims",
                 sum(ai**2*vi**2 for ai, vi in zip(amplitudes, speeds, strict=True))/2 > 0)
    print("Derived example initial phase:", phase)
    print("Derived example gradient phase:", sum(di**2 for di in derivatives))
    print("Actual Euler mode-family jets and same-field geometry remain explicit inputs.")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
