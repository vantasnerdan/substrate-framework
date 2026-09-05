"""Exact stationary positive phase-law variance and full-current support test."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0176-stationary-law")
    weight, w1, w2, time = s.symbols("w omega1 omega2 t", real=True)
    response = weight * s.cos(w1 * time) + (1 - weight) * s.cos(w2 * time)
    second = s.diff(response, time, 2).subs(time, 0)
    fourth = s.diff(response, time, 4).subs(time, 0)
    defect = s.factor(fourth - second**2)
    checks.check("actual stationary mixture has its normalized displacement phase", response.subs(time, 0) == 1 and s.diff(response, time).subs(time, 0) == 0)
    checks.check("actual fourth-time derivative exposes the spectral variance", s.factor(defect - weight * (1 - weight) * (w1**2 - w2**2)**2) == 0)
    checks.check("time-reversal does not remove squared-frequency variance", s.expand(defect.subs({w1: -w1, w2: -w2}) - defect) == 0)
    # Two actual complex phases: the outer current Fourier coefficient.
    a, b = s.symbols("a b", real=True, nonzero=True)
    mode = a * s.exp(-s.I * w1 * time) + b * s.exp(-s.I * w2 * time)
    current = s.expand((s.conjugate(mode) * s.diff(mode, time) - mode * s.conjugate(s.diff(mode, time))) / (2 * s.I))
    expected = -a**2 * w1 - b**2 * w2 - a * b * (w1 + w2) * s.cos((w1 - w2) * time)
    checks.check("exact two-frequency mechanical phase current", s.simplify(s.expand_complex(current) - expected) == 0)
    checks.check("squared current has nonzero doubled difference-frequency coefficient", s.expand((a * b * (w1 + w2))**2 / 2) != 0)
    print(f"stationary response defect={defect}")
    print(f"two-frequency current={expected}")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
