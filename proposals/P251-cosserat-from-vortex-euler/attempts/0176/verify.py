"""Exact positive-marker band edge and exposed all-time phase coefficients."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0176")
    p, z, h, r, c, om, energy = s.symbols("P z h r c Omega E", positive=True)
    time = s.symbols("t", real=True)
    scale = p**s.Rational(3, 2)
    shift = z * (scale - 1)
    weight_plus, weight_minus = (r + 1) / 2, (r - 1) / 2
    denominator = weight_plus * (shift + h * scale) + weight_minus * (shift - h * scale)
    numerator = weight_plus * (z + h) * (shift + h * scale) + weight_minus * (z - h) * (shift - h * scale)
    mean = s.cancel(numerator / denominator)
    ratio = (shift + r * h * scale) / (r * shift + h * scale)
    checks.check("actual nodal quotient from two positive weights", s.simplify(mean - z - h * ratio) == 0)
    checks.check("finite physical mean at reference carrier", s.simplify(mean.subs(p, 1) - z - r * h) == 0)
    mean_prime = s.simplify(s.diff(mean, p).subs(p, 1))
    checks.check("first fixed-radius carrier derivative", mean_prime == 3 * z * (1 - r**2) / 2)
    checks.check("actual singular second susceptibility", s.simplify(s.limit(h * s.diff(mean, p, 2).subs(p, 1), h, 0) - s.Rational(9, 2) * r * (r**2 - 1) * z**2) == 0)
    sigma = -2 * om + energy * c / s.sqrt(p)
    gamma_initial = sigma - c * p * mean / 2
    band_polynomial = 3 * z * r**2 / 2 - h * r - energy - 5 * z / 2
    gamma_prime = s.simplify(s.diff(gamma_initial, p).subs(p, 1))
    checks.check("finite-h band equation derived from actual phase", s.simplify(gamma_prime - c * band_polynomial / 2) == 0)
    root = (h + s.sqrt(h**2 + 6 * z * (energy + 5 * z / 2))) / (3 * z)
    checks.check("positive branch solves finite-h band equation", s.simplify(band_polynomial.subs(r, root)) == 0)
    checks.check("band equation has simple r derivative", s.diff(band_polynomial, r) == 3 * z * r - h)
    checks.check("seed r above one follows from equation at one", s.expand(band_polynomial.subs(r, 1)) == -energy - h - z)
    leading_curvature = s.diff(s.diff(gamma_initial**2, p, 2), c).subs({p: 1, c: 0})
    checks.check("band curvature has positive leading coefficient", s.simplify(s.limit(h * leading_curvature, h, 0) - 9 * om * r * (r**2 - 1) * z**2) == 0)

    angle = c * p * h * time / 2
    aa = r * shift + h * scale
    bb = shift + r * h * scale
    real_part = aa * s.cos(angle)
    imag_part = -bb * s.sin(angle)
    frequency = s.simplify((real_part * s.diff(imag_part, time) - imag_part * s.diff(real_part, time)) / (real_part**2 + imag_part**2))
    dnorm = s.cos(angle)**2 + ratio**2 * s.sin(angle)**2
    expected_frequency = -c * p * h * ratio / (2 * dnorm)
    checks.check("full interference phase is obtained from actual complex numerator", s.simplify(frequency - expected_frequency) == 0)
    amplitude_rate = s.simplify((real_part * s.diff(real_part, time) + imag_part * s.diff(imag_part, time)) / (real_part**2 + imag_part**2))
    expected_rate = c * p * h * (ratio**2 - 1) * s.sin(angle) * s.cos(angle) / (2 * dnorm)
    checks.check("full amplitude connection derived independently", s.simplify(amplitude_rate - expected_rate) == 0)
    t2 = c**3 * p**3 * h**3 * ratio * (ratio**2 - 1) / 8
    derivative_t2 = s.simplify(s.diff(t2, p).subs(p, 1))
    band_energy = 3 * z * r**2 / 2 - h * r - 5 * z / 2
    expected_t2 = -3 * c**3 * h**2 * (r**2 - 1) * (energy + 2 * z) / 8
    checks.check("exact initial band edge has nonzero t-squared carrier drift", s.simplify(derivative_t2 - expected_t2.subs(energy, band_energy)) == 0)
    checks.check("initial band edge does not remove amplitude-rate derivative", s.simplify(s.diff(amplitude_rate, time).subs({p: 1, time: 0}) - c**2 * h**2 * (r**2 - 1) / 4) == 0)

    eps, tail, negative_weight = s.symbols("e tau A", positive=True)
    locations = (-s.Integer(1), eps, s.Integer(1))
    weights = (negative_weight, s.Integer(1), tail)
    moments = [sum(w * loc**order for w, loc in zip(weights, locations, strict=True)) for order in range(5)]
    variance_numerator = s.expand(moments[1] * moments[3] - moments[2]**2)
    repaired_weight = tail * eps * (1 - eps)**2 / (eps * (1 + eps)**2 + 4 * tail)
    checks.check("third positive annulus solves actual signed response variance", s.factor(variance_numerator.subs(negative_weight, repaired_weight)) == 0)
    carrier_strength = s.factor(moments[0] * moments[2] / moments[1]**2)
    tscale = s.symbols("b", positive=True)
    repaired_strength = carrier_strength.subs(negative_weight, repaired_weight).subs(tail, tscale * eps**2)
    checks.check("variance-free positive marker retains a tunable initial slope", s.limit(repaired_strength, eps, 0) == 1 + 2 * tscale)
    response_mean = moments[2] / moments[1]
    third_central = s.factor(moments[4] / moments[1] - 3 * response_mean * moments[3] / moments[1] + 2 * response_mean**3)
    third_repaired = s.factor(third_central.subs(negative_weight, repaired_weight))
    checks.check("zero variance leaves a nonzero third phase cumulant", third_repaired != 0)
    # This identity uses a POSITIVE material measure, not a signed probability.
    mu = s.symbols("mu", real=True)
    central2 = sum(w * loc * (loc - mu)**2 for w, loc in zip(weights, locations, strict=True))
    central3 = sum(w * loc * (loc - mu)**3 for w, loc in zip(weights, locations, strict=True))
    positive_form = sum(w * loc**2 * (loc - mu)**2 for w, loc in zip(weights, locations, strict=True))
    checks.check("variance and third cumulant imply an exact positive rigidity form", s.expand(central3 + mu * central2 - positive_form) == 0)

    print(f"finite band equation: {band_polynomial}")
    print(f"initial-band time drift gamma_P,t2={expected_t2}")
    print(f"three-annulus positive weight={repaired_weight}")
    print(f"remaining third signed cumulant={third_repaired}")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
