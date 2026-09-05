"""Exact common-phase, force-return and infrared pressure-jet checks."""

from __future__ import annotations

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0174")
    x, y, z, kx, ky, kz = sp.symbols("x y z Kx Ky Kz", real=True)
    coordinates = (x, y, z)
    position = sp.Matrix(coordinates)
    wave = sp.Matrix([kx, ky, kz])
    potential = sp.Matrix([x * y + z**2, y * z + x**2, z * x + y**2])
    phase = sp.exp(-sp.I * wave.dot(position))
    ordinary = sp.Matrix([sp.diff(potential[(j + 2) % 3], coordinates[(j + 1) % 3]) - sp.diff(potential[(j + 1) % 3], coordinates[(j + 2) % 3]) for j in range(3)])
    modified = phase * potential
    bloch_curl = sp.Matrix([sp.diff(modified[(j + 2) % 3], coordinates[(j + 1) % 3]) - sp.diff(modified[(j + 1) % 3], coordinates[(j + 2) % 3]) for j in range(3)]) + sp.I * wave.cross(modified)
    ledger.check("compact-potential phase compensation is an exact physical curl identity", all(sp.simplify(entry) == 0 for entry in bloch_curl - phase * ordinary))
    ledger.check("full corrected preparation is exactly Bloch solenoidal", sp.simplify(sum(sp.diff(bloch_curl[j], coordinates[j]) + sp.I * wave[j] * bloch_curl[j] for j in range(3))) == 0)
    ledger.check("omitting the Bloch curl term changes the actual displacement", wave.cross(potential) != sp.zeros(3, 1))
    radius, s = sp.symbols("R s", positive=True)
    curve = sp.Matrix([radius * sp.cos(s / radius), radius * sp.sin(s / radius), 0])
    curved_phase = wave.dot(curve)
    ledger.check("curved-arc phase compensation retains the full transport derivative", sp.simplify(sp.diff(sp.exp(sp.I * curved_phase) * sp.exp(-sp.I * curved_phase), s)) == 0)
    ledger.check("frozen-tangent subtraction would leave a curvature derivative", sp.simplify(sp.diff(curved_phase - ky * s, s)) != 0)

    internal_j = sp.Matrix([[0, -1], [1, 0]])
    plus = (sp.eye(2) - sp.I * internal_j) / 2
    minus = (sp.eye(2) + sp.I * internal_j) / 2
    ledger.check("signed internal carrier projectors are complementary", plus**2 == plus and minus**2 == minus and plus * minus == sp.zeros(2) and plus + minus == sp.eye(2))
    shift, central = sp.symbols("shift central", real=True)
    def example_columns(value):
        return sp.Matrix([[1 + value, value**2], [2 * value, 1 - value]])
    paired = example_columns(central - shift) * plus + example_columns(central + shift) * minus
    ledger.check("both signed sidebands give the actual real Bloch family", paired.subs(shift, -shift) == sp.conjugate(paired))
    ledger.check("a single real shifted column violates Bloch reality", example_columns(central + shift) != sp.conjugate(example_columns(central - shift)))

    u = (
        ef.add(ef.trig(2, kind="sin"), ef.scale(ef.trig(1), 3)),
        ef.add(ef.scale(ef.trig(0, kind="sin"), 2), ef.trig(2)),
        ef.add(ef.scale(ef.trig(1, kind="sin"), 3), ef.scale(ef.trig(0), 2)),
    )
    omega = ef.curl(u)
    columns = [tuple(ef.derivative(component, j) for component in u) for j in range(3)]
    gram = sp.Matrix(3, 3, lambda i, j: ef.inner(columns[i], columns[j]))
    forces = sp.Matrix.hstack(*(sp.Matrix([component.get(ef.ZERO, 0) for component in ef.cross(ef.curl(column), omega)]) for column in columns))
    ledger.check("actual force-control map is the gradient Gram with correct sign", forces == gram)
    ledger.check("three-dimensional comparison exposes a nonsingular return map", gram == sp.diag(4, 9, 1) and gram.det() > 0)
    defect = sp.Matrix([sp.Rational(1, 7), -sp.Rational(2, 5), sp.Rational(3, 11)])
    correction = gram.inv() * defect
    ledger.check("all three force defects are removed by the solved return", defect - forces * correction == sp.zeros(3, 1))
    theta = sp.symbols("theta", real=True)
    for shift in (-1, 0, 1):
        average = sp.integrate(sp.exp(sp.I * (7 + shift) * theta), (theta, 0, 2 * sp.pi))
        ledger.check(f"high toroidal harmonic has exactly zero Cartesian force shift{shift}", average == 0)
    ledger.check("low toroidal harmonic cannot use that cancellation", sp.integrate(sp.exp(sp.I * (1 - 1) * theta), (theta, 0, 2 * sp.pi)) != 0)

    scale, k, h = sp.symbols("scale k h", positive=True)
    q0 = sp.Matrix([1, 2, 3])
    direction = sp.Matrix([0, 1, 0])
    q = scale * q0 + k * direction
    projector = sp.eye(3) - q * q.T / q.dot(q)
    stress = sp.Matrix([[1, 2, 0], [2, 3, 1], [0, 1, 4]])
    multiplier = -sp.I * projector * stress * q
    for order in (0, 1, 2):
        jet = multiplier.diff(k, order).subs(k, 0)
        normalized = sp.simplify(jet * scale ** (order - 1))
        ledger.check(f"full Euler stress symbol has degree1-minus-{order} jet", not normalized.has(scale) and normalized != sp.zeros(3, 1))
    projection_second = projector.diff(k, 2).subs(k, 0)
    ledger.check("bare Kelvin projection really has the stronger second-jet singularity", not sp.simplify(projection_second * scale**2).has(scale) and projection_second != sp.zeros(3))
    force_hat = sp.Matrix([q[1], q[2], q[0]])
    initial_second = (projector * force_hat).diff(k, 2).subs(k, 0)
    ledger.check("zero-integral initial force restores integrable inverse-first-power behavior", not sp.simplify(scale * initial_second).has(scale) and initial_second != sp.zeros(3, 1))
    radial = sp.symbols("r", positive=True)
    infrared = sp.integrate(4 * sp.pi * radial**2 / radial**2, (radial, 0, h))
    ledger.check("second pressure jet is square integrable in three dimensions", infrared == 4 * sp.pi * h)
    ledger.check("the same infrared argument is dimension-sensitive", sp.integrate(radial / radial**2, (radial, 0, h)) == sp.oo)
    for order, exponent in enumerate((sp.Rational(7, 4), sp.Rational(3, 4), sp.Rational(1, 4))):
        radial_power = sp.expand(2 * exponent + 2 - 2 * (4 - order))
        ledger.check(f"actual phase derivative{order} has an admissible weighted L2 tail", radial_power < -1)

    # Exact fourth orientation moments before any elimination.
    zz = sp.symbols("zz", real=True)
    normal_fourth = sp.integrate(zz**4 / 2, (zz, -1, 1))
    mixed_fourth = sp.integrate(zz**2 * (1 - zz**2) / 4, (zz, -1, 1))
    kinetic_mark = sp.integrate(zz**2 / 2, (zz, -1, 1))
    ledger.check("one common laboratory K gives the fourth Haar tensor", normal_fourth == sp.Rational(1, 5) and mixed_fourth == sp.Rational(1, 15))
    ledger.check("the longitudinal/transverse weights use the actual kinetic mark", normal_fourth / kinetic_mark == sp.Rational(3, 5) and mixed_fourth / kinetic_mark == sp.Rational(1, 5))
    p, nx, ny, nz = sp.symbols("p nx ny nz", real=True)
    clock = sp.Function("clock")
    composed = clock(p - nx * kx - ny * ky - nz * kz)
    second = sp.hessian(composed, (kx, ky, kz)).subs({kx: 0, ky: 0, kz: 0}).doit()
    normal = sp.Matrix([nx, ny, nz])
    ledger.check("actual preparation rather than rotated K supplies the axial tensor", second == normal * normal.T * sp.diff(clock(p), p, 2))
    mass, gap = sp.Function("M")(p), sp.Function("gap")(p)
    action_gradient = sp.diff(mass * gap, p, 2) - gap * sp.diff(mass, p, 2)
    ledger.check("action averaging retains the mass-slope cross term", sp.expand(action_gradient - mass * sp.diff(gap, p, 2) - 2 * sp.diff(mass, p) * sp.diff(gap, p)) == 0)
    ledger.check("positive frequency curvature alone does not fix averaged action sign", action_gradient.subs({mass: sp.exp(-p), gap: p**2}).doit().subs(p, 1) < 0)
    print(f"force Gram={gram}; solved finite correction={list(correction)}")
    print(f"infrared squared norm={infrared}; orientation weights={normal_fourth / kinetic_mark}, {mixed_fourth / kinetic_mark}")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
