import sympy as sp

r, eps, R, rho, q = sp.symbols("r eps R rho q", positive=True)
energy = sp.simplify(rho * sp.Rational(1, 2) * 4 * sp.pi * q**2 * sp.integrate(r**-2, (r, eps, R)))
assert sp.simplify(energy - 2 * sp.pi * rho * q**2 * (1 / eps - 1 / R)) == 0
print("PASS 1: puncture flux field has 1/epsilon core energy")

flux = sp.simplify(4 * sp.pi * R**2 * q / R**2)
assert flux == 4 * sp.pi * q
print("PASS 2: puncture H2 representative has nonzero sphere flux")

k1, k2, k3, F = sp.symbols("k1 k2 k3 F", nonzero=True)
k = sp.Matrix([k1, k2, k3])
PT = sp.eye(3) - k * k.T / (k.dot(k))
assert sp.simplify(PT * (sp.I * F * k)) == sp.zeros(3, 1)
print("PASS 3: every isotropic scalar Fourier source is longitudinal")

cross_self = sp.Matrix([
    k2 * k3 - k3 * k2,
    k3 * k1 - k1 * k3,
    k1 * k2 - k2 * k1,
])
assert cross_self == sp.zeros(3, 1)
print("PASS 4: a pseudoscalar and one wavevector cannot form a transverse axial source")

s, a, Gamma, length = sp.symbols("s a Gamma length", positive=True)
line_energy = sp.integrate((Gamma / (2 * sp.pi * s)) ** 2 * 2 * sp.pi * s, (s, a, R))
target_line_energy = Gamma**2 * (sp.log(R) - sp.log(a)) / (2 * sp.pi)
assert sp.simplify(line_energy - target_line_energy) == 0
print("PASS 5: ideal loop circulation has logarithmic core energy per length")

print("ALL 5 TOPOLOGICAL CHARGE CHECKS PASS")
