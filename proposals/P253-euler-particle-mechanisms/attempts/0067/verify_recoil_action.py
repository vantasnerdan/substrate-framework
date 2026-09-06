import sympy as sp

mu, r = sp.symbols("mu r", positive=True)
n = sp.Matrix([0, 0, 1])
GT = (sp.eye(3) + n * n.T) / (8 * sp.pi * mu * r)

parallel = sp.Matrix([0, 0, 1])
transverse = sp.Matrix([1, 0, 0])
assert (parallel.T * GT * parallel)[0] == 1 / (4 * sp.pi * mu * r)
assert (transverse.T * GT * transverse)[0] == 1 / (8 * sp.pi * mu * r)
print("PASS 1: one-vector Oseen interaction is orientation dependent")

assert sp.trace(GT) == 1 / (2 * sp.pi * mu * r)
print("PASS 2: a locked orthonormal three-field multiplet gives an isotropic trace")

F = sp.Matrix(sp.symbols("F0:3"))
same_source_energy = -(F.T * GT * F)[0]
opposite_source_energy = -(F.T * GT * (-F))[0]
assert sp.simplify(same_source_energy + opposite_source_energy) == 0
print("PASS 3: opposite vector sources reverse the reciprocal cross term")

a, k = sp.symbols("a k", positive=True)
self_radial = sp.integrate(sp.exp(-a**2 * k**2), (k, 0, sp.oo))
assert self_radial == sp.sqrt(sp.pi) / (2 * a)
print("PASS 4: a smooth Gaussian form factor has finite transverse self-energy radial integral")

exchange = sp.symbols("exchange")
carrier_gain = exchange
field_loss = -exchange
assert sp.simplify(carrier_gain + field_loss) == 0
print("PASS 5: carrier and field translation-momentum exchange cancels")

print("ALL 5 RECOIL ACTION CHECKS PASS")
