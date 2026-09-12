#!/usr/bin/env python
import sympy as sp

from substrate_framework.euler_charge_multipole import (
    internal_stress_source,
    pressure_quadrupole_far_field,
    radial_monopole_flux,
    transverse_static_green,
    transverse_scalar_source,
    vortex_dipole_cross_energy,
    vortex_impulse_far_velocity,
)

x, y, z, rho = sp.symbols("x y z rho", positive=True)
r = sp.sqrt(x**2 + y**2 + z**2)
coords = (x, y, z)
G = 1 / (4 * sp.pi * r)
hessian = sp.Matrix(3, 3, lambda i, j: sp.diff(G, coords[i], coords[j]))
expected_hessian = sp.Matrix(3, 3, lambda i, j: (3 * coords[i] * coords[j] / r**2 - int(i == j)) / (4 * sp.pi * r**3))
assert all(sp.simplify(v) == 0 for v in hessian - expected_hessian)
print("PASS 1: Newton Hessian and physical sign are derived")

I1 = sp.Matrix(sp.symbols("a0:3", real=True))
I2 = sp.Matrix(sp.symbols("b0:3", real=True))
M1 = sp.Matrix(3, 3, lambda i, j: sum(sp.LeviCivita(k, i, j) * I1[k] for k in range(3)))
M2 = sp.Matrix(3, 3, lambda i, j: sum(sp.LeviCivita(k, i, j) * I2[k] for k in range(3)))
contraction = sp.Matrix(3, 3, lambda i, j: sum(M1[i, k] * M2[j, k] for k in range(3)))
expected = sp.eye(3) * I1.dot(I2) - I2 * I1.T
assert contraction == expected
print("PASS 2: antisymmetric vorticity moments contract to the dipole tensor")

axis_energy = vortex_dipole_cross_energy(rho, [0, 0, 1], [0, 0, 1], [0, 0, z])
assert sp.simplify(axis_energy - rho / (2 * sp.pi * z**3)) == 0
assert sp.simplify(-sp.diff(axis_energy, z) - 3 * rho / (2 * sp.pi * z**4)) == 0
print("PASS 3: cross energy is d^-3 and translation force is d^-4")

axis_velocity = vortex_impulse_far_velocity([0, 0, 1], [0, 0, z])
assert axis_velocity == sp.Matrix([0, 0, 1 / (2 * sp.pi * z**3)])
print("PASS 4: impulse velocity has the derived dipole sign")

m1, m2, m3 = sp.symbols("m1 m2 m3", real=True)
pressure = pressure_quadrupole_far_field(rho, sp.diag(m1, m2, m3), [0, 0, z])
assert sp.simplify(pressure - rho * (2 * m3 - m1 - m2) / (4 * sp.pi * z**3)) == 0
assert pressure_quadrupole_far_field(rho, sp.eye(3), [0, 0, z]) == 0
print("PASS 5: pressure begins at the trace-free quadrupole with physical density")

q = sp.symbols("q", real=True)
assert radial_monopole_flux(q) == 4 * sp.pi * q
print("PASS 6: a radial r^-2 monopole carries forbidden nonzero spherical flux")

K1, K2, K3 = sp.symbols("K1 K2 K3", real=True)
Sigma = sp.Matrix(3, 3, sp.symbols("s0:9"))
source = internal_stress_source([K1, K2, K3], Sigma)
assert source.subs({K1: 0, K2: 0, K3: 0}) == sp.zeros(3, 1)
assert source == sp.I * Sigma * sp.Matrix([K1, K2, K3])
G3 = transverse_static_green(2, [0, 0, 3])
G6 = transverse_static_green(2, [0, 0, 6])
assert G6 == G3 / 2
print("PASS 7: internal stress has zero force monopole and transverse Green is r^-1")

c1, c2, c3 = sp.symbols("c1 c2 c3", real=True)
C = sp.Matrix([c1, c2, c3])
Rx = sp.diag(1, -1, -1)
Ry = sp.diag(-1, 1, -1)
Rz = sp.diag(-1, -1, 1)
solutions = sp.solve(list((Rx-sp.eye(3))*C)+list((Ry-sp.eye(3))*C)+list((Rz-sp.eye(3))*C), [c1,c2,c3], dict=True)
assert solutions == [{c1: 0, c2: 0, c3: 0}]
print("PASS 8: isotropy forbids a nonzero scalar-to-transverse-vector source")

factor = sp.symbols("factor", real=True)
assert transverse_scalar_source([K1, K2, K3], factor) == sp.zeros(3, 1)
print("PASS 9: a local isotropic scalar source is annihilated by the transverse projector")

print("ALL 9 EXACT CHARGE-MULTIPOLE CHECKS PASS")
