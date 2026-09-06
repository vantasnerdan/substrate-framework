import sympy as s

from substrate_framework.euler_scale_causality import (
    axisymmetric_swirl_pressure_quadrupole,
    euler_similarity_weights,
    pressure_quadrupole,
)


A, B, rho, m = s.symbols("A B rho m", positive=True)
x, y, z = s.symbols("x y z", positive=True)
w = euler_similarity_weights(A, B)
assert w.vorticity == A * B
assert w.energy == A**2 / B**3
assert w.vorticity_impulse == A / B**3
assert w.kks_action == A / B**4
assert w.helicity == A**2 / B**2
assert w.circulation == A / B
assert w.topological_charge == 1
print("PASS exact fixed-density Euler similarity weights")

target, action = s.symbols("target action", positive=True)
same_topology = euler_similarity_weights(target / action, s.Integer(1))
assert s.simplify(action * same_topology.kks_action-target) == 0
assert same_topology.topological_charge == 1
print("PASS continuous action rescaling inside one topology class")

r2 = x**2+y**2+z**2
green = 1/(4*s.pi*s.sqrt(r2))
M = s.diag(m, m, 0)
direct = rho*sum(M[i,j]*s.diff(green, (x,y,z)[i], (x,y,z)[j])
                 for i in range(3) for j in range(3))
api = pressure_quadrupole(M, s.Matrix([x,y,z]), density=rho)
assert s.simplify(api-direct) == 0
print("PASS positive-inverse Euler pressure sign and Hessian contraction")

axisymmetric = axisymmetric_swirl_pressure_quadrupole(m, x, z, density=rho)
expected = rho*m*(x**2-2*z**2)/(4*s.pi*(x**2+z**2)**s.Rational(5,2))
assert s.simplify(axisymmetric-expected) == 0
axis = s.simplify(axisymmetric.subs(x,0))
assert axis == -rho*m/(2*s.pi*z**3)
assert s.diff(axis,z) == 3*rho*m/(2*s.pi*z**4)
print("PASS compact-axisymmetric swirl pressure and acceleration tails")

R = s.symbols("R", positive=True)
assert s.integrate(1/R, (R, 1, s.oo)) == s.oo
assert s.integrate(1/R**2, (R, 1, s.oo)) == 1
print("PASS r^-3 absolute-sum and r^-4 first-moment logarithmic boundaries")

K, c = s.symbols("K c", positive=True)
phi_t, phi_x = s.symbols("phi_t phi_x", real=True)
lagrangian = K*(phi_t**2-c**2*phi_x**2)/2
assert s.diff(lagrangian, phi_t) == K*phi_t
assert s.diff(lagrangian, phi_x) == -K*c**2*phi_x
principal_ratio = s.simplify(-s.diff(lagrangian, phi_x, phi_x)
                             / s.diff(lagrangian, phi_t, phi_t))
assert principal_ratio == c**2
print("PASS action prefactor cancels from equations while scaling symplectic action")
