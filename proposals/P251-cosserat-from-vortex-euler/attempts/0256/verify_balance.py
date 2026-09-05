"""Exact corroboration of0256 identities; not a free-boundary solver."""

import sympy as s

eps, x, phi_x, g = s.symbols("eps x phi_x g", real=True)
curvature = eps / (1 + eps*x) * phi_x + ((1 + eps*x)**2 - 1)*g
forcing = s.diff(curvature, eps).subs(eps, 0)
assert s.expand(forcing - phi_x - 2*x*g) == 0
assert s.expand(forcing - phi_x - x*g) != 0

# An independent compact C2 radial diagnostic has all integration-by-parts
# boundary terms zero. It is not claimed to have the proposed Bessel core or
# the flat smooth physical edge. Its polynomial integrals expose every factor.
rho, mu = s.symbols("rho mu", positive=True)
profile = (1-rho**2)**3
derivative = s.diff(profile, rho)
g_profile = mu*profile + s.diff(profile, rho, 2) + derivative/rho
G_profile = s.integrate(s.expand(g_profile*derivative), rho)
G_profile -= G_profile.subs(rho, 1)
mass = 2*s.pi*s.integrate(rho*profile**2, (rho, 0, 1))
x_gradient = s.pi*s.integrate(rho*derivative**2, (rho, 0, 1))
G_integral = 2*s.pi*s.integrate(rho*G_profile, (rho, 0, 1))
projection_direct = s.pi*s.integrate(
    rho*derivative**2 + 2*rho**2*g_profile*derivative, (rho, 0, 1)
)
assert s.simplify(G_integral-mu*mass/2) == 0
assert s.simplify(projection_direct-(x_gradient-mu*mass)) == 0
balance = s.solve(projection_direct, mu)[0]
assert s.simplify(projection_direct.subs(mu, balance)) == 0
assert projection_direct.subs(mu, balance/2) > 0
assert projection_direct.subs(mu, 2*balance) < 0
weighted_g = 2*s.pi*s.integrate(rho*g_profile*profile, (rho, 0, 1))
assert s.simplify((weighted_g+mu*mass).subs(mu, balance)) == 0

# Scaling derives the finite-radius identities independently of the radial
# diagnostic, retaining the negative swirl term EF from the functional.
k = s.symbols("k", positive=True)
Er, Ez, EF, EB = s.symbols("Er Ez EF EB", real=True)
radial_row = s.diff(Er/k**2+Ez+EF+k**2*EB, k).subs(k, 1)
axial_row = s.diff(k*Er+Ez/k+k*EF+k*EB, k).subs(k, 1)
virial = s.solve([radial_row, axial_row], [EB, Ez])
assert virial[EB] == Er
assert virial[Ez] == 2*Er+EF
assert s.simplify((Ez-EF-2*Er).subs(virial)) == 0
assert s.simplify((Ez+EF-2*Er).subs(virial)) != 0

print("first curvature forcing:", forcing)
print("diagnostic mass, x-gradient:", mass, x_gradient)
print("direct translation pairing:", s.factor(projection_direct))
print("balanced diagnostic lambda^2:", balance)
print("finite-radius virial solution:", virial)
print("PASS: exact projection, Pohozaev, signed-transition and virial identities")
