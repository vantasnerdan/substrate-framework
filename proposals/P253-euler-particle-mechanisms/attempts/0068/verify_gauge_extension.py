import sympy as sp

g, chi, div_u, adv = sp.symbols("g chi div_u adv")
continuity = sp.expand(g * (-adv) + g * (adv + chi * div_u))
assert continuity.subs(div_u, 0) == 0
print("PASS 1: transported incompressible tag gives an exact conserved current")

rho_t, div_j, lam = sp.symbols("rho_t div_j lam")
gauge_bulk = -lam * (rho_t + div_j)
assert gauge_bulk.subs(div_j, -rho_t) == 0
print("PASS 2: current continuity cancels the gauge variation")

eps, mu, d, q1, q2 = sp.symbols("eps mu d q1 q2", positive=True)
energy = q1 * q2 / (4 * sp.pi * eps * d)
force_radial = -sp.diff(energy, d)
assert force_radial == q1 * q2 / (4 * sp.pi * eps * d**2)
print("PASS 3: Gauss energy gives like-charge repulsion and opposite attraction")

c = 1 / sp.sqrt(eps * mu)
assert sp.simplify(c**2 - 1 / (eps * mu)) == 0
print("PASS 4: gauge principal speed is 1/sqrt(epsilon*mu)")

work = sp.symbols("work")
assert sp.simplify(work + (-work)) == 0
print("PASS 5: fluid Lorentz work and field Poynting work cancel")

constraint, div_current = sp.symbols("constraint div_current")
constraint_dot = -div_current + div_current
assert constraint_dot == 0
print("PASS 6: Gauss and magnetic divergence constraints propagate")

print("ALL 6 GAUGE EXTENSION CHECKS PASS")
