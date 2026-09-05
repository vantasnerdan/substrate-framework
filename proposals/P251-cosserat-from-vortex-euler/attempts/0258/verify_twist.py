"""Exact Bessel series and physical flux factors for the0258 core twist."""

import sympy as s

rho, lam, c, R = s.symbols("rho lam c R", positive=True)
j0 = s.besselj(0, lam*rho).series(rho, 0, 8).removeO()
j1 = s.besselj(1, lam*rho).series(rho, 0, 9).removeO()
q = (j1/(rho*j0)).series(rho, 0, 6).removeO()
flux_action = c*rho*j1
direct_flux = lam*c*s.integrate(j0*rho, (rho, 0, rho))
assert s.series(direct_flux-flux_action, rho, 0, 8).removeO() == 0
twist = s.limit(s.diff(q, rho)/s.diff(flux_action, rho), rho, 0)
assert twist == lam**2/(8*c)
vorticity_twist = s.limit(s.diff(q, rho)/s.diff(lam*flux_action, rho), rho, 0)
assert vorticity_twist == lam/(8*c)
assert s.simplify(twist-vorticity_twist) != 0
return_angle_twist = 2*s.pi*R*twist
assert return_angle_twist == s.pi*R*lam**2/(4*c)

# Substitution into the actual finite-ring period/flux quadrature's straight
# limit, independently retaining its toroidal angular-velocity division by r.
phi = c*j0
speed = c*lam*j1
delta_phi = lam*phi*(2*s.pi*rho)/(R*speed)
normalized_rotation = s.cancel(2*s.pi/(R*delta_phi))
assert s.cancel(normalized_rotation-j1/(rho*j0)) == 0
wrong_delta_phi = lam*phi*(2*s.pi*rho)/speed
wrong_normalized = s.cancel(2*s.pi/(R*wrong_delta_phi))
assert s.cancel(wrong_normalized-normalized_rotation) != 0

print("normalized rotation series:", q)
print("velocity flux action series:", s.series(flux_action, rho, 0, 6))
print("dq/dJ_velocity at core:", twist)
print("dq/dJ_vorticity at core:", vorticity_twist)
print("d(return angle)/dJ_velocity at core:", return_angle_twist)
print("PASS: exact Bessel coefficient, physical flux and toroidal period factors")
