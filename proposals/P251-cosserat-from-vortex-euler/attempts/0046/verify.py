"""Exact field, circulation and intrinsic four-coordinate Routh checks."""

import sympy as s

checks = []


def check(name, condition):
    passed = bool(condition)
    checks.append(passed)
    print(f"{'PASS' if passed else 'FAIL'} {name}")


r, theta, z = s.symbols("r theta z", real=True, positive=True)
lam, amp, anis = s.symbols("lambda A B", positive=True)
f0, f2 = s.Function("f0")(r), s.Function("f2")(r)
psi = amp*f0+anis*f2*s.cos(2*theta)
ode = {s.diff(f0, r, 2): -s.diff(f0, r)/r-lam**2*f0,
       s.diff(f2, r, 2): -s.diff(f2, r)/r+(4/r**2-lam**2)*f2}


def curl(v):
    return s.Matrix([s.diff(v[2], theta)/r-s.diff(v[1], z),
                     s.diff(v[0], z)-s.diff(v[2], r),
                     (s.diff(r*v[1], r)-s.diff(v[0], theta))/r])


u = s.Matrix([s.diff(psi, theta)/r, -s.diff(psi, r), lam*psi])
check("full cylinder Beltrami equation from radial ODE",
      s.simplify((curl(u)-lam*u).subs(ode)) == s.zeros(3, 1))
divergence = s.diff(r*u[0], r)/r+s.diff(u[1], theta)/r+s.diff(u[2], z)
check("full cylinder incompressibility", s.simplify(divergence) == 0)
check("material tube level surfaces", s.simplify(u[0]*s.diff(psi, r)+u[1]*s.diff(psi, theta)/r) == 0)
check("sidewall tangency at the declared J2 root", s.simplify(u[0].subs(f2, 0)) == 0)
actual = amp*s.besselj(0, lam*r)+anis*s.besselj(2, lam*r)*s.cos(2*theta)
jet = amp+lam**2*r**2*(-amp/4+anis*s.cos(2*theta)/8)
check("analytic noncircular core jet derived from Bessel functions",
      s.simplify(s.series(actual, r, 0, 4).removeO()-jet) == 0)
vk = -u.diff(theta)
check("physical global rotation tangent remains Beltrami",
      s.simplify((curl(vk)-lam*vk).subs(ode)) == s.zeros(3, 1))
check("rotation is not an axisymmetric relabeling",
      s.simplify(vk[2]-2*lam*anis*f2*s.sin(2*theta)) == 0 and vk[2] != 0)

# Boundary circulation cancellation is the exact differential of a closed
# one-form; this checks the local integrand before integrating the torus.
vtheta, vz, potential = [s.Function(name)(theta, z) for name in ("vt", "vz", "p")]
wedge = vtheta*s.diff(potential, z)-vz*s.diff(potential, theta)
boundary_div = s.diff(potential*vtheta, z)-s.diff(potential*vz, theta)
check("boundary curl pairing equals divergence for a closed velocity one-form",
      s.expand(wedge-boundary_div-potential*(s.diff(vz, theta)-s.diff(vtheta, z))) == 0)
check("distinct axial harmonics are exactly orthogonal",
      s.integrate(s.cos(3*z)*s.cos(5*z), (z, 0, 2*s.pi)) == 0)
check("mean/common and internal harmonics are exactly orthogonal",
      s.integrate(s.cos(3*z), (z, 0, 2*s.pi)) == 0)

bb, bi, hb, hq, hs, hcross = s.symbols("Bb Bi hb hq hs hcross", real=True, nonzero=True)
bodydot, qdot, q, sb, si = s.symbols("Bdot qdot q sb si", real=True)
lagrangian = bb*sb*bodydot-hb*sb**2/2+bi*si*qdot-(hq*q*q+2*hcross*q*si+hs*si**2)/2
momenta = s.solve([s.diff(lagrangian, sb), s.diff(lagrangian, si)], (sb, si))
reduced = s.expand(lagrangian.subs(momenta))
target = bb**2*bodydot**2/(2*hb)+bi**2*qdot**2/(2*hs)-(hq-hcross**2/hs)*q**2/2
check("full four-coordinate Euler Routh action",
      s.simplify(reduced-target+bi*hcross*q*qdot/hs) == 0)
check("unreduced mixed Hessian changes the locking coefficient",
      s.simplify(s.diff(reduced, q, 2)+hq-hcross**2/hs) == 0)
ib, iq, stiffness = s.symbols("IB IQ K", positive=True)
psidot, betadot = s.symbols("Psidot betadot", real=True)
factor = 1+iq/ib
qd = (psidot-betadot)/factor
bd = betadot+qd
kinetic = s.expand((ib*bd**2+iq*qd**2)/2)
check("physical section map diagonalizes both derived inertias",
      s.simplify(kinetic-(ib**2*psidot**2+ib*iq*betadot**2)/(2*(ib+iq))) == 0)
check("absolute-angle inertia is not copied from relative-angle inertia",
      s.simplify(s.diff(kinetic, psidot, 2)-ib**2/(ib+iq)) == 0)
check("affine-cage optical frequency retains the common rotor",
      s.simplify((stiffness/factor**2)/(ib**2/(ib+iq))-stiffness/(ib+iq)) == 0)
print(f"{sum(checks)}/{len(checks)} checks passed")
raise SystemExit(0 if all(checks) else 1)
