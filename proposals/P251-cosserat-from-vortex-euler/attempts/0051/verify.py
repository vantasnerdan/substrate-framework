"""Exact weak pressure-bond balances, with the angular boundary term."""

import sympy as s

from substrate_framework.verification import CheckLedger


ledger = CheckLedger("P251-0051-material-pressure-spin")
t = s.symbols("t", real=True)
force = s.Matrix(s.symbols("F1:4", real=True))
branch = s.Matrix(s.symbols("R1:4", real=True))
moment = s.Matrix(s.symbols("m1:4", real=True))
couple = branch.cross(force)
other_moment = -moment+couple
coefficients = s.symbols("a0:6", real=True)
probe = sum(value*t**index for index, value in enumerate(coefficients))
line_probe = s.integrate(probe, (t, 0, 1))
weak_force = -force*s.integrate(s.diff(probe, t), (t, 0, 1))
expected_force = force*(probe.subs(t, 0)-probe.subs(t, 1))
ledger.check("arbitrary degree-five test trace: exact bond force divergence",
             s.simplify(weak_force-expected_force) == s.zeros(3, 1))
weak_couple = -(moment*s.integrate(s.diff(probe, t), (t, 0, 1))
                -couple*s.integrate(t*s.diff(probe, t), (t, 0, 1)))
expected_torques = moment*probe.subs(t, 0)+other_moment*probe.subs(t, 1)
ledger.check("full distributed couple includes center-force moment",
             s.simplify(weak_couple+couple*line_probe-expected_torques) == s.zeros(3, 1))
sigma = force*branch.T
axial = s.Matrix([sum(s.LeviCivita(i, j, k)*sigma[j, k]
                     for j in range(3) for k in range(3)) for i in range(3)])
ledger.check("axial-stress convention fixed from tensor components", axial == -couple)
ledger.check("wrong angular reaction sign is exposed",
             s.simplify(weak_couple-couple*line_probe-expected_torques) != s.zeros(3, 1))
ledger.check("premature opposite-torque assignment loses noncentral force moment",
             other_moment+moment == couple and couple != s.zeros(3, 1))

x, y, z = s.symbols("x y z", real=True)
coordinates = (x, y, z)
position = s.Matrix(coordinates)
stress = s.Matrix(3, 3, lambda i, j: s.Function(f"s{i}{j}")(*coordinates))
div_stress = s.Matrix([sum(s.diff(stress[i, j], coordinates[j])
                          for j in range(3)) for i in range(3)])
angular_flux = s.Matrix.hstack(*(position.cross(stress[:, j]) for j in range(3)))
div_angular = s.Matrix([sum(s.diff(angular_flux[i, j], coordinates[j])
                           for j in range(3)) for i in range(3)])
ax_stress = s.Matrix([sum(s.LeviCivita(i, j, k)*stress[j, k]
                         for j in range(3) for k in range(3)) for i in range(3)])
ledger.check("orbital and spin stress torques cancel in total angular balance",
             s.simplify(position.cross(div_stress)-div_angular-ax_stress) == s.zeros(3, 1))

velocity = s.Matrix([s.Function(f"u{i}")(*coordinates) for i in range(3)])


def curl(field):
    return s.Matrix([s.diff(field[2], y)-s.diff(field[1], z),
                     s.diff(field[0], z)-s.diff(field[2], x),
                     s.diff(field[1], x)-s.diff(field[0], y)])


ledger.check("vorticity impulse differs from physical spin by exact boundary flux",
             s.simplify(curl(position.dot(position)*velocity)
                        -2*position.cross(velocity)-position.dot(position)*curl(velocity)) == s.zeros(3, 1))
# Solid rotation in a unit ball: both volume moments are exact radial integrals.
omega, rho = s.symbols("Omega rho", positive=True)
u = s.Matrix([-omega*y, omega*x, 0])
r = s.symbols("r", positive=True)
second_radial_moment = 4*s.pi*s.integrate(r**4, (r, 0, 1))
spin = rho*omega*2*second_radial_moment/3
impulse = -rho*omega*second_radial_moment
boundary = rho*omega*(2*s.pi-s.Rational(2, 3)*s.pi)*1
ledger.check("solid rotation exposes a nonzero angular boundary contribution",
             s.simplify(spin-impulse-boundary) == 0 and boundary != 0)
ledger.check("using vorticity impulse alone even gives the wrong solid-rotation sign",
             spin.is_positive and impulse.is_negative)
raise SystemExit(ledger.finish())
