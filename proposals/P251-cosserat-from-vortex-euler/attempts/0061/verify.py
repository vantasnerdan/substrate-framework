"""Exact relative-frame kinetic action, with its nontrivial curvature correction."""

import sympy as s

from substrate_framework.euler_orbit import micropolar_kinetic_normal_form
from substrate_framework.verification import CheckLedger

ledger = CheckLedger("P251/0061")
rho, j, alpha, c, mu = s.symbols("rho j alpha C mu", positive=True)
k, ud, pd = s.symbols("k ud pd", real=True)
for helicity in (-1, 1):
    relative_rate = pd-helicity*k*ud/2
    lagrangian_kinetic = rho*ud**2/2+j*relative_rate**2/2
    mass = s.hessian(lagrangian_kinetic, [ud, pd])
    ledger.check(f"h={helicity}: complete relative-frame mass determinant is rho*j",
                 s.simplify(mass.det()-rho*j) == 0 and mass[0, 0].is_positive)
    result = micropolar_kinetic_normal_form(rho, j, alpha, c, j/4, 0, -j/2, k, helicity)
    ledger.check(f"h={helicity}: kinetic cross produces C_eff=C-alpha*j/rho",
                 s.simplify(result.transverse_curvature-c+alpha*j/rho) == 0)
    potential = s.Matrix([[(mu+alpha)*k*k, -2*alpha*helicity*k],
                          [-2*alpha*helicity*k, 4*alpha+c*k*k]])
    pulled = result.field_map.T*potential*result.field_map
    ledger.check(f"h={helicity}: actual potential pullback retains the same correction",
                 s.simplify(s.expand(pulled[1, 1]).coeff(k, 2)-c+alpha*j/rho) == 0)
    # Exact U,q decoupling is not new physical coupling by a variable rename.
    coordinate = s.Matrix([[1, 0], [-helicity*k/2, 1]])
    exact_mass = coordinate.T*s.diag(rho, j)*coordinate
    exact_potential = coordinate.T*s.diag(mu*k*k, 4*alpha+c*k*k)*coordinate
    omega2 = s.Symbol("omega2", real=True)
    determinant = s.det(exact_potential-omega2*exact_mass)
    expected = (mu*k*k-rho*omega2)*(4*alpha+c*k*k-j*omega2)
    ledger.check(f"h={helicity}: exact separable-action modes remain separable in U,q",
                 s.simplify(determinant-expected) == 0)

ledger.check("dropping the relative kinetic cross changes the curvature coefficient",
             s.simplify(c-(c-alpha*j/rho)) != 0)
raise SystemExit(ledger.finish())
