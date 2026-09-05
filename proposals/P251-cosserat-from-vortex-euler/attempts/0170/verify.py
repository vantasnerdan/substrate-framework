"""Exact pressure/current cancellation and physical action checks; no truncation."""

from __future__ import annotations

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0170")
    k, q = sp.symbols("k q", real=True, nonzero=True)
    kap = sp.Matrix([0, 1, 0])
    wave = sp.Matrix([q, 0, 0])
    shifted = wave + k * kap
    projector = sp.eye(3) - shifted * shifted.T / shifted.dot(shifted)
    gradient_return = projector * (sp.I * wave)
    ledger.check("full shifted Leray gradient cancellation", all(sp.simplify(entry) == 0 for entry in gradient_return + sp.I * k * projector * kap))
    derivative = projector.diff(k).subs(k, 0)
    ledger.check("bare projector derivative has a real inverse-gap singularity", derivative[0, 1] == -1 / q)
    hessian_column = -wave * wave[0]
    ledger.check("actual Hessian-source derivative retains a small-wave numerator", derivative * hessian_column == q * kap)
    harmonic_projector = sp.eye(3) - kap * kap.T
    ledger.check("harmonic gradient identity uses the slow transverse projection", harmonic_projector * kap == sp.zeros(3, 1))

    u = (
        ef.add(ef.trig(2, kind="sin"), ef.scale(ef.trig(1), 3)),
        ef.add(ef.scale(ef.trig(0, kind="sin"), 2), ef.trig(2)),
        ef.add(ef.scale(ef.trig(1, kind="sin"), 3), ef.scale(ef.trig(0), 2)),
    )
    pressure = ef.scale(ef.add(*(ef.mul(component, component) for component in u)), -sp.Rational(1, 2))
    grad = tuple(ef.derivative(pressure, j) for j in range(3))
    acceleration = ef.transport(u, u)
    ledger.check("nonconstant-pressure ABC test field solves stationary Euler", all(not ef.add(acceleration[j], grad[j]) for j in range(3)))
    ledger.check("test pressure is not the elementary-wave constant", any(grad))
    chi = ef.leray((grad[1], grad[2], {}))
    # A genuinely transverse real test field, independent of the pressure data.
    chi_t = ef.add(ef.trig(1), ef.trig(2)), ef.trig(2), ef.trig(0)
    ledger.check("first-cell test data are solenoidal and mean zero", not ef.divergence(chi) and not ef.divergence(chi_t) and all(component.get(ef.ZERO, 0) == 0 for component in (*chi, *chi_t)))
    hess_chi = tuple(ef.add(*(ef.mul(ef.derivative(grad[i], j), chi[j]) for j in range(3))) for i in range(3))
    a_chi_t = ef.transport(u, chi_t)
    a2_chi = ef.transport(u, ef.transport(u, chi))
    chi_tt = ef.leray(tuple(ef.add(ef.scale(a_chi_t[i], -2), ef.scale(a2_chi[i], -1), ef.scale(hess_chi[i], -1)) for i in range(3)))
    bracket = tuple(ef.add(left, ef.scale(right, -1)) for left, right in zip(ef.transport(u, chi), ef.transport(chi, u), strict=True))
    w = tuple(ef.add(chi_t[i], bracket[i]) for i in range(3))
    w_t = tuple(ef.add(chi_tt[i], a_chi_t[i], ef.scale(ef.transport(chi_t, u)[i], -1)) for i in range(3))
    euler_rhs = ef.leray(tuple(ef.scale(ef.add(left, right), -1) for left, right in zip(ef.transport(u, w), ef.transport(w, u), strict=True)))
    ledger.check("full Jacobi equation factors into the actual first-order Euler/Lin system", all(not ef.add(w_t[i], ef.scale(euler_rhs[i], -1)) for i in range(3)))

    # Build an actual second cell satisfying div(zeta)=kappa.chi.
    zeta = tuple({mode: -sp.I * mode[j] * coefficient / sum(entry**2 for entry in mode) for mode, coefficient in chi[1].items() if mode != ef.ZERO} for j in range(3))
    ledger.check("constructed second cell satisfies the Bloch divergence jet", not ef.add(ef.divergence(zeta), ef.scale(chi[1], -1)))
    hessian_mean = sp.Matrix([sum(ef.mul(ef.derivative(grad[i], j), zeta[j]).get(ef.ZERO, 0) for j in range(3)) for i in range(3)])
    pressure_return = sp.Matrix([ef.mul(grad[i], chi[1]).get(ef.ZERO, 0) for i in range(3)])
    ledger.check("unknown second cell drops out through the full pressure identity", hessian_mean == -pressure_return)
    a = u[1]
    a_chi = ef.transport(u, chi)
    convective_mean = sp.Matrix([ef.mul(a, component).get(ef.ZERO, 0) for component in a_chi])
    first_pressure = sp.Matrix([ef.mul(grad[1], component).get(ef.ZERO, 0) for component in chi])
    ledger.check("transport integration retains the other pressure row", convective_mean == first_pressure)
    ledger.check("dropping pressure rows is exposed by the actual test field", harmonic_projector * (first_pressure + pressure_return) != sp.zeros(3, 1))

    # Derive the mean material acceleration directly from all unprojected terms.
    a2_zeta = ef.transport(u, ef.transport(u, zeta))
    aa = ef.transport(u, ({}, a, {}))[1]
    raw_second = sp.Matrix([
        -a2_zeta[i].get(ef.ZERO, 0)
        -hessian_mean[i]
        +2 * ef.mul(a, chi_t[i]).get(ef.ZERO, 0)
        +ef.mul(aa, chi[i]).get(ef.ZERO, 0)
        +2 * ef.mul(a, a_chi[i]).get(ef.ZERO, 0)
        for i in range(3)
    ])
    predicted_material = sp.Matrix([2 * ef.mul(a, chi_t[i]).get(ef.ZERO, 0) for i in range(3)]) + first_pressure + pressure_return
    ledger.check("full material mean acceleration derives both pressure signs", raw_second == predicted_material)
    physical_current_t = sp.Matrix([ef.add(ef.mul(a, chi_t[i]), ef.scale(ef.mul(u[i], chi_t[1]), -1)).get(ef.ZERO, 0) for i in range(3)])
    physical_stress = sp.Matrix([ef.add(ef.mul(a, chi_t[i]), ef.mul(u[i], chi_t[1])).get(ef.ZERO, 0) for i in range(3)])
    ledger.check("actual mean-current subtraction changes the factor-two stress", raw_second - physical_current_t == physical_stress + first_pressure + pressure_return)
    ledger.check("substituting material mean for physical mean gives a nonzero error", harmonic_projector * physical_current_t != sp.zeros(3, 1))

    t, rho, epsilon, speed = sp.symbols("t rho epsilon a0", positive=True)
    # A nonautonomous perturbation makes the required connections observable.
    f = 1 + k**2 * (-speed * t**2 / 2 + epsilon * t**3)
    g = t + k**2 * (-speed * t**3 / 6 + epsilon * t**4)
    wronskian = sp.expand(f * sp.diff(g, t) - g * sp.diff(f, t))
    mass = rho / wronskian
    stiffness = rho * (sp.diff(f, t) * sp.diff(g, t, 2) - sp.diff(g, t) * sp.diff(f, t, 2)) / wronskian**2
    mass_jet = sp.series(mass / rho, k, 0, 3).removeO()
    stiffness_jet = sp.series(stiffness / rho, k, 0, 3).removeO()
    ledger.check("physical mass carries the true perturbed Wronskian", sp.expand(mass_jet - (1 - 2 * epsilon * k**2 * t**3)) == 0)
    ledger.check("physical stiffness retains its nonautonomous second jet", sp.expand(stiffness_jet - k**2 * (speed - 6 * epsilon * t)) == 0)
    residual = sp.diff(mass, t) * sp.diff(f, t) + mass * sp.diff(f, t, 2) + stiffness * f
    ledger.check("Wronskian action reproduces the actual observed history", sp.simplify(residual) == 0)
    ledger.check("discarding the time-dependent mass connection changes the equation", sp.simplify(mass * sp.diff(f, t, 2) + stiffness * f) != 0)
    print(f"pressure rows: {list(first_pressure)}, {list(pressure_return)}")
    print(f"material/physical current difference: {list(physical_current_t)}")
    print(f"mass/rho second jet: {mass_jet}; stiffness/rho: {stiffness_jet}")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
