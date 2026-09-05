"""Actual full-phase pairing and coherent physical mean action."""

import sympy as s

from substrate_framework.euler_fourier import (
    add, cross, curl, inner, leray, material_kelvin_operator, mul, scale,
    transport, trig,
)
from substrate_framework.euler_phase import moving_phase_pullback, physical_scalar_chart
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0154-coherent-physical-acoustic-action")
    rho = s.Symbol("rho", positive=True)
    k = s.Symbol("k", nonzero=True, real=True)
    plus = (3*k/5, 0, 4*k/5)
    minus = tuple(-value for value in plus)
    profile = {plus: 1/s.sqrt(2), minus: 1/s.sqrt(2)}
    polarization = [s.Rational(4, 5), 0, -s.Rational(3, 5)]
    eta = tuple(scale(profile, value) for value in polarization)
    background = (trig(2), trig(2, kind="sin"), {})
    initial_rate = material_kelvin_operator(background, eta)
    advective = transport(background, eta)
    pi_d = tuple(scale(component, rho) for component in leray(
        tuple(add(initial_rate[i], advective[i]) for i in range(3))))
    pi_v = tuple(scale(component, rho) for component in eta)
    dot = add(*(scale(background[i], polarization[i]) for i in range(3)))
    # The real derivative of exp(+-iKx), not one complex sign reused twice.
    target = tuple(scale(component, rho) for component in leray(tuple(
        mul(dot, {plus: -s.I*plus[i]/s.sqrt(2), minus: s.I*plus[i]/s.sqrt(2)})
        for i in range(3))))
    checks.check("canonical Fourier material momentum retains the finite-K pressure return",
                 all(s.simplify(value) == 0 for i in range(3)
                     for value in add(pi_d[i], scale(target[i], -1)).values()))
    checks.check("the real physical macro column has unit mean-square normalization",
                 s.simplify(inner(eta, eta)) == 1)
    checks.check("full finite-K displacement/common-V symplectic pairing is exactly rho",
                 s.simplify(inner(eta, pi_v)) == rho
                 and s.simplify(inner(eta, pi_d)) == 0)
    euler_d = leray(cross(eta, curl(background)))
    checks.check("material cotangent momentum is not renamed Eulerian momentum",
                 any(s.simplify(value) != 0 for i in range(3)
                     for value in add(pi_d[i], scale(euler_d[i], -rho)).values()))

    # Exact moving pullback identity, independently of any coordinate clock.
    omega = rho*s.Matrix([[0, 1], [-1, 0]])
    h11, h12, h22 = s.symbols("h11 h12 h22", real=True)
    h = s.Matrix([[h11, h12], [h12, h22]])
    embedding = s.Matrix(2, 2, s.symbols("e11 e12 e21 e22", real=True))
    ambient_generator = -omega.inv()*h
    moving = moving_phase_pullback(omega, h, embedding, ambient_generator*embedding)
    checks.check("actual solution-column connection cancels its coefficient Hamiltonian",
                 moving.hamiltonian == s.zeros(2)
                 and moving.generator == s.zeros(2))
    checks.check("the propagated phase embedding has zero ambient residual",
                 moving.residual == s.zeros(2)
                 and moving.symplectic_rate == s.zeros(2))

    t = s.Symbol("t", real=True)
    f, g = s.Function("f", real=True)(t), s.Function("g", real=True)(t)
    row = s.Matrix([[f, g]])
    rate = row.diff(t)
    chart = physical_scalar_chart(
        omega, s.zeros(2), row, angle_rate=rate,
        angle_acceleration=row.diff(t, 2), generator_rate=s.zeros(2), spin=rho*rate)
    wronskian = f*s.diff(g, t)-g*s.diff(f, t)
    numerator = s.diff(f, t)*s.diff(g, t, 2)-s.diff(g, t)*s.diff(f, t, 2)
    checks.check("the physical row derives its exact Wronskian and inertia",
                 s.simplify(chart.wronskian-wronskian) == 0
                 and s.simplify(chart.mass-rho/wronskian) == 0)
    checks.check("the full moving stiffness is derived rather than supplied",
                 s.simplify(chart.stiffness-rho*numerator/wronskian**2) == 0)
    checks.check("the exact physical momentum row and canonical mismatch remain separate",
                 chart.spin_inertia == rho and chart.spin_connection == 0
                 and s.simplify(chart.angle_spin_bracket-wronskian) == 0)
    expected_h = s.Matrix([[chart.stiffness, chart.mass_rate/2],
                           [chart.mass_rate/2, chart.mass]])
    checks.check("the physical Hamiltonian retains both moving connection entries",
                 s.simplify(chart.hamiltonian-expected_h) == s.zeros(2))

    displacement = s.Function("X")(t)
    action = chart.mass*s.diff(displacement, t)**2/2-chart.stiffness*displacement**2/2
    variation = s.diff(s.diff(action, s.diff(displacement, t)), t)-s.diff(action, displacement)
    equation = chart.mass*s.diff(displacement, t, 2) \
        +chart.mass_rate*s.diff(displacement, t)+chart.stiffness*displacement
    checks.check("variation of the actual chart action gives its full physical equation",
                 s.simplify(variation-equation) == 0)

    # 0151's derived rows, used only through their proven spatial order.
    a = s.Symbol("a", positive=True)
    jet_row = s.Matrix([[1-a*k*k*t*t/2, t-a*k*k*t**3/6]])
    jet = physical_scalar_chart(
        omega, s.zeros(2), jet_row, angle_rate=jet_row.diff(t),
        angle_acceleration=jet_row.diff(t, 2), generator_rate=s.zeros(2),
        spin=rho*jet_row.diff(t))

    def second_jet(expression):
        return s.simplify(s.series(expression, k, 0, 3).removeO())

    checks.check("both actual initial phases cancel the second-order Wronskian correction",
                 second_jet(jet.wronskian) == 1 and second_jet(jet.mass) == rho)
    checks.check("the same averaged phase action yields positive acoustic stiffness",
                 second_jet(jet.stiffness) == rho*a*k*k
                 and second_jet(jet.mass_rate) == 0)
    checks.check("physical and canonical momenta agree at the claimed spatial order",
                 second_jet(jet.spin_inertia-jet.mass) == 0
                 and jet.spin_connection == 0)
    wrong_row = s.Matrix([[1-a*k*k*t*t/2, t]])
    wrong_w = s.det(wrong_row.col_join(wrong_row.diff(t)))
    checks.check("omitting the genuine common-V response changes the kinetic coefficient",
                 second_jet(rho/wrong_w) != rho)

    # Exposing algebraic counterexample to averaging after elimination.
    # This is not an Euler oscillator input or a new physical model.
    w = s.Symbol("w", positive=True)
    avg_row = s.Matrix([[(1+s.cos(w*t))/2, (t+s.sin(w*t)/w)/2]])
    avg_w = s.trigsimp(s.det(avg_row.col_join(avg_row.diff(t))))
    target_w = (2+2*s.cos(w*t)+w*t*s.sin(w*t))/4
    checks.check("average physical chart differs from individually averaged unit masses",
                 s.simplify(avg_w-target_w) == 0 and s.simplify(avg_w-1) != 0)
    checks.check("wrong-order reduction can hide a real physical-chart degeneracy",
                 s.simplify(avg_w.subs(t, s.pi/w)) == 0)
    print("Actual finite-K initial symplectic pairing:", s.simplify(inner(eta, pi_v)))
    print("Derived second-jet physical mass:", second_jet(jet.mass))
    print("Derived second-jet physical stiffness:", second_jet(jet.stiffness))
    print("Scope: actual coherent averaged phase action, fixed-time positive second jet;")
    print("not unrestricted autonomous Euler invariance or completed coupled Cosserat closure.")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
