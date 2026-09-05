"""Scoped action/scale interface checks; Euler sources are checked separately."""

import sympy as sp

from substrate_framework.euler_phase import physical_scalar_chart
from substrate_framework.homogenization import sphere_fourth_moment_isotropic
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0219-action-window")
    moment_order, n_power, k_power, off_order = 45,4,22,6
    powers = {
        "physical smooth-band error":moment_order+1-2*k_power,
        "full cubic action remainder":k_power-6-3*n_power,
        "control energy second jet":n_power-3,
        "first-gradient physical return":off_order*n_power-k_power-sp.Rational(3,2),
        "control phase second jet":2*n_power-3,
        "macro-scaled prepared velocity H4":k_power-sp.Rational(3,2)-sp.Rational(9,2)*n_power,
    }
    for name, power in powers.items():
        print(f"derived h exponent, {name}: {power}")
        checks.check(f"the SAME explicit diagonal makes {name} vanish", power > 0)
    nu, parallel, perpendicular, inertia = sp.symbols("nu b_parallel b_perp j", positive=True)
    k = sp.Symbol("k",real=True)
    tensor = sphere_fourth_moment_isotropic()
    angular_parallel = sp.Matrix(3,3,lambda i,j:tensor[i,j,2,2])
    angular_perpendicular = sp.eye(3)/3-angular_parallel
    raw = 2*nu*inertia*k**2*(parallel*angular_parallel+perpendicular*angular_perpendicular)
    normalized = sp.simplify(3*raw/inertia)
    expected_t = 2*nu*(parallel+4*perpendicular)/5
    expected_l = 2*nu*(3*parallel+2*perpendicular)/5
    checks.check("the common-laboratory-K Haar contraction gives both physical curvatures", normalized == k**2*sp.diag(expected_t,expected_t,expected_l))
    checks.check("both derived optical curvatures are strictly positive", expected_t.is_positive is True and expected_l.is_positive is True)
    time = sp.Symbol("t",real=True)
    mass, frequency = sp.symbols("M Gamma",positive=True)
    observer = sp.Matrix([[sp.cos(frequency*time),sp.sin(frequency*time)/frequency]])
    phase = mass*sp.Matrix([[0,1],[-1,0]])
    chart = physical_scalar_chart(phase,sp.zeros(2),observer,
        angle_rate=observer.diff(time),angle_acceleration=observer.diff(time,2),
        generator_rate=sp.zeros(2),spin=mass*observer.diff(time))
    checks.check("the canonical physical chart retains the genuinely matched spin inertia", sp.simplify(chart.mass-mass) == 0 and sp.simplify(chart.spin_inertia-mass) == 0)
    checks.check("the limiting prepared action has no concealed time connection", chart.mass_rate == 0 and chart.spin_connection == 0 and sp.simplify(chart.stiffness-mass*frequency**2) == 0)
    gradient_mass = sp.Symbol("J2",real=True)
    actual_mass = inertia+gradient_mass*k**2
    actual_stiffness = actual_mass*(nu**2+expected_t*k**2)
    numerator = sp.series(actual_stiffness-nu**2*actual_mass,k,0,3).removeO()
    checks.check("subtracting the full gradient mass leaves the positive curvature numerator", sp.simplify(numerator-inertia*expected_t*k**2) == 0)
    checks.check("dropping gradient inertia would leave a false stiffness contribution", sp.expand(sp.expand(actual_stiffness-nu**2*inertia).coeff(k,2)-inertia*expected_t-nu**2*gradient_mass) == 0)
    return int(checks.finish())


if __name__ == "__main__":
    raise SystemExit(main())
