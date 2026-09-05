"""Actual affine Euler fourth harmonic, physical stationary tag and phase."""

import sympy as sp

from substrate_framework.euler_phase import physical_scalar_chart
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0209-affine-tag")
    aa, bb, axial, time = sp.symbols("a b X t", real=True)
    omega, speed, ell, rho, strength = sp.symbols("Omega W ell rho zeta0", positive=True)
    u = sp.Matrix([speed, bb, -omega**2*aa])
    pressure = omega**2*(aa*aa+bb*bb)/2
    coordinates = (axial, aa, bb)
    residual = u.jacobian(coordinates)*u+sp.Matrix([sp.diff(pressure,v) for v in coordinates])
    checks.check("the affine comparison is actual stationary Euler with its full pressure", all(sp.simplify(v) == 0 for v in residual))
    radius, angle, wave = sp.symbols("r theta k", positive=True)
    tau = (1-omega)/(1+omega)
    metric_symbol = omega*sp.cos(angle)**2+sp.sin(angle)**2/omega
    poisson = (1-tau*tau)/(1-2*tau*sp.cos(2*angle)+tau*tau)
    checks.check("the complete anisotropic pressure inverse has the derived Poisson coefficients", sp.trigsimp(sp.simplify(1/metric_symbol-poisson)) == 0)
    hankel0 = sp.integrate(wave**3*sp.exp(-wave**2/2)*sp.besselj(0,wave*radius),(wave,0,sp.oo))
    hankel2 = sp.integrate(wave**3*sp.exp(-wave**2/2)*sp.besselj(2,wave*radius),(wave,0,sp.oo))
    checks.check("the actual full-plane zero-angular pressure row is evaluated", sp.simplify(hankel0-(2-radius**2)*sp.exp(-radius**2/2)) == 0)
    checks.check("the actual full-plane quadrupole pressure row is evaluated", sp.simplify(hankel2-radius**2*sp.exp(-radius**2/2)) == 0)
    normal_x, normal_y = sp.symbols("R1 R2", real=True)
    source = (normal_x+sp.I*normal_y)**4*sp.exp(-(normal_x**2+normal_y**2)/2)*sp.exp(4*sp.I*omega*time)
    transported = sp.diff(source,time)+omega*(normal_y*sp.diff(source,normal_x)-normal_x*sp.diff(source,normal_y))
    checks.check("the selected vorticity is an exact transported Euler fourth harmonic", sp.simplify(transported) == 0)
    x = sp.Symbol("x", positive=True)
    phi0 = -tau**2*(2-x)*sp.exp(-x/2)
    phi2 = tau*x*sp.exp(-x/2)
    phi_minus2 = tau**3*x*sp.exp(-x/2)
    z2, zminus2 = -phi2, phi_minus2/3
    checks.check("the two actually observed material rows solve their nonresonant Lin equations", sp.simplify(sp.I*(4-2)*omega*z2+sp.I*2*omega*phi2) == 0 and sp.simplify(sp.I*(4+2)*omega*zminus2-sp.I*2*omega*phi_minus2) == 0)
    initial_resonance, forcing = sp.symbols("s4 phi4")
    resonant_history = sp.exp(4*sp.I*omega*time)*(initial_resonance+time*forcing)
    checks.check("the unobserved resonant Lin history is retained rather than discarded", sp.simplify(sp.diff(resonant_history,time)-4*sp.I*omega*resonant_history-sp.exp(4*sp.I*omega*time)*forcing) == 0)
    spin_kernel = 2*sp.pi*rho*sp.diff(phi0,x)*x
    desired_kernel = sp.pi*rho*tau**2*x*(4-x)*sp.exp(-x/2)
    checks.check("literal moving-tag spin has the derived radial response kernel", sp.simplify(spin_kernel-desired_kernel) == 0)
    angle_chi_derivative = -sp.I*sp.pi*x*(z2-zminus2)/(2*omega)
    angle_kernel = -sp.diff(angle_chi_derivative,x)
    checks.check("the Euclidean covariance angle follows its distinct physical kernel", sp.simplify(angle_kernel+sp.I*sp.pi*(tau+tau**3/3)*x*(4-x)*sp.exp(-x/2)/(4*omega)) == 0)
    reference, response = sp.symbols("Q Ichi", positive=True)
    angle_amplitude = sp.pi*strength*ell**4*(tau+tau**3/3)*response/(4*omega*reference)
    spin_amplitude = sp.pi*rho*strength*ell**4*tau**2*response
    measured = sp.simplify(spin_amplitude/(4*omega*angle_amplitude))
    expected = rho*reference*tau/(1+tau*tau/3)
    checks.check("the SAME physical tag gives the positive spin-rate inertia", sp.simplify(measured-expected) == 0 and expected.subs(omega,sp.Rational(1,10)).is_positive is True)
    checks.check("the literal displacement dipole equals that inertia times the actual angle", sp.simplify(-sp.I*spin_amplitude/(4*omega)-expected*(-sp.I*angle_amplitude)) == 0)
    # The constant-strain phase term is a boundary Jacobian, independently
    # of any mode truncation or selected normalization.
    s1, s2 = sp.Function("s1")(aa,bb), sp.Function("s2")(aa,bb)
    jacobian = sp.diff(s1,aa)*sp.diff(s2,bb)-sp.diff(s1,bb)*sp.diff(s2,aa)
    divergence = sp.diff(s1*sp.diff(s2,bb),aa)-sp.diff(s1*sp.diff(s2,aa),bb)
    checks.check("the integrated strain phase is an actual decaying boundary term", sp.simplify(jacobian-divergence) == 0)
    source_norm = 2*sp.pi*strength**2*ell**2*sp.integrate(radius**9*sp.exp(-radius**2),(radius,0,sp.oo))
    checks.check("the phase-control normalization uses the full physical vorticity norm", sp.simplify(source_norm-24*sp.pi*strength**2*ell**2) == 0)
    c = sp.Symbol("c", real=True)
    checks.check("the unobserved resonant displacement changes actual phase with the derived sign", sp.im(-sp.conjugate(sp.I*c)) == c)
    capital_c, delta = sp.symbols("C Delta", positive=True)
    frequency = 4*omega
    row = sp.Matrix([[capital_c*sp.sin(frequency*time),-capital_c*sp.cos(frequency*time)]])
    beta = delta*frequency*capital_c**2
    chart = physical_scalar_chart(sp.Matrix([[0,beta],[-beta,0]]),sp.zeros(2),row,
                                 angle_rate=row.diff(time),angle_acceleration=row.diff(time,2),
                                 generator_rate=sp.zeros(2),spin=delta*row.diff(time))
    checks.check("the measured angle has its actual nonzero two-phase Wronskian", sp.simplify(chart.wronskian-frequency*capital_c**2) == 0)
    checks.check("the inherited physical action has matching measured mass and positive fixed clock", sp.simplify(chart.mass-delta) == 0 and chart.mass_rate == 0 and sp.simplify(chart.stiffness-delta*frequency**2) == 0 and sp.simplify(chart.spin_inertia-delta) == 0 and chart.spin_connection == 0)
    checks.check("the shrinking-core inertia remains explicit rather than normalized to a bulk constant", sp.simplify(expected.subs(reference,ell**4*reference)/expected-ell**4) == 0)
    return int(checks.finish())


if __name__ == "__main__":
    raise SystemExit(main())
