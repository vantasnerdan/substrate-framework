"""Actual passive Euler sector, phase-energy repair and physical mean chart."""

import sympy as sp

from substrate_framework.euler_phase import physical_scalar_chart
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0196-passive-phase")
    x,y,z,time = sp.symbols("x y z time",real=True)
    a,b = sp.symbols("A B",positive=True)
    psi = b*sp.cos(y)+a*sp.cos(z)
    u = sp.Matrix([psi,a*sp.sin(z),-b*sp.sin(y)])
    scalar = sp.Function("g")(y,z)
    perturbation = sp.Matrix([scalar,0,0])
    coordinates = (x,y,z)
    full_advective_generator = -(perturbation.jacobian(coordinates)*u+u.jacobian(coordinates)*perturbation)
    divergence = sum(sp.diff(full_advective_generator[i],coordinates[i]) for i in range(3))
    transport_g = a*sp.sin(z)*sp.diff(scalar,y)-b*sp.sin(y)*sp.diff(scalar,z)
    checks.check("an arbitrary smooth passive axial field has the full pressure-free Euler generator", full_advective_generator == sp.Matrix([-transport_g,0,0]) and divergence == 0)
    c,delta,local = sp.symbols("c delta local",real=True)
    speed_squared = b*b-(c-a*sp.cos(z))**2
    local_speed = sp.expand(sp.series(speed_squared.subs({c:b-a-delta,z:sp.pi+local}),local,0,4).removeO())
    checks.check("the actual saddle has the positive logarithmic-period quadratic form", sp.expand(local_speed.subs(delta,0)).coeff(local,2) == a*b and sp.diff(local_speed,delta).subs({delta:0,local:0}) == 2*b)
    period_integrand = 1/sp.sqrt(speed_squared)
    checks.check("actual streamline frequency monotonicity uses the differentiated period integrand", sp.simplify(sp.diff(period_integrand,c)-(c-a*sp.cos(z))/speed_squared**sp.Rational(3,2)) == 0)
    theta,omega,G = sp.symbols("theta omega G",real=True,nonzero=True)
    fraction = sp.Symbol("fraction",real=True)
    g = G*sp.cos(theta)
    h = fraction*G*sp.sin(theta)/omega
    th = omega*sp.diff(h,theta)
    material_rate = g-th
    hamiltonian_correction = sp.expand(material_rate**2-th**2)/2
    checks.check("the physical half-inverse configuration cancels the complete VV energy pointwise", sp.simplify(hamiltonian_correction.subs(fraction,sp.Rational(1,2))) == 0)
    checks.check("dropping the half normalization changes the actual energy", sp.simplify(hamiltonian_correction.subs(fraction,1)+g*g/2) == 0)
    checks.check("the actual internal VV phase row vanishes by angular orthogonality", sp.integrate(g*h,(theta,0,2*sp.pi)) == 0)
    s1,mu,w = sp.symbols("s1 mu w",real=True,nonzero=True)
    selected_h = h.subs(fraction,sp.Rational(1,2))
    phase_integrand = sp.integrate(s1*sp.sin(theta)*selected_h,(theta,0,2*sp.pi))/(2*sp.pi)
    checks.check("the complete DV phase carries the negative-frequency quarter moment", sp.simplify(a*mu*phase_integrand.subs(G,10*w/(a*mu*s1))/10-w/(4*omega)) == 0)
    checks.check("the configuration norm retains its actual inverse-frequency growth", sp.simplify(sp.integrate(selected_h**2,(theta,0,2*sp.pi))/(2*sp.pi)-G*G/(8*omega*omega)) == 0)
    packet_output = w*sp.sin(omega*time)
    raw_current = w*sp.cos(omega*time)/(4*omega)+time*packet_output/2
    current_from_output = time*packet_output/2-sp.integrate(packet_output,(time,0,time))/4
    checks.check("the literal packet momentum current retains precisely its initial phase moment", sp.simplify(raw_current-current_from_output-w/(4*omega)) == 0)
    frequencies = sp.symbols("w1 w2 w3",positive=True)
    matrix = sp.Matrix([[q**power for q in frequencies] for power in (-1,1,3)])
    vandermonde = sp.prod(frequencies[j]**2-frequencies[i]**2 for i in range(3) for j in range(i+1,3))/sp.prod(frequencies)
    checks.check("the simultaneous phase-output determinant is the actual generalized Vandermonde", sp.factor(matrix.det()-vandermonde) == 0)
    epsilon = sp.Symbol("epsilon",positive=True)
    # Concrete exact interpolation of t+2t^3, including the extra phase row.
    finite = matrix.subs(dict(zip(frequencies,(epsilon,2*epsilon,3*epsilon),strict=True)))
    weights = finite.inv()*sp.Matrix([0,1,-12])
    response = sum(weights[i]*sp.sin((i+1)*epsilon*time) for i in range(3))
    checks.check("finite actual frequencies solve output and initial phase simultaneously", sp.simplify(sum(weights[i]/((i+1)*epsilon) for i in range(3))) == 0 and sp.expand(sp.series(response,time,0,5).removeO()) == time+2*time**3)
    epsilon_remainder = sp.simplify(sp.series(response,epsilon,0,3).removeO()-(time+2*time**3))
    checks.check("the controlled output remainder is genuinely second order in frequency", epsilon_remainder != 0 and sp.simplify(epsilon_remainder/epsilon**2).has(epsilon) is False)

    rho,stiffness,spatial = sp.symbols("rho stiffness spatial",positive=True)
    remainder = sp.Function("H")(time)
    f = 1-stiffness*spatial*time**2/2
    gg = time-stiffness*spatial*time**3/6+spatial*remainder
    row = sp.Matrix([[f,gg]])
    phase = sp.Matrix([[0,rho],[-rho,0]])
    chart = physical_scalar_chart(phase,sp.zeros(2),row,
                                 angle_rate=row.diff(time),angle_acceleration=row.diff(time,2),
                                 generator_rate=sp.zeros(2),spin=rho*row.diff(time))

    def jet(expression):
        return sp.simplify(sp.series(expression,spatial,0,2).removeO())

    checks.check("the actual physical mean chart retains its controlled Wronskian connection", jet(chart.wronskian-1-spatial*sp.diff(remainder,time)) == 0)
    checks.check("canonical mass and its complete time connection follow the same action", jet(chart.mass-rho*(1-spatial*sp.diff(remainder,time))) == 0 and jet(chart.mass_rate+rho*spatial*sp.diff(remainder,time,2)) == 0)
    checks.check("the actual restoring stiffness is unchanged through the second spatial jet", jet(chart.stiffness-rho*stiffness*spatial) == 0)
    checks.check("measured physical mean momentum is not silently equated to moving canonical mass", sp.simplify(chart.spin_inertia-rho) == 0 and chart.spin_connection == 0 and jet(rho-chart.mass-rho*spatial*sp.diff(remainder,time)) == 0)
    actual_initial_energy = sp.diag(rho*stiffness*spatial,rho)
    inverse = chart.coordinates.inv()
    observed_energy = inverse.T*actual_initial_energy*inverse
    expected_energy = sp.diag(rho*stiffness*spatial,rho*(1-2*spatial*sp.diff(remainder,time)))
    checks.check("the conserved actual Euler energy retains the separate physical-chart connection", all(jet(value) == 0 for value in observed_energy-expected_energy))
    return int(checks.finish())


if __name__ == "__main__":
    raise SystemExit(main())
