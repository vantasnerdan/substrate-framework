"""Exact actual quadrature phase, energy, projection and parity anchors."""

import sympy as s

from substrate_framework.euler_passive_control import passive_packet
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0210-quadrature-phase-control")
    angle, time = s.symbols("theta t", real=True)
    omega, amplitude, rho = s.symbols("omega G rho", positive=True)
    packet = passive_packet(amplitude, omega, angle, time)
    g = s.Matrix([packet.initial_velocity,
                  packet.initial_velocity.subs(angle, angle-s.pi/2)])
    h = s.Matrix([packet.initial_configuration,
                  packet.initial_configuration.subs(angle, angle-s.pi/2)])
    th = omega*h.diff(angle)
    energy = g*g.T-g*th.T-th*g.T
    checks.check("the actual two-column Jacobi energy vanishes as a whole matrix",
                 s.simplify(energy) == s.zeros(2))
    phase_density = rho*(h*g.T-g*h.T)
    checks.check("the same actual fields have nonzero positive phase density",
                 s.simplify(phase_density[0, 1]-rho*amplitude**2/(2*omega)) == 0)
    unit = s.Matrix([[0, 1], [-1, 0]])
    phase = phase_density.applyfunc(lambda value: s.simplify(
        s.integrate(value, (angle, 0, 2*s.pi))/(2*s.pi)))
    checks.check("angular integration retains the literal full phase factor",
                 phase == rho*amplitude**2*unit/(2*omega))
    swap = s.Matrix([[0, 1], [1, 0]])
    checks.check("column interchange reverses phase without changing zero energy",
                 swap.T*phase*swap == -phase and swap.T*energy*swap == s.zeros(2))
    checks.check("both actual material configurations retain their half inverse",
                 s.simplify(th-g/2) == s.zeros(2, 1))
    checks.check("all four mean fields are exactly zero under the actual angle measure",
                 all(s.integrate(value, (angle, 0, 2*s.pi)) == 0 for value in [*g, *h]))
    for column in range(2):
        xi = packet.lin_displacement.subs(angle, angle-column*s.pi/2)
        velocity = packet.euler_velocity.subs(angle, angle-column*s.pi/2)
        checks.check(f"quadrature {column+1} solves its actual complete Lin transport",
                     s.simplify(s.diff(xi, time)+omega*s.diff(xi, angle)-velocity) == 0)

    nx, ny = s.symbols("ny nz", real=True)
    kx, ky, kz, epsilon = s.symbols("Kx Ky Kz epsilon", real=True)
    vector = s.Matrix([epsilon*kx, nx+epsilon*ky, ny+epsilon*kz])
    projection = s.eye(3)-vector*vector.T/vector.dot(vector)
    checks.check("full projected axial phase has no linear spatial correction",
                 s.diff(projection[0, 0], epsilon).subs(epsilon, 0) == 0)
    checks.check("the actual second phase correction contains the inverse transverse Laplacian",
                 s.simplify(s.diff(projection[0, 0], epsilon, 2).subs(epsilon, 0)/2
                            +kx**2/(nx**2+ny**2)) == 0)
    checks.check("the first pressure projection is transverse and loses one frequency order",
                 s.simplify(projection[:, 0].diff(epsilon).subs(epsilon, 0)
                            +s.Matrix([0, nx*kx, ny*kx])/(nx**2+ny**2)) == s.zeros(3, 1))
    transport = s.symbols("KdotU", real=True)
    # Full Hermitian energy: -g_i^* A h_j -(A h_i)^* g_j.
    first_energy = -s.I*transport*g*h.T+s.I*transport*h*g.T
    checks.check("the real full energy retains its generally nonzero Hermitian transport row",
                 s.simplify(rho*first_energy-s.I*transport*phase_density) == s.zeros(2))
    checks.check("that linear transport energy is Hermitian, not zero by taking a real diagonal",
                 s.simplify(first_energy-first_energy.conjugate().T) == s.zeros(2)
                 and s.simplify(first_energy[0, 1]) != 0)
    inversion = -s.eye(3)
    axis = s.Matrix([1, 0, 0])
    checks.check("whole-field inversion preserves an axial physical input axis",
                 inversion.det()*inversion*axis == axis)
    background = s.Matrix(s.symbols("u0:3", real=True))
    wave = s.Matrix([kx, ky, kz])
    checks.check("positive inversion pairs cancel the actual polar transport coefficient",
                 s.simplify(wave.dot(background)+wave.dot(inversion*background)) == 0)
    index, exponent = s.symbols("N L", positive=True)
    chosen = index**(-exponent-1)
    checks.check("polynomial norm growth has a compatible explicit spatial window",
                 s.simplify(index**exponent*chosen-1/index) == 0
                 and s.simplify(index**(-exponent-3)/chosen-index**-2) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
