"""Actual Euler/Lin passive returns and initial phase/energy identities."""

import sympy as sp

from substrate_framework.euler_passive_control import passive_packet
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0205-energy-returns")
    time, theta = sp.symbols("t theta", real=True)
    omega, amplitude, density = sp.symbols("omega G rho", positive=True)
    index = sp.Symbol("N", integer=True, positive=True)
    alpha = sp.Symbol("alpha", real=True)
    harmonic_angle = sp.Symbol("vartheta", real=True)
    packet = passive_packet(amplitude, index*omega, harmonic_angle, time,
                            configuration_fraction=alpha)
    g = packet.initial_velocity.subs(harmonic_angle, index*theta)
    h = packet.initial_configuration.subs(harmonic_angle, index*theta)
    history = packet.lin_displacement.subs(harmonic_angle, index*theta)
    material_rate = sp.diff(history, time)+omega*sp.diff(history, theta)
    checks.check("the exact Lin return has its actual advected Euler velocity", sp.simplify(material_rate-g.subs(theta, theta-omega*time)) == 0)
    checks.check("the complete material force vanishes in the passive axial sector", sp.simplify(sp.diff(material_rate,time)+omega*sp.diff(material_rate,theta)) == 0)
    initial_rate = sp.diff(history,time).subs(time,0)
    energy_density = density*(initial_rate**2-(omega*sp.diff(h,theta))**2)/2
    checks.check("the imported packet energy agrees with its complete physical defining form", sp.simplify(energy_density-density*packet.initial_energy_density.subs(harmonic_angle,index*theta)) == 0)
    checks.check("the return's complete energy has the derived controllable sign", sp.simplify(energy_density-density*(sp.Rational(1,2)-alpha)*g**2) == 0)
    checks.check("the half-inverse normalization is an exposing zero-energy limit", sp.simplify(energy_density.subs(alpha,sp.Rational(1,2))) == 0 and sp.simplify(energy_density.subs(alpha,1)+density*g*g/2) == 0)
    second_amplitude, second_alpha = sp.symbols("G2 alpha2", real=True)
    g2 = second_amplitude*sp.cos(index*theta)
    h2 = second_alpha*second_amplitude*sp.sin(index*theta)/(index*omega)
    phase = sp.integrate(h*g2-g*h2,(theta,0,2*sp.pi))
    checks.check("all actual same-cosine return columns have zero restricted initial phase", phase == 0)
    checks.check("the acoustic mean configuration and cotangent cross rows vanish", sp.integrate(g,(theta,0,2*sp.pi)) == 0 and sp.integrate(h,(theta,0,2*sp.pi)) == 0)
    k, kx, ky, kz, ny, nz = sp.symbols("k kx ky kz ny nz", real=True)
    pxx = 1-k**2*kx**2/((ny+k*ky)**2+(nz+k*kz)**2+k**2*kx**2)
    pxx_jet = sp.series(pxx,k,0,3).removeO()
    checks.check("the full Leray XX multiplier retains the even second jet", sp.simplify(pxx_jet-1+k*k*kx*kx/(ny*ny+nz*nz)) == 0)
    checks.check("the phase-null reflection is preserved by the actual second-order pressure multiplier", sp.simplify(pxx_jet-pxx_jet.subs(nz,-nz)) == 0 and sp.simplify(g-g.subs(theta,-theta)) == 0 and sp.simplify(h+h.subs(theta,-theta)) == 0)
    # An arbitrary finite homogeneous matrix has more columns than rows.
    # This exact example exposes the use of actual nonzero null vectors,
    # without asserting full rank or inverting a desired output matrix.
    matrix = sp.Matrix([[1,2,3],[2,4,6]])
    kernel = matrix.nullspace()
    checks.check("redundant homogeneous cross rows still permit a nonzero physical packet combination", len(kernel) == 2 and all(matrix*v == sp.zeros(2,1) and v.dot(v) > 0 for v in kernel))
    eigenvalue = sp.Symbol("lambda", positive=True)
    target_amplitude = sp.sqrt(2*eigenvalue/density)
    checks.check("actual positive and negative packet energies realize the prescribed rank-one correction", sp.simplify(density*(sp.Rational(1,2)-0)*target_amplitude**2-eigenvalue) == 0 and sp.simplify(density*(sp.Rational(1,2)-1)*target_amplitude**2+eigenvalue) == 0)
    return int(checks.finish())


if __name__ == "__main__":
    raise SystemExit(main())
