"""Measured tag inertia and full KKS normalization, with exposing signs."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0131-measured-moment-normalization")
    omega, rho, mass_fraction = s.symbols("Omega rho mu", positive=True)
    sigma = s.Symbol("sigma", nonzero=True, real=True)
    azimuth = s.Symbol("m", positive=True, integer=True)
    mark, z0, zc, rm, rp, rg = s.symbols("d Z0 Zc Rm Rp Rg", nonzero=True, real=True)
    q0 = rho*mass_fraction*s.pi*mark*rm*z0
    deformation = rho*mass_fraction*azimuth*s.pi*zc*rg/(sigma*(2*omega+sigma))
    spin = rho*mass_fraction*azimuth*s.pi*mark*zc*rp/sigma
    measured = s.simplify(spin*azimuth*q0/(sigma*deformation))
    expected = rho*mass_fraction*azimuth*s.pi*mark**2*z0*rp*rm*(2*omega+sigma)/(sigma*rg)
    checks.check("measured inertia follows from the independent actual shape and spin rows",
                 s.simplify(measured-expected) == 0)
    checks.check("tag fraction scales measured inertia but leaves the angle row unchanged",
                 s.diff(s.simplify(deformation/q0), mass_fraction) == 0
                 and s.simplify(s.diff(measured, mass_fraction)-measured/mass_fraction) == 0)

    delta, radius, value, weight1, weightm = s.symbols("delta r P weight1 weightm", positive=True)
    # At a lower-branch impermeable wall c=1/(1+delta), the combination
    # P'+mP/r is exactly -delta*m*P/r. Continuity extends the sign to
    # a sufficiently thin interior annulus; this is the pointwise limit.
    wall_limit = s.simplify(measured.subs({sigma: -2*omega/(1+delta),
                                         rp: value*weight1,
                                         rg: -delta*azimuth*value*weightm/radius}))
    positive_expression = rho*mass_fraction*s.pi*mark**2*z0*rm*radius*weight1/weightm
    checks.check("wall-adjacent material marker has a positive measured spin-rate coefficient",
                 s.simplify(wall_limit-positive_expression) == 0)

    # Derive the KKS scalar from the complex Euler energy pairing. Pressure
    # contributes no boundary term on the declared divergence-free wall mode.
    re = s.Matrix(s.symbols("a:3", real=True))
    im = s.Matrix(s.symbols("b:3", real=True))
    rotation = s.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
    complex_velocity = re+s.I*im
    pairing = s.expand((s.conjugate(complex_velocity).T*rotation*complex_velocity)[0])
    checks.check("Coriolis pairing fixes the cross-product sign before integration",
                 s.simplify(pairing+2*s.I*re.cross(im)[2]) == 0)
    norm, integrated_cross = s.symbols("norm integrated_cross", real=True)
    solved_cross = s.solve(sigma*norm+4*omega*integrated_cross, integrated_cross)[0]
    beta = s.simplify(2*rho*omega*solved_cross/sigma**2)
    kinetic = rho*norm/2
    checks.check("complete KKS scalar equals minus rotating kinetic Hessian over frequency",
                 s.simplify(beta+kinetic/sigma) == 0)
    checks.check("wrong vorticity factor would halve the KKS normalization",
                 s.simplify(rho*omega*solved_cross/sigma**2+kinetic/sigma) != 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
