"""Exact full-Euler cell action and its pressure/divergence correction."""

import sympy as s

from substrate_framework import euler_fourier as f
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0144-full-cell-Jacobi")
    psi = f.mul(f.trig(0), f.trig(1))
    base = (f.scale(f.derivative(psi, 1), -1), f.derivative(psi, 0),
            f.scale(psi, s.sqrt(2)))
    pressure = f.scale(f.add(*(f.mul(x, x) for x in base)), -s.Rational(1, 2))
    acceleration = f.transport(base, base)
    ledger.check("comparison field is actual stationary incompressible Euler",
                 not f.divergence(base) and all(not f.add(acceleration[j], f.derivative(pressure, j))
                                               for j in range(3)))
    omega = f.curl(base)
    ledger.check("comparison is an actual columnar constant-curl field, not a supplied oscillator",
                 all(not f.add(omega[j], f.scale(base[j], s.sqrt(2))) for j in range(3)))

    time, eps = s.symbols("t epsilon", real=True)
    q = s.Function("Q", real=True)(time)
    r = s.Function("R", real=True)(time)
    direction = (1, 0, 1)
    mean = ({}, {f.ZERO: q}, {})
    stream = f.mul(f.trig(0, kind="sin"), f.trig(1))
    chi = (f.scale(f.derivative(stream, 1), -r),
           f.scale(f.derivative(stream, 0), r),
           f.scale(f.trig(1, 2, "sin"), r))
    longitudinal = f.add(*(f.scale(chi[j], direction[j]) for j in range(3)))
    potential = {g: -value/s.sympify(sum(x*x for x in g))
                 for g, value in longitudinal.items() if g != f.ZERO}
    second = tuple(f.derivative(potential, j) for j in range(3))
    ledger.check("second Bloch divergence correction is derived by periodic inverse Laplacian",
                 not f.divergence(chi) and not f.add(f.divergence(second), f.scale(longitudinal, -1)))

    hessian = [[f.derivative(f.derivative(pressure, j), i) for j in range(3)] for i in range(3)]
    hsecond = tuple(f.add(*(f.mul(hessian[i][j], second[j]) for j in range(3))) for i in range(3))
    ugp = f.mul(mean[1], f.derivative(pressure, 1))
    pressure_border = f.inner(mean, hsecond)
    expected_border = -f.mul(ugp, longitudinal).get(f.ZERO, 0)
    ledger.check("actual pressure border has the required nonzero sign and magnitude",
                 s.simplify(pressure_border-expected_border) == 0 and s.simplify(pressure_border) != 0)

    eta = tuple(f.add(mean[j], f.scale(chi[j], s.I*eps), f.scale(second[j], eps**2))
                for j in range(3))
    adot = f.add(*(f.scale(base[j], direction[j]) for j in range(3)))
    material_velocity = tuple(f.add({g: s.diff(c, time) for g, c in eta[j].items()},
                                   f.transport(base, eta)[j], f.scale(f.mul(adot, eta[j]), s.I*eps))
                              for j in range(3))
    def conjugate(vector):
        real_derivatives = {s.conjugate(s.diff(value, time)): s.diff(value, time)
                            for value in (q, r)}
        return tuple({tuple(-x for x in g): s.conjugate(c).xreplace(real_derivatives)
                      for g, c in item.items()} for item in vector)
    heta = tuple(f.add(*(f.mul(hessian[i][j], eta[j]) for j in range(3))) for i in range(3))
    direct = s.expand((f.inner(conjugate(material_velocity), material_velocity)
                       - f.inner(conjugate(eta), heta))/2).coeff(eps, 2)
    dchi = tuple(f.add({g: s.diff(c, time) for g, c in chi[j].items()}, f.transport(base, chi)[j])
                 for j in range(3))
    hchi = tuple(f.add(*(f.mul(hessian[i][j], chi[j]) for j in range(3))) for i in range(3))
    kgp = f.add(*(f.scale(f.derivative(pressure, j), direction[j]) for j in range(3)))
    forcing = tuple(f.add(f.scale(f.mul(adot, {g: s.diff(c, time) for g, c in mean[j].items()}), -2),
                         f.mul(kgp, mean[j]), f.scale(ugp, direction[j])) for j in range(3))
    covariance = f.mul(adot, adot).get(f.ZERO, 0)
    reduced = (f.inner(dchi, dchi)-f.inner(chi, hchi)+covariance*q*q)/2+f.inner(chi, forcing)
    boundary = s.diff(f.inner(tuple(f.mul(adot, x) for x in mean), chi), time)
    print("Exposing action residual:", s.simplify(direct-reduced-boundary))
    ledger.check("full Euler material action yields the cell response after its exact time boundary",
                 s.simplify(direct-reduced-boundary) == 0)
    wrong_forcing = tuple(f.add(forcing[j], f.scale(ugp, -direction[j])) for j in range(3))
    ledger.check("omitting the divergence-pressure term changes the actual action",
                 s.simplify(f.inner(chi, forcing)-f.inner(chi, wrong_forcing)) != 0)

    # Frozen-in compensation: an actual divergence-free axial shear map.
    x, y, z, k, gamma = s.symbols("x y z k gamma", real=True, nonzero=True)
    displacement = s.Matrix([s.cos(k*z), s.sin(k*z), 0])
    compensation = s.Matrix([0, 0, -gamma])
    axes = (x, y, z)
    pushed = s.Matrix([sum(compensation[j]*s.diff(displacement[i], axes[j])
                             - displacement[j]*s.diff(compensation[i], axes[j]) for j in range(3))
                       for i in range(3)])
    ledger.check("uniform compensating vorticity actually bends under a three-dimensional Kelvin map",
                 sum(s.diff(displacement[j], axes[j]) for j in range(3)) == 0
                 and pushed == -gamma*displacement.diff(z) and pushed != s.zeros(3, 1))
    ledger.check("its missing response vanishes in the genuinely planar limit only",
                 pushed.subs(k, 0) == s.zeros(3, 1))
    ell, carrier = s.symbols("ell carrier", integer=True, nonzero=True)
    axial_cross = s.integrate(s.exp(s.I*ell*z), (z, 0, 2*s.pi))
    ledger.check("distinct integer axial sectors have zero full-fluid quadratic cross action",
                 axial_cross == 0)
    print("Derived covariance for the exposing Euler field:", covariance)
    print("Nonzero second-order pressure border:", s.simplify(pressure_border))
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
