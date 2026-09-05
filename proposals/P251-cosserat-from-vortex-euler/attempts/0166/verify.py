"""Exact triangular-core Euler, optical pressure, KKS and material-moment checks."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0166")
    x, y = sp.symbols("x y", real=True)
    lam, amp, rho, omega, axial_period = sp.symbols("lambda Psi rho Omega L_s", positive=True)
    waves = [lam * sp.Matrix([1, 0]), lam * sp.Matrix([-sp.Rational(1, 2), sp.sqrt(3) / 2]), lam * sp.Matrix([-sp.Rational(1, 2), -sp.sqrt(3) / 2])]
    position = sp.Matrix([x, y])
    phases = [wave.dot(position) for wave in waves]
    psi = amp * sum(sp.cos(phase) for phase in phases)
    base = sp.Matrix([-sp.diff(psi, y), sp.diff(psi, x), -lam * psi])
    vorticity = sp.Matrix([sp.diff(base[2], y), -sp.diff(base[2], x), sp.diff(base[1], x) - sp.diff(base[0], y)])
    ledger.check("actual triangular field is incompressible", sp.simplify(sp.diff(base[0], x) + sp.diff(base[1], y)) == 0)
    ledger.check("all components have the same constant curl eigenvalue", all(sp.trigsimp(entry) == 0 for entry in vorticity - lam * base))
    pressure = -base.dot(base) / 2
    acceleration = base.jacobian((x, y)) * base[:2, :]
    stationary = acceleration + sp.Matrix([sp.diff(pressure, x), sp.diff(pressure, y), 0])
    ledger.check("full stationary Euler pressure including axial field", all(sp.trigsimp(entry) == 0 for entry in stationary))
    core_hessian = sp.hessian(psi, (x, y)).subs({x: 0, y: 0})
    ledger.check("actual circular core Hessian", core_hessian == -3 * amp * lam**2 * sp.eye(2) / 2)
    ledger.check("axial core speed is derived, not prescribed", base[2].subs({x: 0, y: 0}) == -3 * lam * amp)
    sixth_harmonic = x**6 - 15 * x**4 * y**2 + 15 * x**2 * y**4 - y**6
    radius2 = x**2 + y**2
    jet = sp.expand(amp * sum(sum((-1) ** order * phase ** (2 * order) / sp.factorial(2 * order) for order in range(4)) for phase in phases))
    predicted = 3 * amp - 3 * amp * lam**2 * radius2 / 4 + 3 * amp * lam**4 * radius2**2 / 64 - amp * lam**6 * (10 * radius2**3 + sixth_harmonic) / 7680
    ledger.check("full sixth-order Cartesian triangular jet", sp.expand(jet - predicted) == 0)
    ledger.check("first nonradial term really is harmonic six", sp.expand(sp.diff(sixth_harmonic, x, 2) + sp.diff(sixth_harmonic, y, 2)) == 0)
    radial_series = 3 * amp * sum((-1) ** j * (lam**2 * radius2 / 4) ** j / sp.factorial(j) ** 2 for j in range(4))
    ledger.check("actual radial jet matches negative-amplitude Lundquist", sp.expand(jet - radial_series + amp * lam**6 * sixth_harmonic / 7680) == 0)
    ledger.check("a quartic anisotropy mutation breaks the field jet", sp.expand(jet - predicted + amp * lam**4 * (x**4 - 6 * x**2 * y**2 + y**4)) != 0)

    # Derive the pressure identity from divergence plus the full axial equation.
    carrier = sp.symbols("k", nonzero=True, real=True)
    vx, vy, pfield = sp.Function("Vx")(x, y), sp.Function("Vy")(x, y), sp.Function("P")(x, y)
    vector = sp.Matrix([vx, vy])
    rotation = sp.Matrix([[0, -1], [1, 0]])
    background = -omega * rotation * position
    axial = lam * omega * radius2 / 2
    divergence = sp.diff(vx, x) + sp.diff(vy, y)
    axial_velocity = sp.I * divergence / carrier
    normal_rhs = -vector.jacobian((x, y)) * background - background.jacobian((x, y)) * vector - sp.I * carrier * axial * vector - sp.Matrix([sp.diff(pfield, x), sp.diff(pfield, y)])
    from_normal = sp.I * (sp.diff(normal_rhs[0], x) + sp.diff(normal_rhs[1], y)) / carrier
    from_axial = -background[0] * sp.diff(axial_velocity, x) - background[1] * sp.diff(axial_velocity, y) - sp.I * carrier * axial * axial_velocity - vector.dot(sp.Matrix([sp.diff(axial, x), sp.diff(axial, y)])) - sp.I * carrier * pfield
    pressure_equation = carrier**2 * pfield - sp.diff(pfield, x, 2) - sp.diff(pfield, y, 2) - 2 * omega * (sp.diff(vy, x) - sp.diff(vx, y)) - 2 * sp.I * carrier * vector.dot(sp.Matrix([sp.diff(axial, x), sp.diff(axial, y)]))
    ledger.check("pressure equation follows from full Euler divergence", sp.simplify((from_normal - from_axial) * carrier / sp.I - pressure_equation) == 0)
    p, delta = sp.symbols("p delta", positive=True)
    ledger.check("trapped carrier has the required quadratic sign", (carrier * axial).subs(carrier, -p) == -p * lam * omega * radius2 / 2)
    ledger.check("opposite carrier gives inverted quadratic potential", sp.simplify((carrier * axial).subs(carrier, p) + (carrier * axial).subs(carrier, -p)) == 0 and (carrier * axial).subs(carrier, p) != (carrier * axial).subs(carrier, -p))
    ell4 = 2 / (lam * p**3)
    ledger.check("pressure and axial well fix the mode width", sp.simplify(omega / (p**2 * sp.sqrt(ell4)) - p * lam * omega * sp.sqrt(ell4) / 2) == 0)
    dimensionless_fourth = sp.simplify((lam**4 * ell4).subs(p, 2 * lam / delta**2))
    ledger.check("first triangular transverse error has order delta-six", dimensionless_fourth == delta**6 / 4)
    axial_sixth = sp.simplify(((p / lam) * (lam**4 * ell4) ** sp.Rational(3, 2)).subs(p, 2 * lam / delta**2))
    ledger.check("first triangular axial error has order delta-seven", axial_sixth == delta**7 / 4)

    real, imaginary = sp.symbols("a b", real=True)
    plus = sp.Matrix([1, sp.I]) / sp.sqrt(2)
    velocity = plus * (real + sp.I * imaginary)
    displacement = sp.I * velocity / (2 * omega)
    velocity3 = sp.Matrix([velocity[0], velocity[1], 0])
    displacement3 = sp.Matrix([displacement[0], displacement[1], 0])
    background3 = sp.Matrix([background[0], background[1], 0])
    lagrangian_velocity = velocity3 + background3.jacobian((x, y)) * displacement
    literal_spin = (displacement3.cross(background3) + sp.Matrix([x, y, 0]).cross(lagrangian_velocity))[2]
    ledger.check("principal tagged spin cancels before pressure pointwise", sp.simplify(literal_spin) == 0)
    kks_density = rho * sp.Matrix([0, 0, -2 * omega]).dot(sp.re(displacement3).cross(sp.im(displacement3)))
    ledger.check("actual core KKS sign and factor", sp.simplify(kks_density + rho * (real**2 + imaginary**2) / (4 * omega)) == 0)
    ledger.check("displacement angular moment is not silently spin", sp.simplify(sp.Matrix([x, y, 0]).cross(displacement3)[2]) != 0)

    radial, w = sp.symbols("r2 w", positive=True)
    laguerre = sp.assoc_laguerre(2, 1, radial)
    norm_poly = sp.Poly(radial * laguerre**2, radial)
    radial_norm = sum(coef * sp.factorial(degree[0]) for degree, coef in norm_poly.terms())
    ledger.check("full radial action norm for the newly selected mode", radial_norm == 3)
    laplace = sum(coef * sp.factorial(degree[0] + 2) * (2 / w) ** (degree[0] + 3) for degree, coef in sp.Poly(laguerre, radial).terms())
    slope = sp.simplify(sp.diff(laplace, w).subs(w, 1) / laplace.subs(w, 1))
    ledger.check("actual quadrupole material phase coefficient", slope == -sp.Rational(19, 3))
    gamma = 2 * omega + omega * delta / 3
    square = sp.series(gamma**2, delta, 0, 2).removeO()
    curvature = sp.expand(delta**2 * sp.diff(square, delta, 2) / 4 + 3 * delta * sp.diff(square, delta) / 4)
    ledger.check("physical clock has positive natural-scale curvature", curvature == omega**2 * delta)
    pressure_poly = sp.expand(laguerre - 2 * sp.diff(laguerre, radial))
    ledger.check("pressure spin uses its new radial polynomial", pressure_poly == radial**2 / 2 - 5 * radial + 9)

    def carrier_row(poly: sp.Expr) -> sp.Expr:
        return sp.expand(-poly / 2 + sp.Rational(3, 2) * radial * (sp.diff(poly, radial) - poly / 2))

    rows = [sp.expand(poly) for poly in [pressure_poly, carrier_row(pressure_poly), carrier_row(carrier_row(pressure_poly)), radial * pressure_poly, carrier_row(radial * pressure_poly), carrier_row(carrier_row(radial * pressure_poly))]]
    coefficients = sp.Matrix([[poly.coeff(radial, j) for j in range(6)] for poly in rows])
    ledger.check("six time/carrier rows independent for n2 m2", coefficients.det() != 0)
    jets = [[sp.factorial(h) if j == h else 0 for j in range(8)] for h in range(2)]
    for poly in rows:
        jets.append([sum(sp.binomial(j, h) * (-sp.Rational(1, 2)) ** (j - h) * sp.factorial(h) * poly.coeff(radial, h) for h in range(j + 1)) for j in range(8)])
    moment_minor = sp.Matrix(jets).det(method="domain-ge")
    ledger.check("actual reference/background extension has nonzero scalar minor", not moment_minor.has(radial) and moment_minor != 0)
    for radial_index, angular_index in ((2, 2), (8, 5)):
        energy = 2 * radial_index + angular_index
        cutoff_poly = sp.expand((radial / 2 - energy - sp.Rational(angular_index, energy)) * sp.assoc_laguerre(radial_index, angular_index - 1, radial))
        cutoff_rows = [cutoff_poly]
        for _ in range(2):
            previous = cutoff_rows[-1]
            cutoff_rows.append(sp.expand(sp.Rational(3, 2) * radial * (sp.diff(previous, radial) - previous / 2)))
        cutoff_matrix = sp.Matrix([[poly.coeff(radial, j) for j in range(radial_index + 4)] for poly in cutoff_rows])
        ledger.check(f"three cutoff carrier functionals independent n{radial_index} m{angular_index}", cutoff_matrix.rank() == 3)
        homogeneous = [sum(coef * sp.factorial(degree[0] + angular_index) * 2 ** (degree[0] + angular_index + 1) for degree, coef in sp.Poly(poly, radial).terms()) for poly in cutoff_rows]
        ledger.check(f"all three ideal homogeneous phase jets vanish n{radial_index} m{angular_index}", all(value == 0 for value in homogeneous))
        two_point = sp.Matrix([[poly.subs(radial, point) for point in (1, 2)] for poly in cutoff_rows[:2]])
        ledger.check(f"one cutoff value does not force its derivative n{radial_index} m{angular_index}", two_point.det() != 0)
    print(f"Omega=3*Psi*lambda^2/2; ell^4={ell4}; delta=sqrt(2*lambda/p)")
    print(f"radial norm={radial_norm}; physical phase slope={slope}; curvature={curvature}")
    print(f"eight-row n2,m2 moment minor={moment_minor}; axial action factor=L_s")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
