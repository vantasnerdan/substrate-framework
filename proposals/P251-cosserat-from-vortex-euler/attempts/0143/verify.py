"""Exact localization, pose-family and virial identities; no spectral numerics."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0143-global-localizable")
    x, y, z, shift = s.symbols("x y z shift", real=True)
    rho = s.Symbol("rho", positive=True)
    coords = (x, y, z)
    # A real stationary swirl exercises the full Cartesian localization
    # identity. It is only a local example, not the imported compact3D field.
    v = s.Matrix([-y, x, 0])
    p = rho*(x*x+y*y)/2
    chi = s.Function("chi")(p)
    u = chi*v
    gradp = s.Matrix([s.diff(p, c) for c in coords])
    ledger.check("actual pressure is a first integral of the local Euler field",
                 s.expand(v.dot(gradp)) == 0)
    ledger.check("pressure localization preserves complete incompressibility",
                 s.simplify(sum(s.diff(u[i], coords[i]) for i in range(3))) == 0)
    ledger.check("pressure localization preserves all Cartesian Euler components",
                 s.simplify(rho*u.jacobian(coords)*u+chi**2*gradp) == s.zeros(3, 1))
    # A non-material cutoff is a genuinely wrong construction, not a
    # tolerance change or an inserted surface pressure.
    wrong = x*v
    ledger.check("an arbitrary spatial cutoff is detected by divergence",
                 s.expand(sum(s.diff(wrong[i], coords[i]) for i in range(3))) == -y)
    ledger.check("wrong pressure without squared cutoff leaves a force residual",
                 s.simplify(rho*u.jacobian(coords)*u+gradp) != s.zeros(3, 1))

    # Differentiate a translated full nonlinear Euler residual, then compare
    # with the exact linearized residual. No stationary eigenvalue is copied.
    fields = s.Matrix([s.Function(f"u{i}")(*coords) for i in range(3)])
    pressure = s.Function("p")(*coords)
    shifted = fields.subs(x, x-shift, simultaneous=True)
    shifted_p = pressure.subs(x, x-shift, simultaneous=True)
    euler = rho*shifted.jacobian(coords)*shifted+s.Matrix(
        [s.diff(shifted_p, c) for c in coords])
    tangent = -fields.diff(x)
    linear = rho*(fields.jacobian(coords)*tangent
                  +tangent.jacobian(coords)*fields)-s.Matrix(
                      [s.diff(pressure, x, c) for c in coords])
    ledger.check("stationary translation derivative is the full Euler linearization",
                 s.simplify((euler.diff(shift).subs(shift, 0)-linear).doit())
                 == s.zeros(3, 1))
    stress = rho*fields*fields.T+pressure*s.eye(3)
    # This divergence identity, integrated with compact support, is the
    # tensor virial theorem used in the proof, including off-diagonal terms.
    virial = []
    for i in range(3):
        for k in range(3):
            actual = sum(s.diff(coords[k]*stress[i, j], coords[j]) for j in range(3))
            source = coords[k]*sum(s.diff(stress[i, j], coords[j]) for j in range(3))
            virial.append(s.expand(actual-source-stress[i, k]))
    ledger.check("full tensor stress virial identity retains pressure and all components",
                 all(value == 0 for value in virial))
    momentum = [s.expand(sum(s.diff(coords[i]*fields[j], coords[j]) for j in range(3))
                         -coords[i]*sum(s.diff(fields[j], coords[j]) for j in range(3))
                         -fields[i]) for i in range(3)]
    ledger.check("compact solenoidal mean-velocity identity is an exact divergence",
                 all(value == 0 for value in momentum))
    print("Scope: exact identities supporting the analytic global-cell construction;")
    print("no positive rigidity or general Euler/Cosserat no-go is asserted")
    raise SystemExit(ledger.finish())


if __name__ == "__main__":
    main()
