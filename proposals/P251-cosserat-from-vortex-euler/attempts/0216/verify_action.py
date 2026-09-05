"""Full three-component Euler/Lin action of actual homogeneous tag data."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0216-full-label-action")
    a, b = s.symbols("a b", real=True)
    aa = s.symbols("A", positive=True)
    psi = s.cos(b)+aa*s.cos(a)
    u = s.Matrix([psi, s.sin(b), -aa*s.sin(a)])
    grad = u.jacobian([s.Symbol("X"), a, b])
    pressure = -u.dot(u)/2
    hp = s.hessian(pressure, [s.Symbol("X"), a, b])

    def transport(v):
        return s.sin(b)*s.diff(v, a)-aa*s.sin(a)*s.diff(v, b)

    def integrate(v):
        return s.simplify(s.integrate(s.expand_trig(s.expand(v)), (a, -s.pi, s.pi), (b, -s.pi, s.pi)))

    streams = [s.sin(a), s.cos(a)*s.sin(b), s.sin(2*a)*s.cos(b)]
    fields = []
    for index, stream in enumerate(streams):
        xi = s.Matrix([0, -s.diff(stream, b), s.diff(stream, a)])
        rate = grad*xi-transport(xi)
        ts = transport(stream)
        expected = s.Matrix([-ts, s.diff(ts, b), -s.diff(ts, a)])
        ledger.check(f"full homogeneous Lin field {index}", s.simplify(rate-expected) == s.zeros(3, 1))
        full_energy = integrate((rate.dot(rate)-transport(xi).dot(transport(xi))+xi.dot(hp*xi))/2)
        target = integrate(ts**2)
        ledger.check(f"full pressure and axial energy identity {index}", s.simplify(full_energy-target) == 0)
        no_axial = full_energy-integrate(rate[0]**2)/2
        ledger.check(f"dropping axial coordinate velocity changes energy {index}", s.simplify(no_axial-full_energy) != 0)
        fields.append(xi)
    phase = integrate(fields[0].dot(grad*fields[1])-fields[1].dot(grad*fields[0]))
    target_phase = -integrate(streams[0]*transport(streams[1]))
    ledger.check("full cotangent phase equals actual transport pairing", s.simplify(phase-target_phase) == 0)
    ledger.check("actual full cotangent phase is nonzero on an exposing pair", phase != 0)
    ledger.check("physical Euler velocity is zero but cotangent is not", grad*fields[0] != s.zeros(3, 1))
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
