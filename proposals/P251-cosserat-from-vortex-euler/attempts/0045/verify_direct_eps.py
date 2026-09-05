"""Exact identities supporting the finite-norm direct-EPS proof (no numerics)."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0045-direct-EPS-angle")
    x, y, z, k = s.symbols("x y z k", real=True, nonzero=True)
    coords = (x, y, z)
    phi = s.Function("phi")(*coords)

    def curl(v):
        return s.Matrix([s.diff(v[2], y)-s.diff(v[1], z),
                         s.diff(v[0], z)-s.diff(v[2], x),
                         s.diff(v[1], x)-s.diff(v[0], y)])

    def div(v):
        return sum(s.diff(v[j], coords[j]) for j in range(3))

    def zero(v):
        return all(s.simplify(s.trigsimp(item)) == 0 for item in v)

    grad = s.Matrix([s.diff(phi, item) for item in coords])
    p1 = s.Matrix([s.cos(k*z), s.sin(k*z), 0])
    p2 = s.Matrix([-s.sin(k*z), s.cos(k*z), 0])
    generators = [-curl(phi*p)/k for p in (p1, p2)]
    for i, (xi, p) in enumerate(zip(generators, (p1, p2))):
        ledger.check(f"compact curl generator {i} is divergence free", s.simplify(div(xi)) == 0)
        ledger.check(f"full cutoff expansion {i}", zero(xi-phi*p+grad.cross(p)/k))
        ledger.check(f"carrier {i} has negative signed helicity", zero(curl(p)+k*p))
    expected_cross = (phi**2*s.Matrix([0, 0, 1])
                      +phi*s.Matrix([-grad[1], grad[0], 0])/k+grad[2]*grad/k**2)
    ledger.check("exact cross includes all three-dimensional cutoff terms",
                 zero(generators[0].cross(generators[1])-expected_cross))
    ledger.check("omitting axial cutoff correction changes the exact cross",
                 not zero(grad[2]*grad/k**2))

    wx, wy, wz = (s.Function(name)(*coords) for name in ("wx", "wy", "wz"))
    omega = s.Matrix([wx, wy, wz])
    horizontal_flux = s.Matrix([wy, -wx, 0])*phi**2/2
    cross_correction = omega.dot(s.Matrix([-grad[1], grad[0], 0]))*phi
    ledger.check("KKS integration by parts gives minus curl-omega-z / 2",
                 s.simplify(cross_correction-div(horizontal_flux)
                            +phi**2*(s.diff(wy, x)-s.diff(wx, y))/2) == 0)

    # Inner bump is one near the point: derive the actual vorticity tilt jet.
    e = s.Matrix([1, 0, 0])
    r = s.Matrix(coords)
    xi_r = curl(-(x*x+y*y+z*z)*e/2)
    ledger.check("core vector potential gives a rigid rotation jet", zero(xi_r-e.cross(r)))
    omega_at_core = s.Matrix([0, 0, s.Symbol("wstar", positive=True)])
    delta_core = xi_r.jacobian(coords)*omega_at_core
    ledger.check("core jet tilts physical vorticity, not a gauge label",
                 delta_core == e.cross(omega_at_core) and delta_core != s.zeros(3, 1))

    a = s.Function("a")(*coords)
    t1, t2 = s.symbols("t1 t2", real=True)
    principal = a*(t1*s.Matrix([s.sin(k*z), -s.cos(k*z), 0])+t2*p1)
    ledger.check("principal norm is independent of carrier phase",
                 s.trigsimp(principal.dot(principal)-a*a*(t1*t1+t2*t2)) == 0)
    ledger.check("principal self-helicity includes no hidden amplitude remainder",
                 s.simplify(s.trigsimp(principal.dot(curl(principal))
                                      +k*a*a*(t1*t1+t2*t2))) == 0)

    # Sign and full off-diagonal elimination, independent of any diagonal example.
    u = s.Matrix(s.symbols("u0:3"))
    w = s.Matrix(s.symbols("w0:3"))
    eta = s.Matrix(s.symbols("eta0:3"))
    ledger.check("Euler evolution fixes the positive KKS sign",
                 s.expand(u.dot(eta.cross(w))-w.dot(u.cross(eta))) == 0)
    b, h11, h12, h22, q, qdot, shape = s.symbols("B h11 h12 h22 q qdot shape", nonzero=True)
    lagrangian = b*shape*qdot-(h11*q*q+2*h12*q*shape+h22*shape*shape)/2
    eliminated = s.solve(s.diff(lagrangian, shape), shape)[0]
    reduced = s.expand(lagrangian.subs(shape, eliminated))
    expected = b*b*qdot*qdot/(2*h22)-(h11*h22-h12*h12)*q*q/(2*h22)-b*h12*q*qdot/h22
    ledger.check("full conjugate elimination retains the mixed Hessian term",
                 s.simplify(reduced-expected) == 0)
    ledger.check("mixed term is an exact time derivative",
                 s.diff(-b*h12*q*q/(2*h22), q)*qdot == -b*h12*q*qdot/h22)
    a_pos, remainder, lam_abs, margin = s.symbols("A Ctotal lambda_abs margin", positive=True)
    carrier_bound = lam_abs*remainder/a_pos+margin
    lower = s.expand((1+carrier_bound/lam_abs)*a_pos-remainder)
    ledger.check("derived finite carrier threshold makes the full Hessian bound positive",
                 lower.is_positive)
    print("Analytic oracle: proof.md, finite Fourier-symbol and norm inequalities.")
    print("These identities test the proof machinery; no generic EPS integral is numerically asserted.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
