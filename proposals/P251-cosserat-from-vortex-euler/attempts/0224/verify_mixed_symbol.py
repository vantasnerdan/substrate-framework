"""Independent frozen-coefficient mixed-c Bloch symbol audit; not PDE signs."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0224-mixed-c-symbol")
    kx, ny, nz = s.symbols("Kx qy qz", real=True)
    ci, cj = s.symbols("ci cj", real=True)
    wx, wy, wz = s.symbols("wx wy wz", real=True)
    wave = s.Matrix([kx, ny, nz])
    omega = s.Matrix([wx, wy, wz])
    normal_square = ny**2+nz**2
    full_square = wave.dot(wave)
    projection = s.eye(3)-wave*wave.T/full_square

    def generator(c):
        return s.Matrix([c, -s.I*nz-kx*c*ny/normal_square,
                         s.I*ny-kx*c*nz/normal_square])

    xi, xj = generator(ci), generator(cj)
    checks.check("actual completed Bloch generator is solenoidal for arbitrary c",
                 s.simplify(wave.dot(xi)) == 0 and s.simplify(wave.dot(xj)) == 0)
    phase = s.simplify(omega.dot(xi.conjugate().cross(xj)))
    checks.check("mixed-c frozen phase has no order-two principal term",
                 s.simplify(phase+s.I*(ci+cj)*omega.dot(wave)) == 0)
    vi = (projection*xi.cross(omega)).applyfunc(s.simplify)
    vj = (projection*xj.cross(omega)).applyfunc(s.simplify)
    energy = s.simplify(vi.conjugate().dot(vj+s.I*wave.cross(vj)))
    principal = omega.dot(wave)**2*(normal_square/full_square-ci-cj)
    lower = ci*cj*omega.dot(wave)**2/normal_square
    checks.check("full Leray helicity cancels cubic order also between distinct c sectors",
                 s.simplify(energy-principal-lower) == 0)
    variables = (kx, ny, nz)

    def dilation(value):
        return s.simplify(sum(q*s.diff(value, q) for q in variables))

    checks.check("surviving frozen mixed energy symbol has degree two and degree zero",
                 s.simplify(dilation(principal)-2*principal) == 0
                 and dilation(lower) == 0)
    checks.check("every first covariant-frequency derivative lowers its degree",
                 all(s.simplify(dilation(s.diff(principal, q))-s.diff(principal, q)) == 0
                     for q in variables))
    checks.check("every second covariant-frequency energy derivative has degree at most zero",
                 all(dilation(s.diff(principal, q, p)) == 0
                     for q in variables for p in variables))
    checks.check("axial and transverse covariant-frequency phase second derivatives vanish here",
                 all(s.diff(phase, q, p) == 0 for q in variables for p in variables))
    wrong = s.Matrix([ci, -s.I*nz, s.I*ny])
    checks.check("omitting the actual completion exposes a nonzero Bloch divergence",
                 s.simplify(wave.dot(wrong)-kx*ci) == 0
                 and wave.dot(wrong) != 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
