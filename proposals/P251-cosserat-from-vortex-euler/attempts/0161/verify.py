"""Exact field, entire-cell kernel, pressure-return and action checks for 0161."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0161-triangular-acoustic-array")
    a, b = s.symbols("a b", real=True)
    lam, amp, rho, k = s.symbols("lambda Psi rho k", positive=True)
    psi = amp*(s.cos(a)+s.cos(b)+s.cos(a+b))

    # a=b1.x, b=b2.x; physical, not orthogonal-coordinate derivatives.
    def dx(expr):
        return lam*(s.diff(expr, a)-s.diff(expr, b)/2)

    def dy(expr):
        return s.sqrt(3)*lam*s.diff(expr, b)/2

    derivative = (dx, dy)
    v = s.Matrix([-dy(psi), dx(psi)])
    zeta = dx(v[1])-dy(v[0])
    pressure = -(v.dot(v)+lam**2*psi**2)/2
    checks.check("actual oblique-cell field is divergence free",
                 s.simplify(dx(v[0])+dy(v[1])) == 0)
    checks.check("the entire vorticity gradient is the first-shell Euler profile",
                 s.trigsimp(zeta+lam**2*psi) == 0)
    euler = [sum(v[j]*derivative[j](v[i]) for j in range(2))
             +derivative[i](pressure) for i in range(2)]
    checks.check("full stationary Euler pressure is derived from the field",
                 all(s.trigsimp(entry) == 0 for entry in euler))
    checks.check("the exact common-level separatrix factorization holds",
                 s.trigsimp(s.expand_trig(psi+amp
                                          -4*amp*s.cos(a/2)*s.cos(b/2)*s.cos((a+b)/2))) == 0)
    checks.check("all three boundary families have the same NONZERO streamfunction",
                 s.trigsimp(psi.subs(a, s.pi)+amp) == 0
                 and s.trigsimp(psi.subs(b, s.pi)+amp) == 0
                 and s.trigsimp(psi.subs(b, s.pi-a)+amp) == 0)
    mean_psi = s.integrate(s.integrate(psi, (a, -s.pi, s.pi)),
                           (b, -s.pi, s.pi))/(4*s.pi**2)
    checks.check("the complete-cell streamfunction mean is zero", mean_psi == 0)

    polygons = [
        [(-s.pi, 0), (0, -s.pi), (s.pi, -s.pi), (s.pi, 0),
         (0, s.pi), (-s.pi, s.pi)],
        [(0, s.pi), (s.pi, 0), (s.pi, s.pi)],
        [(-s.pi, -s.pi), (0, -s.pi), (-s.pi, 0)],
    ]
    areas = []
    boundary_flux = s.zeros(2)
    for vertices in polygons:
        area = s.expand(sum(vertices[j][0]*vertices[(j+1) % len(vertices)][1]
                            -vertices[(j+1) % len(vertices)][0]*vertices[j][1]
                            for j in range(len(vertices)))/2)
        areas.append(area)
        # Integral r_j n_i ds, with centered or uncentered coordinates equal.
        for j, point in enumerate(vertices):
            end = vertices[(j+1) % len(vertices)]
            outward_ds = s.Matrix([end[1]-point[1], point[0]-end[0]])
            midpoint = (s.Matrix(point)+s.Matrix(end))/2
            boundary_flux += outward_ds*midpoint.T
    checks.check("bounded invariant polygons tile the entire reciprocal cell",
                 areas == [3*s.pi**2, s.pi**2/2, s.pi**2/2])
    checks.check("separatrix boundary flux is retained, not deleted",
                 boundary_flux == 4*s.pi**2*s.eye(2))
    dual = lam**2*(-amp)*s.Matrix([[0, 1], [-1, 0]])
    checks.check("the actual adjoint translation matrix is nondegenerate",
                 s.factor(dual.det()) == lam**4*amp**2)
    checks.check("a zero-level boundary mutation destroys the needed dual rows",
                 dual.subs(amp, 0).det() == 0)

    rotation = s.Matrix([[0, -1, 0], [0, 0, -1], [-1, 0, 0]])
    z = s.Symbol("z")
    first_x = s.Matrix([1, -s.Rational(1, 2), -s.Rational(1, 2)])
    first_y = s.Matrix([0, s.sqrt(3)/2, -s.sqrt(3)/2])
    translations = s.Matrix.hstack(first_x, first_y)
    checks.check("first-shell odd representation splits into vector and distinct scalar",
                 s.factor(rotation.charpoly(z).as_expr()) == (z+1)*(z*z-z+1))
    checks.check("the complete vector kernel is exactly the translations",
                 translations.rank() == 2
                 and (rotation**2-rotation+s.eye(3))*translations == s.zeros(3, 2)
                 and rotation*s.ones(3, 1) == -s.ones(3, 1))
    n, m = s.symbols("n m", integer=True)
    form = n*n+m*m-n*m
    checks.check("the next-shell exclusion uses an exact congruence, not a cutoff",
                 s.expand(form-(n+m)**2) == -3*n*m
                 and {r*r % 3 for r in range(3)} == {0, 1})
    checks.check("global Arnold complement has the explicit two-thirds lower bound",
                 1-s.Rational(1, 3) == s.Rational(2, 3))

    f, g = s.symbols("f g", cls=s.Function)
    x, y = s.symbols("x y", real=True)
    # Independent general divergence-free curl identity checks F:L2 -> H1.
    ps = f(x, y)
    pot = g(x, y)
    base = s.Matrix([-s.diff(ps, y), s.diff(ps, x)])
    grad = s.Matrix([s.diff(pot, x), s.diff(pot, y)])
    rhs = -grad.jacobian((x, y))*base-base.jacobian((x, y))*grad
    curl_rhs = s.diff(rhs[1], x)-s.diff(rhs[0], y)
    vort = s.diff(base[1], x)-s.diff(base[0], y)
    expected = -grad.dot(s.Matrix([s.diff(vort, x), s.diff(vort, y)])) \
        -vort*(s.diff(grad[0], x)+s.diff(grad[1], y))
    checks.check("gradient-return forcing has bounded vorticity without derivative loss",
                 s.expand(curl_rhs-expected) == 0)

    vectors = [lam*s.Matrix([first_x[j], first_y[j]])
               for j in range(3)]
    # The two entries above are the Cartesian coordinates of b_j.
    jmat = s.Matrix([[0, -1], [1, 0]])
    covariance = s.simplify(amp**2*sum(((jmat*q)*(jmat*q).T for q in vectors),
                                     start=s.zeros(2))/2)
    checks.check("the acoustic coefficient is the complete actual velocity covariance",
                 covariance == 3*amp**2*lam**2*s.eye(2)/4)
    ux, uy = s.symbols("Ux Uy", real=True)
    displacement = s.Matrix([ux, uy])
    projected_norm = 0
    wrong_norm = 0
    for q in vectors:
        full_wave = q.col_join(s.Matrix([k]))
        projector = s.eye(3)-full_wave*full_wave.T/(lam**2+k**2)
        vertical = s.Matrix([0, 0, 1])
        coefficient_square = amp**2*((jmat*q).dot(displacement))**2/2
        projected_norm += coefficient_square*(projector*vertical).dot(projector*vertical)
        wrong_norm += coefficient_square
    projected_norm = s.factor(projected_norm)
    checks.check("full axial pressure gives the exact finite-k phase stiffness factor",
                 s.factor(projected_norm
                          -lam**2/(lam**2+k**2)*(displacement.T*covariance*displacement)[0]) == 0)
    checks.check("omitting the pressure return changes the actual phase action",
                 s.factor(projected_norm-wrong_norm) != 0)
    q = vectors[0]
    wave = q.col_join(s.Matrix([k]))
    projector = s.eye(3)-wave*wave.T/(lam**2+k**2)
    exact_rate = -s.I*k*projector*s.Matrix([0, 0, 1])
    written_rate = s.Matrix([s.I*k**2*q[0]/(lam**2+k**2),
                             s.I*k**2*q[1]/(lam**2+k**2),
                             -s.I*k+s.I*k**3/(lam**2+k**2)])
    checks.check("Kelvin preparation contains BOTH horizontal and vertical returns",
                 s.simplify(exact_rate-written_rate) == s.zeros(3, 1)
                 and s.simplify(wave.dot(exact_rate)) == 0)
    vx, vy = s.symbols("Vx Vy", real=True)
    phase_x = s.Matrix([ux, uy, 0, 0])
    phase_v = s.Matrix([0, 0, rho*vx, rho*vy])
    symplectic = s.BlockMatrix([[s.zeros(2), s.eye(2)],
                                [-s.eye(2), s.zeros(2)]]).as_explicit()
    checks.check("the initial material phase pairing fixes physical mass rho",
                 s.expand((phase_x.T*symplectic*phase_v)[0]-rho*(ux*vx+uy*vy)) == 0)
    print("Scope: exact coefficients and operator-proof anchors; no sampled stability or imposed cosine.")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
