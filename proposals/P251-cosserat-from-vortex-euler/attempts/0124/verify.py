"""Exact affine-Euler, physical-moment and scale-exponent oracle."""

import sympy as s


def main():
    count = 0

    def check(label, expression):
        nonlocal count
        entries = list(expression) if isinstance(expression, s.MatrixBase) else [expression]
        assert all(s.simplify(v) == 0 for v in entries), (label, expression)
        count += 1
        print(f"PASS {label}")

    av = s.symbols("a0:8", real=True)
    a = s.Matrix([[av[0], av[1], av[2]], [av[3], av[4], av[5]],
                  [av[6], av[7], -av[0] - av[4]]])
    h0, h1, h2, h3, h4 = s.symbols("h0:5", real=True)
    h = s.Matrix([[h0, h1, h2], [h1, h3, h4],
                  [h2, h4, -s.trace(a * a) - h0 - h3]])
    adot = -a * a - h
    check("affine incompressibility", s.trace(a))
    check("pressure Poisson preserves incompressibility", s.trace(adot))
    check("affine acceleration is exact pressure gradient", adot + a * a + h)
    check("affine acceleration curl vanishes", adot + a * a - (adot + a * a).T)
    omega = s.Matrix([a[2, 1] - a[1, 2], a[0, 2] - a[2, 0], a[1, 0] - a[0, 1]])
    omega_dot = s.Matrix([adot[2, 1] - adot[1, 2], adot[0, 2] - adot[2, 0], adot[1, 0] - adot[0, 1]])
    check("actual affine vorticity transport", omega_dot - a * omega)

    ss, tt, aa, bb = s.symbols("shear tilt scale shear2", nonzero=True, real=True)
    m = s.Matrix([[1, ss, tt], [0, aa, bb], [0, 0, 1 / aa]])
    q = s.Matrix(s.symbols("q0:3", real=True))
    w = s.Matrix(s.symbols("w0:3", real=True))
    k = s.Matrix(s.symbols("k0:3", real=True))
    check("material volume preservation fixture", m.det() - 1)
    check("nonorthogonal cross product covariance", (m * q).cross(m * w) - m.inv().T * q.cross(w))
    g = m.inv() * m.inv().T
    kt = m.inv().T * k
    pressure = s.eye(3) - kt * kt.T / (kt.dot(kt))
    material_pressure = g - (g * k) * (g * k).T / (k.dot(g * k))
    check("complete affine Leray material symbol", m.inv() * pressure * m.inv().T - material_pressure)
    check("material symbol preserves divergence", k.T * material_pressure)

    om, zz, xr, yr, sigma, sz = s.symbols("Omega z x y sigma_r sigma_z", real=True)
    qx, qy = s.symbols("qx qy", real=True)
    pos = s.Matrix([xr, yr, zz])
    axis = s.Matrix([0, 0, 1])
    q3 = s.Matrix([qx, qy, 0])
    bg = om * axis.cross(pos)
    xi = q3 * zz
    v = -2 * om * axis.cross(q3) * zz
    material_spin = xi.cross(bg) + pos.cross(v + bg.jacobian(pos) * xi)
    covariance = {xr**2: sigma, yr**2: sigma, zz**2: sz, xr*yr: 0, xr*zz: 0, yr*zz: 0}
    averaged_spin = material_spin.applyfunc(lambda e: s.expand(e).subs(covariance, simultaneous=True))
    check("full material observation recovers actual shear spin", averaged_spin - om * sz * q3)
    check("physical core angle rate", axis.cross(-om * axis.cross(q3)) - om * q3)

    theta = s.Rational(11, 32)
    check("interpolation reaches H eleven-fourths", 8 * theta - s.Rational(11, 4))
    power_a = 3 * (1 - theta) + theta
    power_ell = s.Rational(1, 2) * (1 - theta) - s.Rational(13, 2) * theta
    check("interpolated radial exponent", power_a - s.Rational(37, 16))
    check("interpolated axial exponent", power_ell + s.Rational(61, 32))
    p, qscale = s.Rational(1, 8), s.Rational(1, 16)
    check("strictly decaying C1 transfer exponent", power_a + (1 + p) * power_ell - s.Rational(43, 256))
    check("strictly decaying physical spin exponent", 3 + (1 + qscale) / 2 - 5 * (1 + p) / 2 - s.Rational(23, 32))
    assert p > qscale > 0
    assert s.Rational(43, 256) > 0 and s.Rational(23, 32) > 0
    print("PASS hierarchy and both physical errors have positive powers")
    count += 1
    check("initial Kelvin mismatch is smaller than dynamic remainder", (2 + s.Rational(3, 2) * (1 + p)) - (3 + (1 + p) / 2) - p)
    check("exterior tail is smaller than dynamic remainder", (2 + 3 * (1 + p)) - (3 + (1 + p) / 2) - (s.Rational(3, 2) + 5 * p / 2))

    eps, core, tail, amp, cu, du, a1, b1, c0 = s.symbols("eps ell R amp C D a1 b1 c0", positive=True)
    root_equation = 2 * cu * eps**2 * tail**3 * a1 * amp**2 - du * eps**2 * tail**3 * b1 * amp + eps**2 * core**3 * c0
    scaled = s.cancel(root_equation.subs(amp, core**3 * s.Symbol("s") / tail**3) / (eps**2 * core**3))
    check("simple-root equation has nonvanishing normalized linear term", scaled - (2 * cu * a1 * (core / tail)**3 * s.Symbol("s")**2 - du * b1 * s.Symbol("s") + c0))
    print(f"{count}/{count} exact checks passed")


if __name__ == "__main__":
    main()
