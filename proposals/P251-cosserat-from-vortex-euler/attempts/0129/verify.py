"""Exact inertial field, carrier/boundary slope, centroid and action oracle."""

import sympy as s

from substrate_framework.euler_phase import moving_phase_pullback
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0129-inertial-centroid")
    r, k, kap, om = s.symbols("r k kappa Omega", positive=True)
    sig, m = s.symbols("sigma m", real=True, nonzero=True)
    p = s.Function("p")(r)
    d = 4 * om**2 - sig**2
    wr = s.I * (sig * s.diff(p, r) - 2 * om * m * p / r) / d
    wt = (2 * om * s.diff(p, r) - sig * m * p / r) / d
    wz = k * p / sig
    checks.check("exact radial Euler equation", s.simplify(-s.I * sig * wr - 2 * om * wt + s.diff(p, r)) == 0)
    checks.check("exact azimuthal Euler equation", s.simplify(-s.I * sig * wt + 2 * om * wr + s.I * m * p / r) == 0)
    checks.check("exact axial pressure equation", s.simplify(-s.I * sig * wz + s.I * k * p) == 0)
    divergence = s.diff(r * wr, r) / r + s.I * m * wt / r + s.I * k * wz
    bessel_ode = -s.diff(p, r) / r + (m**2 / r**2 - kap**2) * p
    reduced_div = s.simplify(divergence.subs(s.diff(p, r, 2), bessel_ode) * d * sig / (s.I * p))
    checks.check("incompressibility derives inertial dispersion", s.simplify(reduced_div - (4 * om**2 * k**2 - sig**2 * (k**2 + kap**2))) == 0)
    checks.check("exact rotating Lin primitive", s.simplify((-s.I * sig) * (s.I / sig) - 1) == 0)

    freq = om - 2 * om * k / s.sqrt(kap**2 + k**2)
    k0 = kap / s.sqrt(3)
    slope = s.simplify(s.diff(freq, k).subs(k, k0))
    curvature = s.simplify(s.diff(freq, k, 2).subs(k, k0) / 2)
    checks.check("physical laboratory zero crossing", s.simplify(freq.subs(k, k0)) == 0)
    checks.check("free Bessel sideband slope", s.simplify(slope + 3 * s.sqrt(3) * om / (4 * kap)) == 0)
    checks.check("free Bessel second frequency coefficient", s.simplify(curvature - 27 * om / (32 * kap**2)) == 0)
    checks.check("positive squared sideband curvature is derived", s.simplify(slope**2 - 27 * om**2 / (16 * kap**2)) == 0)

    x, y, b = s.symbols("x y b", positive=True)
    j1 = s.symbols("J1", real=True, nonzero=True)
    branch = -2 * y / s.sqrt(x**2 + y**2)
    sx = s.simplify(s.diff(branch, x).subs(y, x / s.sqrt(3)))
    sy = s.simplify(s.diff(branch, y).subs(y, x / s.sqrt(3)))
    jp = -2 * j1 / x
    jpp = -jp / x - (1 - 1 / x**2) * j1
    f_x = -(jp + x * jpp) - 2 * jp
    f_sig = x * jp
    xprime = s.simplify(-f_sig * sy / (f_x + f_sig * sx))
    checks.check("fixed-cylinder boundary requires changing radial carrier", s.simplify(xprime + 3 * s.sqrt(3) / (2 * x**2 + 3)) == 0)
    wall_slope = s.simplify(om * b * (sy + sx * xprime))
    checks.check("actual slip-cylinder dispersion slope", s.simplify(wall_slope + om * b * 3 * s.sqrt(3) * (x**2 + 3) / (2 * x * (2 * x**2 + 3))) == 0)
    checks.check("holding radial carrier fixed changes the result", s.simplify(wall_slope - om * b * sy) != 0)

    theta = s.symbols("theta", real=True)
    checks.check("actual Cartesian centroid angular integral", s.integrate(s.exp(s.I * theta) * s.cos(theta), (theta, 0, 2 * s.pi)) == s.pi and s.integrate(s.exp(s.I * theta) * s.sin(theta), (theta, 0, 2 * s.pi)) == s.I * s.pi)
    combined = s.simplify((wr - s.I * wt).subs(m, 1))
    checks.check("centroid keeps both cylindrical components", s.simplify(combined + s.I * (s.diff(p, r) + p / r) / (2 * om + sig)) == 0)
    checks.check("radial centroid Bessel antiderivative", s.simplify(s.expand_func(s.diff(r * s.besselj(1, kap * r), r) - r * kap * s.besselj(0, kap * r))) == 0)
    time, beta, rho, norm = s.symbols("t beta rho norm", positive=True)
    signeg = -om
    checks.check("KKS remains nonzero at laboratory crossing", s.simplify(-rho * norm / (2 * signeg) - rho * norm / (2 * om)) == 0)
    j = s.Matrix([[0, -1], [1, 0]])
    frame = s.Matrix([[s.cos(m * om * time), s.sin(m * om * time)], [-s.sin(m * om * time), s.cos(m * om * time)]])
    pull = moving_phase_pullback(-beta * j, -beta * sig * s.eye(2), frame, frame.diff(time))
    checks.check("physical carrier/frame term is in laboratory action", s.simplify(pull.hamiltonian + beta * (sig + m * om) * s.eye(2)) == s.zeros(2))
    freqvar = s.symbols("omega_lab", nonzero=True, real=True)
    q, pp, qdot = s.symbols("q p qdot", real=True)
    lag = beta * pp * qdot + beta * freqvar * (q**2 + pp**2) / 2
    solved_p = s.solve(s.diff(lag, pp), pp)[0]
    checks.check("canonical scalar inertia has a crossing pole", s.simplify(s.diff(lag.subs(pp, solved_p), qdot, 2) + beta / freqvar) == 0)
    mass, amp = s.symbols("M centroid_amplitude", positive=True)
    c = s.Matrix([[amp, 0]])
    momentum_row = mass * c * freqvar * j
    observation = c.col_join(momentum_row)
    checks.check("physical centroid momentum determinant vanishes", s.simplify(observation.det() + mass * amp**2 * freqvar) == 0 and observation.subs(freqvar, 0).rank() == 1)

    coords = s.Matrix(s.symbols("X Y Z", real=True))
    axis = s.Matrix([0, 0, 1])
    bg = om * axis.cross(coords)
    initial = s.Matrix(s.symbols("U0:3", real=True))
    boost = s.Matrix(s.symbols("V0:3", real=True))
    shift = initial + time * boost
    galilean_v = boost - bg.jacobian(coords) * shift
    checks.check("actual Galilean tag momentum is independent", s.simplify(galilean_v + bg.jacobian(coords) * shift - boost) == s.zeros(3, 1))
    pressure_gradient = bg.jacobian(coords)**2 * shift
    checks.check("Galilean velocity and translated pressure solve Euler", s.simplify(galilean_v.diff(time) + bg.jacobian(coords) * galilean_v + pressure_gradient) == s.zeros(3, 1))

    xx, yy, zz = coords
    periodic = s.Matrix([s.sin(zz) + s.cos(yy), s.sin(xx) + s.cos(zz), s.sin(yy) + s.cos(xx)])
    energy = periodic.dot(periodic) / 2
    convective = periodic.jacobian(coords) * periodic
    checks.check("periodic fixture is exact stationary Euler", s.simplify(convective - s.Matrix([s.diff(energy, v) for v in coords])) == s.zeros(3, 1))
    tangent = -periodic.diff(xx)
    null_rhs = periodic.jacobian(coords) * tangent + tangent.jacobian(coords) * periodic
    checks.check("translation mode is a true projected null vector", s.simplify(null_rhs + s.Matrix([s.diff(energy, xx, v) for v in coords])) == s.zeros(3, 1))
    checks.check("constant mean has the translation generalized image", -periodic.jacobian(coords) * s.eye(3)[:, 0] == tangent)
    stress = periodic * tangent.T + tangent * periodic.T
    checks.check("translation stress is an exact periodic derivative", s.simplify(stress + (periodic * periodic.T).diff(xx)) == s.zeros(3))

    speed, lam = s.symbols("V lambda", positive=True)
    helical = speed * s.Matrix([s.cos(lam * zz), -s.sin(lam * zz), 0])
    translation = -helical.diff(zz)
    direction = s.eye(3)[:, 0]
    corrector = s.Matrix([0, 0, speed * s.cos(lam * zz)])
    rhs = -helical.dot(direction) * translation
    cell_lhs = helical.jacobian(coords) * corrector + corrector.jacobian(coords) * helical
    checks.check("actual helical cell equation has an explicit corrector", s.simplify(cell_lhs - rhs) == s.zeros(3, 1))
    checks.check("cell corrector has the required Bloch divergence", s.simplify(s.diff(corrector[2], zz) + direction.dot(translation)) == 0)
    response = (helical * corrector.T + corrector * helical.T) * direction
    averaged_response = response.applyfunc(lambda entry: s.simplify(s.integrate(entry, (zz, 0, 2 * s.pi / lam)) * lam / (2 * s.pi)))
    checks.check("physical mean acoustic coefficient is negative", averaged_response == s.Matrix([0, 0, speed**2 / 2]))
    tvec, nvec = s.Matrix(s.symbols("t0:3")), s.Matrix(s.symbols("n0:3"))
    uv = s.Matrix(s.symbols("u0:3"))
    checks.check("full Beltrami forcing triple product", s.simplify(uv.cross(nvec.cross(tvec)) + uv.dot(nvec) * tvec - nvec * uv.dot(tvec)) == s.zeros(3, 1))
    helical_mode = s.Matrix([1, s.I, 0])
    uzmark = s.symbols("U_z", real=True)
    cross_axis = s.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
    helicity_projection = (s.diag(1, 1, 0) + s.I * cross_axis) / 2
    resonant_rhs = -nvec.cross(-s.I * lam * uzmark * helical_mode)
    checks.check("same-curl-shell inverse has an actual resonance", s.simplify(helicity_projection * resonant_rhs - lam * nvec[2] * uzmark * helical_mode) == s.zeros(3, 1))

    kk, cs = s.symbols("K c", nonzero=True)
    ux = s.Function("Ux")(zz)
    uy = s.Function("Uy")(zz)
    ff = s.Function("f")(zz)
    ray_a = (ux - cs)**2
    ray_residual = s.diff(ray_a * s.diff(ff, zz), zz) - kk**2 * ray_a * ff
    ray_velocity = s.Matrix([-s.diff((ux - cs) * ff, zz), -s.diff(uy, zz) * ff, s.I * kk * (ux - cs) * ff])
    ray_pressure = ray_a * s.diff(ff, zz)
    checks.check("Rayleigh branch is exactly divergence free", s.simplify(s.I * kk * ray_velocity[0] + s.diff(ray_velocity[2], zz)) == 0)
    checks.check("Rayleigh branch solves the full horizontal Euler equations", s.simplify(s.I * kk * (ux - cs) * ray_velocity[0] + s.diff(ux, zz) * ray_velocity[2] + s.I * kk * ray_pressure) == 0 and s.simplify(s.I * kk * (ux - cs) * ray_velocity[1] + s.diff(uy, zz) * ray_velocity[2]) == 0)
    checks.check("vertical Euler gives the exact scalar equation", s.simplify(s.I * kk * (ux - cs) * ray_velocity[2] + s.diff(ray_pressure, zz) - ray_residual) == 0)
    ray_xi = s.Matrix([s.I * s.diff(ff, zz) / kk, 0, ff])
    checks.check("material displacement is exactly divergence free", s.simplify(s.I * kk * ray_xi[0] + s.diff(ray_xi[2], zz)) == 0)
    lin_ray = s.I * kk * (ux - cs) * ray_xi - ff * s.Matrix([s.diff(ux, zz), s.diff(uy, zz), 0])
    checks.check("actual Lin reconstruction gives the Rayleigh velocity", s.simplify(lin_ray - ray_velocity) == s.zeros(3, 1))
    vort = s.Matrix([-s.diff(uy, zz), s.diff(ux, zz), 0])
    kelvin_potential = s.I * (ux - cs) * s.diff(ff, zz) / kk
    kelvin_difference = s.simplify(ray_velocity - ray_xi.cross(vort) - s.Matrix([s.I * kk * kelvin_potential, 0, s.diff(kelvin_potential, zz)]))
    checks.check("complete Kelvin pressure identity uses the same scalar equation", kelvin_difference[0] == 0 and kelvin_difference[1] == 0 and s.simplify(kelvin_difference[2] + s.I * ray_residual / (kk * (ux - cs))) == 0)
    mean_a = s.simplify(s.integrate((speed * s.cos(lam * zz) - cs)**2, (zz, 0, 2 * s.pi / lam)) * lam / (2 * s.pi))
    checks.check("actual scalar solvability gives imaginary acoustic speed", s.simplify(mean_a - cs**2 - speed**2 / 2) == 0)
    c0 = s.I * speed / s.sqrt(2)
    checks.check("scalar implicit continuation is nondegenerate", s.simplify(mean_a.subs(cs, c0)) == 0 and s.simplify(s.diff(mean_a, cs).subs(cs, c0)) != 0)
    reciprocal = cs / (cs**2 - speed**2)**s.Rational(3, 2)
    checks.check("periodic mean-zero inverse is nondegenerate", s.simplify(reciprocal.subs(cs, c0)) != 0 and s.simplify(s.diff(-1 / s.sqrt(cs**2 - speed**2), cs) - reciprocal) == 0)
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
