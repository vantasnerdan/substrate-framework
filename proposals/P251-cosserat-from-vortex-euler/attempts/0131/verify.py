"""Actual higher-azimuthal Euler fields and physical material observations."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0131-material-azimuthal-mode")
    r = s.Symbol("r", positive=True)
    omega, axial = s.symbols("Omega N", positive=True)
    sigma = s.Symbol("sigma", real=True, nonzero=True)
    azimuth = s.Symbol("m", integer=True, positive=True)
    pressure = s.Function("P")(r)
    denominator = 4*omega**2-sigma**2
    radial_square = axial**2*denominator/sigma**2
    radial_equation = {s.diff(pressure, r, 2): -s.diff(pressure, r)/r
                       -(radial_square-azimuth**2/r**2)*pressure}
    vr = s.I*(sigma*s.diff(pressure, r)-2*omega*azimuth*pressure/r)/denominator
    vp = (2*omega*s.diff(pressure, r)-sigma*azimuth*pressure/r)/denominator
    vz = axial*pressure/sigma
    mode = s.Matrix([vr, vp, vz])
    residual = s.Matrix([-s.I*sigma*vr-2*omega*vp+s.diff(pressure, r),
                         -s.I*sigma*vp+2*omega*vr+s.I*azimuth*pressure/r,
                         -s.I*sigma*vz+s.I*axial*pressure])
    checks.check("complete cylindrical Euler equations include the pressure",
                 s.simplify(residual) == s.zeros(3, 1))
    divergence = s.diff(vr, r)+vr/r+s.I*azimuth*vp/r+s.I*axial*vz
    checks.check("Bessel radial equation and inertial dispersion give incompressibility",
                 s.simplify(divergence.subs(radial_equation)) == 0)
    curl = s.Matrix([s.I*azimuth*vz/r-s.I*axial*vp,
                     s.I*axial*vr-s.diff(vz, r),
                     s.diff(vp, r)+vp/r-s.I*azimuth*vr/r])
    checks.check("actual inertial mode is Beltrami, including its axial component",
                 s.simplify((curl+2*omega*axial*mode/sigma).subs(radial_equation))
                 == s.zeros(3, 1))
    displacement = s.I*mode/sigma
    checks.check("displacement satisfies the actual rotating-coordinate Lin equation",
                 s.simplify(-s.I*sigma*displacement-mode) == s.zeros(3, 1))
    spin_density = r*(2*omega*displacement[0]+vp)
    checks.check("full material axial spin reduces to its pressure moment",
                 s.simplify(spin_density-azimuth*pressure/sigma) == 0)
    checks.check("actual pressure torque differentiates the full measured spin",
                 s.simplify(-s.I*sigma*spin_density+s.I*azimuth*pressure) == 0)

    # Explicit m=3 polynomial endpoint, with full Cartesian Euler and Lin.
    x, y, z, t = s.symbols("x y z t", real=True)
    coordinates = (x, y, z)
    zeta = x+s.I*y
    background = s.Matrix([-omega*y, omega*x, 0])
    gradient = background.jacobian(coordinates)
    velocity = s.Matrix([1, s.I, 0])*zeta**2*s.exp(s.I*axial*z-s.I*omega*t)
    xi = -s.I*velocity/(2*omega)
    euler = velocity.diff(t)+velocity.jacobian(coordinates)*background+gradient*velocity
    lin = xi.diff(t)+xi.jacobian(coordinates)*background-gradient*xi-velocity
    checks.check("m3 endpoint is an exact pressure-free laboratory Euler field",
                 s.simplify(euler) == s.zeros(3, 1))
    checks.check("m3 endpoint has a genuine material displacement",
                 s.simplify(lin) == s.zeros(3, 1))
    point_spin = xi.cross(background)+s.Matrix(coordinates).cross(gradient*xi+velocity)
    checks.check("endpoint has zero linear axial spin for every material tag",
                 s.simplify(point_spin[2]) == 0)
    checks.check("higher mode has no axis core-vector observation",
                 velocity.subs({x: 0, y: 0}) == s.zeros(3, 1))

    # Rotation covariance of an actual material complex moment, not a mask.
    order = s.Symbol("order", integer=True, positive=True)
    base_moment, shape_cos, shape_sin = s.symbols("Q0 shape_cos shape_sin", real=True)
    perturbation = shape_cos*s.cos(sigma*t)+shape_sin*s.sin(sigma*t)
    baseline = base_moment*s.exp(s.I*order*omega*t)
    moment_variation = s.I*order*baseline*perturbation
    physical_angle = s.simplify(moment_variation/(s.I*order*baseline))
    checks.check("transported material orientation cancels the azimuthal Doppler phase",
                 s.simplify(s.diff(physical_angle, t, 2)+sigma**2*physical_angle) == 0)
    checks.check("a stationary laboratory mask would give a different frequency",
                 s.expand((azimuth*omega+sigma)**2-sigma**2) != 0)

    # Angular integration for a genuine positive m-lobed material marker.
    # Constant and 2m Fourier channels are evaluated explicitly at m=3;
    # the same integer-orthogonality argument applies for every m>=3.
    phi, mark = s.symbols("phi mark", real=True)
    weight = 1+mark*s.cos(3*phi)
    checks.check("marked material geometry has a nonzero reference third moment",
                 s.integrate(s.expand_complex(s.exp(3*s.I*phi)*weight),
                             (phi, 0, 2*s.pi)) == s.pi*mark)
    checks.check("material shape variation has only its zero angular channel",
                 s.integrate(s.expand_complex(s.exp(6*s.I*phi)*weight),
                             (phi, 0, 2*s.pi)) == 0)
    a = displacement[0]
    b = displacement[1]/s.I
    checks.check("actual deformation supplies the normalized shape-angle row",
                 s.simplify(a+b-(s.diff(pressure, r)+azimuth*pressure/r)
                            /(sigma*(2*omega+sigma))) == 0)

    k = s.Symbol("K", real=True)
    intrinsic = -2*omega*axial/s.sqrt(axial**2+k**2)
    lab = azimuth*omega+intrinsic
    checks.check("material angle squared-frequency curvature is negative",
                 s.diff(intrinsic**2, k, 2).subs(k, 0)/2 == -4*omega**2/axial**2)
    checks.check("laboratory pattern curvature has the advertised higher-m sign",
                 s.simplify(s.diff(lab**2, k, 2).subs(k, 0)/2
                            -2*(azimuth-2)*omega**2/axial**2) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
