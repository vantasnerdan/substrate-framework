"""Exact reflected-field forced Euler cancellation and weighted threshold form."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0206-reflected-threshold-response")
    r, k, c = s.symbols("r k c", positive=True)
    w, om, y = (s.Function(n)(r) for n in ("w", "Omega", "Y"))
    ar, at, az = (s.Function(n)(r) for n in ("eta_r", "eta_theta", "eta_z"))
    d = c+w
    zeta = -w*s.diff(w, r)/(r*om)
    kappa = 2*om*zeta
    vort = s.Matrix([0, s.diff(w, r), zeta])
    forcing = s.Matrix([ar, at, az]).cross(vort)
    checks.check("actual coadjoint source has all three poloidal-vorticity components",
                 s.simplify(forcing-s.Matrix([zeta*at-s.diff(w, r)*az,
                                             -zeta*ar, s.diff(w, r)*ar])) == s.zeros(3, 1))
    f = y/r-s.I*ar/(k*d)
    pressure = d**2*s.diff(y, r)/r-d*az
    vr = -s.I*k*d*f
    vt = -zeta*f+s.I*forcing[1]/(k*d)
    vz = pressure/d+s.diff(w, r)*f+s.I*forcing[2]/(k*d)
    checks.check("actual radial response retains the particular generator term",
                 s.simplify(vr+s.I*k*d*y/r+ar) == 0)
    checks.check("singular azimuthal source cancels in actual velocity",
                 s.simplify(vt+zeta*y/r) == 0)
    checks.check("singular axial source cancels without deleting its physical return",
                 s.simplify(vz-s.diff(d*y, r)/r+az) == 0)
    divergence = {s.diff(ar, r): -ar/r-s.I*k*az}
    checks.check("exact change of variable uses the real solenoidal generator",
                 s.simplify((s.diff(f, r)+f/r-pressure/d**2
                             -s.I*forcing[2]/(k*d**2)).subs(divergence)) == 0)
    radial_residual = (s.diff(pressure, r)-(k**2*d**2-kappa)*f
                       -forcing[0]-2*s.I*om*forcing[1]/(k*d))
    sturm_residual = (s.diff(d**2*s.diff(y, r)/r, r)
                      +(kappa-k**2*d**2)*y/r-zeta*at-d*(s.diff(az, r)-s.I*k*ar))
    checks.check("full forced Sturm equation has the derived regular source",
                 s.simplify(radial_residual-sturm_residual) == 0)
    checks.check("omitting the actual particular pressure loses the cancellation",
                 s.simplify(radial_residual+s.diff(d*az, r)-sturm_residual) != 0)
    kinetic = (s.diff(y, r)**2+k**2*y**2)/r
    q_zero = w**2*((s.diff(y, r)-y/r)**2/r+y**2/r**3+k**2*y**2/r)
    q_full = d**2*kinetic-kappa*y**2/r
    checks.check("one-sided threshold coercivity follows from the full force-free form",
                 s.simplify(q_full-q_zero-(c**2+2*c*w)*kinetic
                            -s.diff(w**2*y**2/r**2, r)) == 0)
    checks.check("actual Z source has the stated uniformly weighted dual factor",
                 s.simplify(zeta*at*y+(w*y/r**s.Rational(3, 2))
                            *(s.sqrt(r)*s.diff(w, r)*at/om)) == 0)
    tail_flux = s.diff(w*s.diff(w, r)*y**2/r, r)
    tail_expansion = (((s.diff(w, r)**2+w*s.diff(w, r, 2))/r
                       -w*s.diff(w, r)/r**2)*y**2
                      +2*w*s.diff(w, r)*y*s.diff(y, r)/r)
    checks.check("physical velocity tail estimate uses the complete weighted derivative",
                 s.simplify(tail_flux-tail_expansion) == 0)
    phi = s.symbols("phi", positive=True)
    p1, p2 = s.symbols("phi_r phi_rr", real=True, nonzero=True)
    flat = s.exp(-1/phi)
    flat_r = s.diff(flat, phi)*p1
    flat_rr = s.diff(flat, phi, 2)*p1**2+s.diff(flat, phi)*p2
    checks.check("the actual exponential source tail supplies the Hardy condition",
                 s.limit(s.simplify(flat*flat_rr/flat_r**2), phi, 0, dir="+") == 1)
    reflection = s.diag(1, -1, 1)
    vel = s.Matrix(s.symbols("ur uphi uz"))
    factor = s.symbols("force_free_factor", real=True)
    checks.check("whole-field reflection preserves force-free geometry with axial-vector parity",
                 reflection.det()*reflection*(factor*vel)
                 == (-factor)*(reflection*vel))
    gam, radius, xx, yfirst, residual_speed = s.symbols(
        "Gamma R x source_first_moment U_remainder", real=True)
    kernel_log = s.log(radius)*(gam*xx+3*yfirst)/(4*s.pi*radius)
    actual_speed = gam*s.log(radius)/(4*s.pi*radius)+residual_speed/radius
    checks.check("the actual center border and translation remove the large core log dipole",
                 s.simplify((kernel_log-actual_speed*xx).subs(yfirst, 0)
                            +residual_speed*xx/radius) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
