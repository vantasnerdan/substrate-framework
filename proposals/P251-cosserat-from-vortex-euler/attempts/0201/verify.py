"""Exact identities for the actual force-free mode and measured section current.

The analytic compactness, principal-root and ring-contour proofs live in the
companions. These checks expose their equations, signs and observable maps;
they are not a numerical eigenvalue certification.
"""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0201-force-free-mode-and-section-current")
    r, k, c, rho = s.symbols("r k c rho", positive=True)
    m, omega = s.symbols("m omega", real=True)
    om, w = s.Function("Omega")(r), s.Function("W")(r)
    f, pressure = s.Function("f")(r), s.Function("P")(r)
    zeta = 2*om+r*s.diff(om, r)
    kappa = 2*om*zeta
    doppler = omega-m*om-k*w
    vr = -s.I*doppler*f
    vt = m*pressure/(r*doppler)-zeta*f
    vz = k*pressure/doppler-s.diff(w, r)*f
    fp = -(1/r+2*m*om/(r*doppler))*f+(m**2/r**2+k**2)*pressure/doppler**2
    pp = (doppler**2-kappa)*f+2*m*om*pressure/(r*doppler)
    replace = {s.diff(f, r): fp, s.diff(pressure, r): pp}
    checks.check("full radial momentum including pressure", s.simplify(
        (-s.I*doppler*vr-2*om*vt+s.diff(pressure, r)).subs(replace)) == 0)
    checks.check("full azimuthal momentum", s.simplify(
        -s.I*doppler*vt+zeta*vr+s.I*m*pressure/r) == 0)
    checks.check("full axial momentum retains W prime", s.simplify(
        -s.I*doppler*vz+s.diff(w, r)*vr+s.I*k*pressure) == 0)
    checks.check("full three-dimensional incompressibility", s.simplify(
        (s.diff(r*vr, r)/r+s.I*m*vt/r+s.I*k*vz).subs(replace)) == 0)
    checks.check("omitting axial shear is exposed by actual Euler momentum", s.simplify(
        -s.I*doppler*(k*pressure/doppler)+s.diff(w, r)*vr+s.I*k*pressure) != 0)

    y = s.Function("y")(r)
    d = c-w
    e = (s.diff(y, r)**2+k**2*y**2)/r
    p0 = d**2*s.diff(y, r)/r
    sturm = s.diff(p0, r)+(kappa-k**2*d**2)*y/r
    checks.check("axisymmetric pressure elimination is the stated Sturm equation",
                 s.simplify(sturm-(s.diff(p0, r)-(k**2*d**2-kappa)*y/r)) == 0)
    xi_z = s.I*s.diff(y, r)/(k*r)
    checks.check("actual displacement is divergence free, not only velocity",
                 s.simplify(s.diff(y, r)/r+s.I*k*xi_z) == 0)
    tail = r*s.besselk(1, k*r)
    checks.check("actual exterior tail solves its full radial equation",
                 s.simplify(s.expand_func(s.diff(tail, r, 2)
                                         -s.diff(tail, r)/r-k**2*tail)) == 0)
    checks.check("Hardy completion retains its exact boundary derivative", s.simplify(
        s.diff(y, r)**2/r-(s.diff(y, r)-y/r)**2/r-y**2/r**3
        -s.diff(y**2/r**2, r)) == 0)
    forcefree_kappa = -2*w*s.diff(w, r)/r
    c0_form = w**2*e-forcefree_kappa*y**2/r
    c0_square = w**2*((s.diff(y, r)-y/r)**2/r+y**2/r**3+k**2*y**2/r)
    checks.check("failed zero/negative-phase trial has an exact positive form",
                 s.simplify(c0_form-c0_square-s.diff(w**2*y**2/r**2, r)) == 0)

    wronskian = s.diff(r*f*pressure, r).subs(replace)
    checks.check("all-poloidal real-frequency energy identity", s.simplify(
        wronskian-r*(m**2/r**2+k**2)*pressure**2/doppler**2
        -r*(doppler**2-kappa)*f**2) == 0)
    matrix = s.Matrix([
        [-(1/r+2*m*om/(r*doppler)), (m**2/r**2+k**2)/doppler**2],
        [doppler**2-kappa, 2*m*om/(r*doppler)],
    ])
    trans = s.Matrix([1, -r*om**2])
    matrix0 = matrix.subs({m: 1, omega: 0, k: 0})
    checks.check("exact Euclidean translation solves the entire radial system",
                 s.simplify(trans.diff(r)-matrix0*trans) == s.zeros(2, 1))
    perturb = matrix.subs({m: 1, omega: k*c}).diff(k).subs(k, 0)
    forced = perturb*trans
    checks.check("complete W-dependent translation first-order Wronskian cancels",
                 s.simplify(r*(forced[1]-trans[1]*forced[0])) == 0)
    outside = -r*(k*c-om)**2
    checks.check("actual exterior derivative selects zero translation group speed",
                 s.simplify(s.diff(outside, k).subs(k, 0)-2*r*om*c) == 0)
    checks.check("translation derivative is sensitive to the exterior normalization",
                 s.diff(outside, k).subs(k, 0) != 0)

    beta_density = -kappa*y**2/(r*d)+s.diff(w, r)*y*s.diff(y, r)/r
    kappa_on_shell = k**2*d**2-r*s.diff(d**2*s.diff(y, r)/r, r)/y
    total = s.diff(d*y*s.diff(y, r)/r, r)
    checks.check("full KKS reduces to a negative positive-weight kinetic form",
                 s.simplify(beta_density.subs(kappa, kappa_on_shell)+d*e-total) == 0)
    checks.check("dropping W prime destroys the inherited KKS identity", s.simplify(
        (-kappa*y**2/(r*d)).subs(kappa, kappa_on_shell)+d*e-total) != 0)
    bcol = -2*om*(y/r)/(k*d)
    ccol = s.diff(y, r)/(k*r)
    checks.check("velocity and displacement KKS rows agree including poloidal vorticity",
                 s.simplify(k*(y/r)*(zeta*bcol+s.diff(w, r)*ccol)*r
                            -beta_density) == 0)

    # Independent conservation derivation, without an on-shell Sturm substitution.
    time, axial = s.symbols("t z", real=True)
    rr = s.Function("xi_r")(r, axial, time)
    tt = s.Function("xi_theta")(r, axial, time)
    zz = s.Function("xi_z")(r, axial, time)
    lin_r = s.diff(rr, time)+w*s.diff(rr, axial)
    lin_t = s.diff(tt, time)+w*s.diff(tt, axial)-r*s.diff(om, r)*rr
    lin_z = s.diff(zz, time)+w*s.diff(zz, axial)-s.diff(w, r)*rr
    bg_l = r**2*om
    torsion = r*(-zeta*rr)+s.diff(bg_l, r)*rr
    checks.check("literal material particle axial spin is exactly zero",
                 s.simplify(torsion) == 0)
    div_xi = s.diff(r*rr, r)/r+s.diff(zz, axial)
    current = r*lin_t-r*s.diff(tt, time)-s.diff(r*w*tt, axial)
    checks.check("literal displacement moment has the radial shape correction",
                 s.simplify(current+r**2*s.diff(om, r)*rr) == 0)
    div_lxi = s.diff(r*bg_l*rr, r)/r+s.diff(bg_l*zz, axial)
    checks.check("Eulerian spin is the full displacement transport of particle momentum",
                 s.simplify(-r*zeta*rr+div_lxi-bg_l*div_xi) == 0)
    checks.check("the section source integral equals the actual axial displacement flux",
                 s.simplify(r*(r*om)*s.diff(y, r)+r*zeta*y
                            -s.diff(r**2*om*y, r)) == 0)
    div_v = s.diff(r*lin_r, r)/r+s.diff(lin_z, axial)
    conserved = (s.diff(r*lin_t, time)+s.diff(r*bg_l*lin_r, r)/r
                 +s.diff(w*r*lin_t+bg_l*lin_z, axial))
    theta_euler = r*(s.diff(lin_t, time)+w*s.diff(lin_t, axial)+zeta*lin_r)
    checks.check("conservative angular momentum includes all axial/radial returns",
                 s.simplify(conserved-theta_euler-bg_l*div_v) == 0)

    # Direct angular integrations of the NONNEGATIVE material marker.
    angle, omc, eps = s.symbols("theta Omega_c epsilon", real=True)
    chi, chip = s.symbols("chi chi_prime", real=True)
    xr, xt = s.symbols("xi_r xi_theta", real=True)
    mark = chi*(1+eps*s.cos(2*angle))
    delta_mark = -xr*chip*(1+eps*s.cos(2*angle))-xt*s.diff(mark, angle)/r
    qbase = s.integrate(s.expand_trig(s.cos(2*angle)*mark), (angle, 0, 2*s.pi))*r**3
    qreal = s.integrate(s.expand_trig(s.cos(2*angle)*delta_mark), (angle, 0, 2*s.pi))*r**3
    qimag = s.integrate(s.expand_trig(s.sin(2*angle)*delta_mark), (angle, 0, 2*s.pi))*r**3
    checks.check("actual painted quadrupole has a nonzero base",
                 s.simplify(qbase-s.pi*eps*chi*r**3) == 0)
    checks.check("radial marker transport changes the actual quadrupole magnitude",
                 s.simplify(qreal+s.pi*eps*chip*xr*r**3) == 0)
    checks.check("actual angular integration gives the physical angle row",
                 s.simplify(qimag-2*s.pi*eps*chi*r**2*xt) == 0)
    checks.check("physical angle has unit response to rigid rotation",
                 s.simplify(qimag.subs(xt, r*angle)/(2*qbase)-angle) == 0)
    transported_mark = s.Function("chi")(r)*(1+eps*s.cos(2*(angle-omc*time)))
    checks.check("base marker obeys actual passive transport, including arbitrary axial W",
                 s.simplify(s.diff(transported_mark, time)
                            +omc*s.diff(transported_mark, angle)
                            +w*s.diff(transported_mark, axial)) == 0)
    cc, ss, freq = s.symbols("C_theta S0 omega_positive", positive=True)
    phase = k*axial-freq*time
    physical_angle = cc*s.sin(phase)
    physical_spin = -ss*s.cos(phase)
    checks.check("literal section spin has positive measured angle-rate overlap",
                 s.simplify(physical_spin-ss/(freq*cc)*s.diff(physical_angle, time)) == 0)
    checks.check("opposite spin sign is detected by physical time differentiation",
                 s.simplify(physical_spin+ss/(freq*cc)*s.diff(physical_angle, time)) != 0)

    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
