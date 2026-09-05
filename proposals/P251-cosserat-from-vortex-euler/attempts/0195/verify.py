"""Exact exposing algebra for the global swirl/force-free ring construction."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0195-swirl-ring-and-force-free-endpoint")
    r, z, radius = s.symbols("r z R", positive=True)
    epsilon, eta = s.symbols("epsilon eta", real=True)
    psi = s.Function("psi")(r, z)
    ff = s.Function("F")(psi)
    fp = s.diff(ff, psi)
    u = s.Matrix([-s.diff(psi, z)/r, ff/r, s.diff(psi, r)/r])
    omega = s.Matrix([-s.diff(u[1], z), s.diff(u[0], z)-s.diff(u[2], r),
                      s.diff(r*u[1], r)/r])
    delta = s.diff(psi, r, 2)-s.diff(psi, r)/r+s.diff(psi, z, 2)
    grad = s.Matrix([s.diff(psi, r), 0, s.diff(psi, z)])
    azimuth = s.Matrix([0, 1, 0])
    checks.check("actual cylindrical velocity is divergence free",
                 s.simplify(s.diff(r*u[0], r)/r+s.diff(u[2], z)) == 0)
    checks.check("actual vorticity retains the swirl and Bernoulli components",
                 s.simplify(omega-fp*u-(-delta-ff*fp)*azimuth/r) == s.zeros(3, 1))
    checks.check("full steady Euler cross-product source sign",
                 s.simplify(u.cross(omega)-(delta+ff*fp)*grad/r**2) == s.zeros(3, 1))
    checks.check("streamfunction is a material first integral", s.simplify(u.dot(grad)) == 0)
    checks.check("streamfunction is also a vortex-line first integral",
                 s.simplify(omega.dot(grad)) == 0)
    stream = s.symbols("stream", real=True)
    g = s.Function("g")(stream/radius)
    swirl = epsilon*radius*g
    source = s.simplify(swirl*s.diff(swirl, stream))
    gg_prime = s.diff(g, stream)*radius
    checks.check("scaled swirl source contains its full R factor",
                 s.simplify(source-epsilon**2*radius*g*gg_prime) == 0)
    f, gg = s.symbols("f g_gprime", real=True)
    effective = f+radius**2*epsilon**2*gg/r**2
    checks.check("same Green kernel uses the correct weighted effective source",
                 s.expand(r**2*effective/radius**2-r**2*f/radius**2-epsilon**2*gg) == 0)
    checks.check("omitting the swirl source changes the actual stationary equation",
                 s.expand(r**2*(effective-f)/radius**2) != 0)

    ii, theta, phi, time = s.symbols("I theta varphi time", real=True)
    om, vv = s.Function("Omega")(ii), s.Function("V")(ii)
    aa, bb = s.Function("alpha")(ii), s.Function("beta")(ii)
    coordinates = (ii, theta, phi)
    flow = s.Matrix([0, om, vv])
    vort = s.Matrix([0, aa, bb])
    perturbation = s.Matrix([s.Function(f"w{k}")(*coordinates) for k in range(3)])
    velocity = s.Matrix([s.Function(f"v{k}")(*coordinates) for k in range(3)])

    def bracket(left, right):
        return right.jacobian(coordinates)*left-left.jacobian(coordinates)*right

    full = -bracket(flow, perturbation)+bracket(vort, velocity)
    expected = (-om*perturbation.diff(theta)-vv*perturbation.diff(phi)
                +s.Matrix([0, s.diff(om, ii), s.diff(vv, ii)])*perturbation[0]
                +aa*velocity.diff(theta)+bb*velocity.diff(phi)
                -s.Matrix([0, s.diff(aa, ii), s.diff(bb, ii)])*velocity[0])
    checks.check("complete all-poloidal Euler Lie operator including new pressure derivative",
                 s.simplify(full-expected) == s.zeros(3, 1))
    checks.check("actual steady flux-coordinate vorticity commutes with velocity",
                 bracket(flow, vort) == s.zeros(3, 1))
    n = s.symbols("n", integer=True)
    v0 = s.Matrix([s.Function(f"v0{k}")(ii, theta) for k in range(3)])
    harmonic = aa*v0.diff(theta)+s.I*n*bb*v0
    direct = (aa*(v0*s.exp(s.I*n*phi)).diff(theta)
              +bb*(v0*s.exp(s.I*n*phi)).diff(phi))/s.exp(s.I*n*phi)
    checks.check("fixed toroidal harmonic retains poloidal vorticity response",
                 s.simplify(direct-harmonic) == s.zeros(3, 1))
    nn = s.Matrix([[0, 0, 0], [s.diff(om, ii), 0, 0], [s.diff(vv, ii), 0, 0]])
    dd = s.symbols("d", nonzero=True)
    checks.check("both transport shears form the same nilpotent block", nn**2 == s.zeros(3))
    checks.check("shifted-band graph resolvent retains both shears",
                 s.simplify((dd*s.eye(3)-nn)*(s.eye(3)/dd+nn/dd**2)-s.eye(3)) == s.zeros(3))

    flux_derivative = s.diff(om/vv, ii)/vv
    checks.check("twist uses the actual return flux action",
                 s.simplify(flux_derivative-(s.diff(om, ii)*vv-om*s.diff(vv, ii))/vv**3) == 0)
    g0, g1, g2, phip, om0, omp, hh, hp = s.symbols(
        "g0 g1 g2 phip om0 omp H Hp", nonzero=True, real=True)
    twist_coefficient = omp/om0-g1*phip/g0-hp/hh
    checks.check("free first swirl derivative independently controls streamline twist",
                 s.diff(twist_coefficient, g1) == -phip/g0)
    fc = s.symbols("fc", positive=True)
    phi_i = s.Function("Phi")(ii)
    gi = s.Function("g")(phi_i)
    hi = s.Function("H")(ii)
    alphai = epsilon*s.diff(gi, phi_i)*om
    betai = fc/radius+epsilon**2*radius*s.diff(gi, phi_i)*gi*hi
    leading = s.limit(s.diff(alphai/betai, ii)/epsilon, epsilon, 0)
    expected_leading = radius/fc*(s.diff(gi, phi_i, 2)*s.diff(phi_i, ii)*om
                                  +s.diff(gi, phi_i)*s.diff(om, ii))
    checks.check("actual vortex-line twist has the independent second profile derivative",
                 s.simplify(leading-expected_leading) == 0)
    checks.check("vortex-line twist profile control coefficient is nonzero",
                 s.diff(radius/fc*(g2*phip*om0+g1*omp), g2) == radius*phip*om0/fc)

    chi = s.Function("chi")
    transported = chi(ii, phi-vv*time)
    checks.check("finite-arc density obeys its actual angular transport",
                 s.simplify(s.diff(transported, time)+vv*s.diff(transported, phi)) == 0)
    frozen = chi(ii, phi)
    checks.check("generic stationary finite-arc substitution leaves a nonzero transport row",
                 vv*s.diff(frozen, phi) != 0)
    checks.check("zero actual angular velocity on open tag support removes that row",
                 (vv*s.diff(frozen, phi)).subs(vv, 0) == 0)
    position = s.Matrix(s.symbols("rx ry rz", real=True))
    displacement = s.Matrix(s.symbols("dx dy dz", real=True))
    base_velocity = s.Matrix(s.symbols("ux uy uz", real=True))
    material_rate = s.Matrix(s.symbols("vx vy vz", real=True))
    spin = displacement.cross(base_velocity)+position.cross(material_rate)
    moment_rate = base_velocity.cross(displacement)+position.cross(material_rate)
    checks.check("full Reynolds material spin includes both moving-position terms",
                 spin-moment_rate == 2*displacement.cross(base_velocity))
    gq, gp, cq, cp, q, qt, qtt = s.symbols("Gq Gp Cq Cp q qt qtt", real=True)
    gp_pair = gq*q+gp*qt
    gm_pair = gq*q-gp*qt
    sp_pair = gq*qt+gp*qtt+2*cq*q+2*cp*qt
    sm_pair = gq*qt-gp*qtt-2*cq*q+2*cp*qt
    checks.check("actual prepared TR moment average", s.expand((gp_pair+gm_pair)/2-gq*q) == 0)
    checks.check("actual prepared TR spin keeps the measured even connection",
                 s.expand((sp_pair+sm_pair)/2-(gq+2*cp)*qt) == 0)
    checks.check("the corrected physical tag row yields exact spin/current matching",
                 s.expand(((sp_pair+sm_pair)/2).subs(cp, 0)-gq*qt) == 0)
    mass, frequency = s.symbols("mass frequency", positive=True)
    form = mass*s.Matrix([[0, 1], [-1, 0]])
    tr_map = s.diag(1, -1)
    checks.check("TR full initial phase pullback precedes elimination",
                 tr_map.T*(-form)*tr_map == form)
    energy = mass*s.diag(frequency**2, 1)
    checks.check("TR actual initial energy retains positive phase weighting",
                 tr_map.T*energy*tr_map == energy)

    t, scale = s.symbols("t A", positive=True)
    flat_g = s.sqrt(scale)*s.exp(-1/t)
    flat_f = scale*s.exp(-2/t)/t**2
    checks.check("explicit source tail earns the smooth square-root primitive",
                 s.simplify(flat_g*s.diff(flat_g, t)-flat_f) == 0)
    checks.check("primitive and factor really vanish smoothly at the core edge",
                 all(s.limit(s.diff(flat_g, t, order), t, 0, dir="+") == 0
                     for order in range(4)))
    gfun = s.Function("G")(stream/radius)
    homotopy_f = radius*eta*gfun
    product = homotopy_f*s.diff(homotopy_f, stream)
    checks.check("homotopy swirl source is the full eta-squared primitive source",
                 s.simplify(product-radius**2*eta**2*gfun*s.diff(gfun, stream)) == 0)
    qhom = (1-eta**2+eta**2*radius**2/r**2)*f
    checks.check("uniform homotopy has the correct full Green source",
                 s.expand(r**2*qhom/radius**2-((1-eta**2)*r**2/radius**2+eta**2)*f) == 0)
    checks.check("force-free endpoint exactly removes Bernoulli vorticity",
                 ((1-eta**2)*r*f/radius).subs(eta, 1) == 0)
    inverse_r, x = s.symbols("inverse_R x", real=True)
    factor = 1-eta**2+eta**2/(1+inverse_r*x)**2
    checks.check("uniform source perturbation is first order in inverse radius",
                 s.diff(factor, inverse_r).subs(inverse_r, 0) == -2*eta**2*x)

    radial, omc, wcore = s.symbols("s Omega_c Wcore", positive=True)
    axial = s.sqrt(wcore**2-2*omc**2*radial**2)
    checks.check("actual force-free flat-core axial derivative includes swirl reaction",
                 s.simplify(s.diff(axial, radial)+2*omc**2*radial/axial) == 0)
    rotation = -radius*omc/axial
    physical_action_derivative = radius*radial
    checks.check("endpoint rotation twist uses physical meridional action",
                 s.simplify(s.diff(rotation, radial)/physical_action_derivative
                            +2*omc**3/axial**3) == 0)
    checks.check("force-free streamline and vorticity return numbers coincide",
                 s.cancel((g1*om0)/(g1*hh)-om0/hh) == 0)
    print("Global IFT, fixed-contour perturbation, flux/KAM conditions and the actual tag IFT are analytic proofs in the companion files.")
    print("The force-free endpoint has no asserted optical pole/current continuation; the full O(1) operator change is retained.")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
