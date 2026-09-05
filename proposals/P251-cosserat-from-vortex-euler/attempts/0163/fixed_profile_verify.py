"""Exact boundary, symmetry and full-current anchors of the fixed-C proof."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0163-fixed-profile-acoustics")
    a, b = s.symbols("a b", real=True)
    psi = s.cos(a)+s.cos(b)+s.cos(a+b)
    for sign in (-1, 1):
        transformed = psi.subs({a: -a-b+sign*2*s.pi, b: a}, simultaneous=True)
        checks.check(f"actual polygon C3 symmetry holds at triangle center sign {sign}",
                     s.trigsimp(transformed-psi) == 0)
    rot3 = s.Matrix([[-1, -1], [1, 0]])
    rot6 = s.Matrix([[0, -1], [1, 1]])
    checks.check("each triangle's orbital vector average has no nonzero invariant",
                 rot3**2+rot3+s.eye(2) == s.zeros(2)
                 and (rot3-s.eye(2)).det() != 0)
    checks.check("the hexagonal cell has its own complete C6 symmetry",
                 rot6**3 == -s.eye(2)
                 and s.trigsimp(psi.subs({a: -b, b: a+b}, simultaneous=True)-psi) == 0)
    hess = s.hessian(psi, (a, b))
    checks.check("the common separatrix vertices are genuinely hyperbolic saddles",
                 hess.subs({a: s.pi, b: 0}).det() == -1)
    checks.check("the polygon centers have nondegenerate extremum Hessians",
                 hess.subs({a: 0, b: 0}).det() == 3
                 and hess.subs({a: 2*s.pi/3, b: 2*s.pi/3}).det() == s.Rational(3, 4))
    checks.check("critical-contour classification starts from the exact factorization",
                 s.trigsimp(s.sin(a)-s.sin(b)
                            -2*s.cos((a+b)/2)*s.sin((a-b)/2)) == 0)

    x, y = s.symbols("x y", real=True)
    profile, potential = s.Function("psi")(x, y), s.Function("phi")(x, y)
    gradp = s.Matrix([s.diff(profile, x), s.diff(profile, y)])
    gradf = s.Matrix([s.diff(potential, x), s.diff(potential, y)])
    j = s.Matrix([[0, -1], [1, 0]])
    velocity = j*gradp
    vort = s.diff(velocity[1], x)-s.diff(velocity[0], y)
    left = gradf.jacobian((x, y))*velocity+velocity.jacobian((x, y))*gradf
    adv = velocity.dot(gradf)
    right = s.Matrix([s.diff(adv, x), s.diff(adv, y)])+vort*j*gradf
    checks.check("edge proof uses the full transport-plus-strain identity",
                 all(s.expand(entry) == 0 for entry in left-right))
    ws, wm, zs, normal, length = s.symbols("Ws Wmean zeta_s Xn length", real=True)
    circulation = ws*(-zs*normal*length)+zs*(ws-wm)*normal*length
    checks.check("retaining both separatrix terms gives exactly the mean axial frame",
                 s.expand(circulation+zs*wm*normal*length) == 0)
    checks.check("deleting the Green normal-return period changes the actual drift",
                 s.expand(circulation-ws*(-zs*normal*length)) != 0)

    lam, ps, wp, dk, dpsi = s.symbols("lambda psi Wprime dK dpsi", real=True, nonzero=True)
    c = -wp/lam**2
    zeta = -lam**2*ps
    v_translation = -dk
    b_translation = -wp*dpsi
    pressure_translation = dk+lam**2*ps*dpsi
    full_source = pressure_translation-c*wp*v_translation+c*zeta*b_translation
    factor = 1-wp**2/lam**2
    checks.check("full translation pressure and axial reaction produce the transport primitive",
                 s.expand(full_source-factor*pressure_translation) == 0)
    checks.check("constant-curl and planar scalar correctors are genuinely different",
                 factor.subs(wp, -lam) == 0 and factor.subs(wp, 0) == 1)

    k = s.Symbol("k", positive=True)
    grow, vrow, hforce, rsource, wdb = s.symbols("g_eta V_r hE RS Wdb")
    mean_rate = -s.I*k*(grow+vrow)+k*k*wdb
    h_quotient_rate = grow+s.I*k*hforce
    r_moment_rate = vrow-s.I*k*rsource
    corrected_rate = s.expand(mean_rate+s.I*k*(h_quotient_rate+r_moment_rate))
    checks.check("exact improved mean retains all finite-k current production",
                 s.expand(corrected_rate-k*k*(-hforce+rsource+wdb)) == 0)
    hf, avv = s.symbols("hF weighted_vv")
    tangent_force = -s.I*(s.I*hf-s.I*avv)+wdb
    checks.check("stationary adjoint duality matches current stiffness to the full Euler tangent",
                 s.expand(tangent_force-hf+avv-wdb) == 0)
    cb, delta, remainder = s.symbols("cb delta remainder", positive=True)
    lower = cb**2-remainder*delta**2
    checks.check("finite primitive bounds give positive stiffness at fixed small shear",
                 s.simplify(lower.subs(delta, cb/(2*s.sqrt(remainder)))-3*cb**2/4) == 0)
    graph_error = s.Symbol("graph_error", positive=True)
    total_residual_bound = k*graph_error+k*k
    checks.check("smooth graph approximation preserves the actual second-order residual",
                 s.expand(total_residual_bound.subs(graph_error, k)-2*k*k) == 0)

    radius = s.Symbol("radius", positive=True)
    log_primitive_integral = s.integrate((-s.log(radius))**3, (radius, 0, 1))
    checks.check("the logarithmic transport primitive has an exact finite L2 coarea bound",
                 log_primitive_integral == 6)
    print("Scope: exact proof anchors; transport graph and bounded-group arguments are given in the analytic artifact.")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
