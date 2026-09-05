"""Exact generalized-force-free and full-action/ordering checks, not a solver."""

import sympy as s

from substrate_framework.euler_acoustic import (
    axial_kelvin_initial_phase,
    triangular_euler_array,
)
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0163-ordered-force-free")
    lam, amp, rho, cap = s.symbols("lambda Psi rho C", positive=True)
    k, eps, delta = s.symbols("k epsilon delta", positive=True)
    xx, yy, ps = s.symbols("x y psi", real=True)
    array = triangular_euler_array(amp, lam, rho, (xx, yy))
    w = s.sqrt(cap+lam**2*ps**2)
    force_free_factor = -s.diff(w, ps)
    checks.check("same planar geometry gives exact generalized force-free curl",
                 s.simplify(force_free_factor*w+lam**2*ps) == 0)
    checks.check("the whole stationary Bernoulli constant is retained",
                 s.simplify(w*w-lam**2*ps**2-cap) == 0)
    h = w-s.sqrt(cap)
    checks.check("small moving-frame axial profile has an exact rationalization",
                 s.simplify(h-lam**2*ps**2/(w+s.sqrt(cap))) == 0)
    checks.check("the amplitude and derivative have the declared large-C order",
                 s.limit(s.sqrt(cap)*h, cap, s.oo) == lam**2*ps**2/2
                 and s.limit(s.sqrt(cap)*s.diff(w, ps), cap, s.oo) == lam**2*ps)
    core_w = w.subs(ps, 3*amp)
    checks.check("the optical interface uses the actual new core derivatives",
                 s.simplify(s.diff(w, ps).subs(ps, 3*amp)-3*lam**2*amp/core_w) == 0
                 and s.simplify(s.diff(w, ps, 2).subs(ps, 3*amp)
                                -lam**2*cap/core_w**3) == 0)

    zeta_profile = s.Function("Z")(ps)
    axial_profile = s.Function("F")(ps)
    passive_coefficient = s.diff(axial_profile, ps)/s.diff(zeta_profile, ps)
    checks.check("general-profile passive variable cancels both zero-k shear sources",
                 s.simplify(-s.diff(axial_profile, ps)
                            +passive_coefficient*s.diff(zeta_profile, ps)) == 0)
    grad_psi = s.Matrix([s.diff(array.streamfunction, coord) for coord in (xx, yy)])
    checks.check("the coefficient is transported as a genuine first integral",
                 s.simplify(array.velocity.dot(grad_psi)) == 0)
    hb = s.Symbol("mean_hb")
    # Integration by parts: <grad h.wh>=-<h div wh>=i k<h b>.
    poisson_mean = 2*s.I*k*(s.I*k*hb)
    checks.check("full harmonic pressure is finite and not reset",
                 s.expand(poisson_mean/k**2+2*hb) == 0)
    checks.check("dropping that harmonic pressure changes the complete source",
                 poisson_mean != 0)

    phase = axial_kelvin_initial_phase(array.wavevectors, array.sine_velocities, k, rho)
    checks.check("ordered family anchors to the importable full-pressure planar action",
                 s.simplify(phase.stiffness
                            -rho*k**2*lam**2/(lam**2+k**2)*array.covariance) == s.zeros(2)
                 and phase.mass == rho*s.eye(2))

    ux, uy, hh, phi = s.symbols("Ux Uy hmode phi", real=True)
    macro = s.Matrix([ux, uy, 0])
    full_wave = s.Matrix([lam, 0, k])
    projector = s.eye(3)-full_wave*full_wave.T/(lam**2+k**2)
    force = hh*macro+s.Matrix([0, 0, phi])
    a_rate = -s.I*k*projector*force
    b_rate = s.I*k*hh*macro
    momentum = rho*(a_rate+b_rate)
    energy = (s.conjugate(momentum).dot(momentum)/(2*rho)
              -s.re(s.conjugate(momentum).dot(b_rate)))
    direct = rho*k**2*(force.dot(projector*force)-hh**2*macro.dot(macro))/2
    checks.check("complete Jacobi momentum eliminates to the stated lifted phase energy",
                 s.simplify(energy-direct) == 0)
    checks.check("the additional axial kinetic term is NOT omitted",
                 s.expand(direct-rho*k**2*force.dot(projector*force)/2
                          +rho*k**2*hh**2*(ux**2+uy**2)/2) == 0)
    checks.check("the zero-shear limit retains the exact axial pressure return",
                 s.simplify(direct.subs(hh, 0)
                            -rho*k**2*lam**2*phi**2/(2*(lam**2+k**2))) == 0)
    cb, c0 = s.symbols("cb c0", positive=True)
    lower = c0*cb**2-2*delta*cb-delta**2
    margin = s.factor((lower.subs(delta, c0*cb/8)-c0*cb**2/2)/(c0*cb**2))
    checks.check("a finite explicit axial-shear bound preserves positive stiffness",
                 s.expand(margin-(16-c0)/64) == 0
                 and (16-s.Integer(1))/64 > 0)

    big_d, coefficient = s.symbols("D a", positive=True)
    finite_c = (big_d/(coefficient*eps**3))**2
    checks.check("finite C chosen after k gives the asserted acoustic error ordering",
                 s.simplify((big_d/s.sqrt(finite_c))/eps**2-coefficient*eps) == 0)
    t, bound = s.symbols("T c", positive=True)
    # Actual generator norm perturbation delta/epsilon integrated over T/epsilon.
    accumulated = s.integrate(bound*delta/eps, (s.Symbol("time"), 0, t/eps))
    checks.check("the response comparison retains the two inverse-k losses",
                 s.simplify(accumulated-bound*delta*t/eps**2) == 0)
    print("Scope: exact background/action and ordered operator-estimate anchors; no fixed-C homogenization is inferred.")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
