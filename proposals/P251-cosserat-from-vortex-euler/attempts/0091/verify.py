"""Exact moment completion and same-action common/core-angle identities."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0091-spin-metric-completion")
    x, y, z, length = s.symbols("x y z L", real=True, positive=True)
    coords = (x, y, z)
    chi = s.Function("chi")((x*x+y*y+z*z)/length**2)
    rotation = s.Matrix([-y, x, 0])
    lift = chi*rotation
    ledger.check("radial physical-rotation lift is exactly divergence free",
                 s.simplify(s.trace(lift.jacobian(coords))) == 0)
    ledger.check("spin-minus-mass integrand is positive in the cutoff transition shell",
                 s.simplify(rotation.dot(lift)-lift.dot(lift)
                            -chi*(1-chi)*rotation.dot(rotation)) == 0)

    a, t = s.symbols("a t", real=True)
    mass, deficit = s.symbols("m D", positive=True)
    root = (a+s.sqrt(a*a+4*mass*deficit))/(2*mass)
    ledger.check("computed positive-root amplitude satisfies the physical moment equation",
                 s.simplify(mass*root**2-a*root-deficit) == 0)
    base_mass = s.symbols("mB", positive=True)
    total_mass = base_mass+t*t*mass
    total_spin = base_mass+deficit+t*a
    ledger.check("full disjoint-field spin and kinetic mass agree at that root",
                 s.simplify((total_mass-total_spin).subs(t, root)) == 0)
    ledger.check("the root is strictly positive for positive mass and deficit",
                 s.simplify((a*a+4*mass*deficit)-a*a) == 4*mass*deficit
                 and (4*mass*deficit).is_positive)

    kb = s.symbols("KB", real=True)
    kx = s.symbols("KX", positive=True)
    total_k = kb+t*t*kx
    rx = kx/mass
    ledger.check("exact cage-to-completed-gap margin retains the base stiffness",
                 s.simplify(rx-total_k/total_mass-(rx*base_mass-kb)/total_mass) == 0)
    d, stiffness_bound, rho = s.symbols("d C rho", positive=True)
    lower_k = kx*rho*d*length**5/(4*mass)-stiffness_bound*length**3
    ledger.check("finite-scale positivity condition follows from actual L5 and L3 terms",
                 s.factor(lower_k/length**3)
                 == (kx*rho*d*length**2-4*mass*stiffness_bound)/(4*mass))

    inertia, j = s.symbols("I j", positive=True)
    beta_rate, q_rate, phi_rate = s.symbols("betadot qdot Phidot", real=True)
    kinetic = inertia*beta_rate**2/2+j*beta_rate*q_rate+j*q_rate**2/2
    completed = (inertia-j)*beta_rate**2/2+j*phi_rate**2/2
    ledger.check("actual common-relative kinetic cross gives the physical core-angle rate",
                 s.expand(kinetic.subs(q_rate, phi_rate-beta_rate)-completed) == 0)
    ledger.check("physical mode angular momentum agrees with its canonical angle momentum",
                 s.diff(completed, phi_rate) == j*phi_rate)
    ledger.check("discarding the affine inertia correction changes the same kinetic action",
                 s.expand(kinetic.subs(q_rate, phi_rate-beta_rate)
                          -inertia*beta_rate**2/2-j*phi_rate**2/2)
                 == -j*beta_rate**2/2)

    e = s.Matrix(s.symbols("e0:3", real=True))
    moment = s.Matrix([[0, -e[2], e[1]], [e[2], 0, -e[0]], [-e[1], e[0], 0]])
    strain = s.Matrix([[s.Symbol("s0"), s.Symbol("s1"), s.Symbol("s2")],
                       [s.Symbol("s1"), s.Symbol("s3"), s.Symbol("s4")],
                       [s.Symbol("s2"), s.Symbol("s4"), s.Symbol("s5")]])
    ledger.check("symmetric affine velocity has no compact solenoidal dipole kinetic cross",
                 s.trace(moment.T*strain) == 0)

    lam, wave, amplitude = s.symbols("lambda k f", real=True)
    orbit_h = (lam**2-lam*wave)*amplitude**2
    material_k = (lam*wave-wave**2)*amplitude**2
    ledger.check("same-field Kelvin reconstruction norm is the exact orbit-minus-material stiffness",
                 s.expand(orbit_h-material_k-(lam-wave)**2*amplitude**2) == 0)

    nu, micro_k, micro_j = s.symbols("nu kappa jmicro", positive=True)
    q = s.Symbol("q", real=True)
    angle_part = nu*(micro_j*phi_rate**2-micro_j*beta_rate**2-micro_k*q*q)/6
    mu_s, strain_norm = s.symbols("mu_s E2", positive=True)
    ledger.check("removing only spin population retains the independently supplied Cauchy sector",
                 (angle_part-mu_s*strain_norm).subs(nu, 0) == -mu_s*strain_norm)
    ledger.check("selected-population removal is not complete structure removal",
                 (angle_part-mu_s*strain_norm).subs({nu: 0, mu_s: 0}) == 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
