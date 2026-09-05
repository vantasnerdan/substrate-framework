"""Exact general joint-action jets and physical centroid mixing."""

import sympy as s

from substrate_framework.micropolar import relative_angle_field_map
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0066-full-joint-jet")
    rho, j, kap = s.symbols("rho j kappa", positive=True)
    k, omega2 = s.symbols("k omega2", real=True)
    mu, mp, b, a, g, c = s.symbols("m_U m_Phi b A g C", real=True)
    delta = mp-b*b/rho
    ceff = c-2*g*b/rho-kap*delta/j
    canonical = relative_angle_field_map((0, 0, k))

    def jet(matrix):
        return matrix.applyfunc(lambda entry: sum(s.expand(entry).coeff(k, n)*k**n
                                                  for n in range(3)))

    for h in (-1, 1):
        e = s.Matrix([1, s.I*h, 0])
        embed = s.zeros(6, 2)
        embed[:3, 0], embed[3:, 1] = e, e
        map2 = s.Matrix([[1, 0], [-h*k/2, 1]])
        ledger.check(f"canonical physical-angle map has correct curl helicity {h}",
                     canonical*embed == embed*map2)
        mass = s.Matrix([[rho+mu*k*k, b*h*k], [b*h*k, j+mp*k*k]])
        stiffness = s.Matrix([[a*k*k, g*h*k], [g*h*k, kap+c*k*k]])
        change = s.Matrix([[1-mu*k*k/(2*rho), -b*h*k/rho],
                           [0, 1-delta*k*k/(2*j)]])
        ledger.check(f"general full mass jet normalizes, helicity {h}",
                     s.simplify(jet(change.T*mass*change)-s.diag(rho, j)) == s.zeros(2))
        target = s.Matrix([[a*k*k, g*h*k], [g*h*k, kap+ceff*k*k]])
        ledger.check(f"general mixed stiffness enters complete C_eff, helicity {h}",
                     s.simplify(jet(change.T*stiffness*change)-target) == s.zeros(2))
        determinant = s.expand((stiffness-omega2*mass).det())
        acoustic = (a-g*g/kap)/rho
        optical = c/j-kap*mp/j**2+(g-kap*b/j)**2/(rho*kap)
        ledger.check(f"acoustic squared-frequency slope from full determinant, helicity {h}",
                     s.simplify(s.expand(determinant.subs(omega2, acoustic*k*k)).coeff(k, 2))
                     == 0)
        ledger.check(f"optical squared-frequency slope from full determinant, helicity {h}",
                     s.simplify(s.expand(determinant.subs(omega2, kap/j+optical*k*k))
                                .coeff(k, 2)) == 0)
        amplitude = (j*g/kap-b)*h*k/rho
        row = (stiffness-omega2*mass)[0, :]*s.Matrix([amplitude, 1])
        ledger.check(f"physical centroid optical response uses g-kappa*b/j, helicity {h}",
                     s.simplify(s.expand(row[0].subs(omega2, kap/j)).coeff(k, 1)) == 0)

    m0, d0, n, ar, r, cr, frame = s.symbols("m0 d0 n a r c frame", real=True)
    transform = s.Matrix([[1, 0], [-frame*k, 1]])
    mr = s.Matrix([[rho+m0*k*k, d0*k], [d0*k, j+n*k*k]])
    kr = s.Matrix([[ar*k*k, r*k], [r*k, kap+cr*k*k]])
    mpull, kpull = transform.T*mr*transform, transform.T*kr*transform
    mexpected = s.Matrix([
        [rho+(m0-2*frame*d0+frame**2*j)*k*k+frame**2*n*k**4,
         (d0-frame*j)*k-frame*n*k**3],
        [(d0-frame*j)*k-frame*n*k**3, j+n*k*k]])
    kexpected = s.Matrix([
        [(ar-2*frame*r+frame**2*kap)*k*k+frame**2*cr*k**4,
         (r-frame*kap)*k-frame*cr*k**3],
        [(r-frame*kap)*k-frame*cr*k**3, kap+cr*k*k]])
    ledger.check("exact general relative pullback retains all third/fourth degree terms",
                 s.simplify(mpull-mexpected) == s.zeros(2)
                 and s.simplify(kpull-kexpected) == s.zeros(2))
    ledger.check("physical centroid mixing is independent of the relative angle convention",
                 s.simplify((r-frame*kap)-kap*(d0-frame*j)/j-(r-kap*d0/j)) == 0)
    substitutions = {m0: 0, d0: 0, n: 0, r: 0, frame: s.Rational(1, 2)}
    exact_m, exact_k = mpull.subs(substitutions), kpull.subs(substitutions)
    full_det = s.expand((exact_k-omega2*exact_m).det())
    factor = (ar*k*k-rho*omega2)*(kap+cr*k*k-j*omega2)
    ledger.check("separable physical action factorizes with every pullback term retained",
                 s.simplify(full_det-factor) == 0)
    truncated_det = s.expand((jet(exact_k)-omega2*jet(exact_m)).det())
    difference = s.expand(truncated_det-full_det)
    ledger.check("discarded curvature pullbacks change no determinant coefficient through k²",
                 all(s.simplify(difference.coeff(k, power)) == 0 for power in range(3)))
    ledger.check("discarded pullbacks do change higher-order determinant coefficients",
                 s.simplify(difference) != 0)
    mixing = (g-kap*b/j).subs({g: -kap/2, b: -j/2})
    ledger.check("separable relative-rate sector has zero original-centroid optical mixing",
                 mixing == 0)
    ledger.check("dropping the physical kinetic cross manufactures centroid mixing",
                 (g-kap*b/j).subs({g: -kap/2, b: 0}) != 0)
    ledger.check("computed relative kinetic cross supplies a genuine leading physical response",
                 s.simplify((g-kap*b/j).subs({g: -kap/2, b: d0-j/2})+kap*d0/j) == 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
