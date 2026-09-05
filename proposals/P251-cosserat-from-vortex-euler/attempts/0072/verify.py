"""Exact induced-mean, cotangent and reaction-centering repair."""

import sympy as s

from substrate_framework.euler_orbit import hermitian_schur_jet
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0072-mean-cotangent-repair")
    rho, j, kap = s.symbols("rho j kappa", positive=True)
    k, scalar = s.symbols("k c", real=True)
    wave = s.Matrix(s.symbols("k0:3", real=True))
    spin = s.Matrix(s.symbols("L0:3", real=True))
    cross_spin = s.Matrix([[0, -spin[2], spin[1]],
                           [spin[2], 0, -spin[0]], [-spin[1], spin[0], 0]])
    first = scalar*s.eye(3)+cross_spin/(2*rho)
    coherent_force = -s.I*first*wave
    expected = -s.I*scalar*wave+s.I*wave.cross(spin)/(2*rho)
    ledger.check("modulated antisymmetric force moment produces curl spin with its exact sign",
                 s.simplify(coherent_force-expected) == s.zeros(3, 1))
    projector = s.eye(3)-wave*wave.T/wave.dot(wave)
    ledger.check("full Leray mean kills only the scalar moment, not the angular moment",
                 s.simplify(projector*coherent_force-s.I*wave.cross(spin)/(2*rho))
                 == s.zeros(3, 1))

    p0, pdot, ell, qdot, u_dot, hs = s.symbols("P0 p L qdot Udot Hs", real=True)
    for h in (-1, 1):
        induced = h*k*ell/(2*rho)
        phi_dot = qdot+h*k*u_dot/2
        theta_old = p0*u_dot+ell*phi_dot
        theta_center = pdot*u_dot+ell*qdot
        ledger.check(f"complete cotangent shift removes the same induced-mean source, helicity {h}",
                     s.simplify(theta_old.subs(p0, pdot-rho*induced)-theta_center) == 0)
        hamiltonian = p0**2/(2*rho)+p0*induced+hs
        centered = pdot**2/(2*rho)+hs-rho*induced**2/2
        ledger.check(f"same Euler kinetic square includes the mandatory mean cross, helicity {h}",
                     s.simplify(hamiltonian.subs(p0, pdot-rho*induced)-centered) == 0)
        lagrangian = theta_old-hamiltonian
        p_solution = s.solve(s.diff(lagrangian, p0), p0)[0]
        eliminated = s.simplify(lagrangian.subs(p0, p_solution))
        target = rho*u_dot**2/2+ell*qdot-hs+rho*induced**2/2
        ledger.check(f"full mean elimination restores physical relative-angle inertia, helicity {h}",
                     s.simplify(eliminated-target) == 0)
        physical_momentum = s.simplify(p_solution+rho*induced)
        ledger.check(f"U remains the physical mean coordinate after centering, helicity {h}",
                     physical_momentum == rho*u_dot)
        wrong = theta_old-p0**2/(2*rho)-hs
        wrong_reduced = wrong.subs(p0, rho*u_dot)
        ledger.check(f"omitting the cross manufactures an extra affine-rate coupling, helicity {h}",
                     s.simplify(s.diff(wrong_reduced-target, ell, u_dot)) == h*k/2)

    # A full finite reaction-matrix jet demonstrates centering BEFORE inversion.
    p = s.Matrix([[4, 1], [1, 3]])
    d = s.Matrix([[1], [2]])
    n = s.Matrix([[1], [-1]])
    hqq = s.Matrix([[2]])+n.T*p.inv()*n
    # Transverse scalar mean: C(k)=k D*/(2rho).
    mean2 = -d*d.T/(4*rho)
    data = hermitian_schur_jet((p, s.zeros(2), mean2),
                               (n, s.zeros(2, 1), s.zeros(2, 1)),
                               (hqq, s.zeros(1), s.zeros(1)))
    j0 = (d.T*p.inv()*d)[0]
    j2 = (d.T*data.inverse_momentum[2]*d)[0]
    ledger.check("mean subtraction changes the full reaction inverse at retained second degree",
                 s.simplify(j2-j0**2/(4*rho)) == 0)
    ledger.check("mixed energy reaction also changes the restoring second-gradient Schur block",
                 s.simplify(data.reduced[2][0]
                            +(n.T*p.inv()*d)[0]**2/(4*rho)) == 0)
    ledger.check("using the uncentered reaction inverse misses a nonzero retained coefficient",
                 j2 != 0)
    b, mphi = -j/2, j*j/(4*rho)
    ledger.check("centered illustrative relative-rate inertia has vanishing mass-normal-form delta",
                 s.simplify(mphi-b*b/rho) == 0)
    cap_c = s.Symbol("C", real=True)
    ceff = cap_c-2*(-kap/2)*b/rho-kap*(mphi-b*b/rho)/j
    ledger.check("complete centered example has C_eff=C-kappa*j/(2rho)",
                 s.simplify(ceff-(cap_c-kap*j/(2*rho))) == 0)
    ledger.check("leading point-mean centroid transfer cancels after the repair",
                 s.simplify(-kap/2-kap*b/j) == 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
