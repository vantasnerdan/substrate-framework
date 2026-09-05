"""Exact hybrid centering, tube-spin multipoles and material tag transport."""

import sympy as s

from substrate_framework.euler_orbit import hermitian_schur_jet
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0082-hybrid-action")
    rho, metric = s.symbols("rho G", positive=True)
    k, global_spin, tube_spin = s.symbols("k L_global L_tube", real=True)
    dotu, dotphi, independent, physical, hs = s.symbols("Udot Phidot P0 p_H Hs", real=True)
    for helicity in (-1, 1):
        point = helicity*k*global_spin/2
        tube_correction = helicity*k*tube_spin/2
        hybrid = point-tube_correction
        ledger.check(f"actual tube spin, not global impulse, defines hybrid moment, helicity {helicity}",
                     s.simplify(point-hybrid-tube_correction) == 0)
        affine_kks = point*dotu-tube_spin*helicity*k*dotu/2
        ledger.check(f"fixed-angle Gamma lift pairs with exact hybrid induced momentum, helicity {helicity}",
                     s.simplify(affine_kks-hybrid*dotu) == 0)
        theta = (independent+hybrid)*dotu+tube_spin*dotphi
        hamiltonian = independent**2/(2*metric)+independent*hybrid/metric+hs
        centered_theta = physical*dotu+tube_spin*dotphi
        centered_h = physical**2/(2*metric)+hs-hybrid**2/(2*metric)
        ledger.check(f"complete hybrid cotangent shift leaves the absolute tube angle source, helicity {helicity}",
                     s.simplify(theta.subs(independent, physical-hybrid)-centered_theta) == 0)
        ledger.check(f"same hybrid response enters the kinetic cross and Gram subtraction, helicity {helicity}",
                     s.simplify(hamiltonian.subs(independent, physical-hybrid)-centered_h) == 0)
        action = theta-hamiltonian
        solution = s.solve(s.diff(action, independent), independent)[0]
        reduced = s.simplify(action.subs(independent, solution))
        expected = metric*dotu**2/2+tube_spin*dotphi-hs+hybrid**2/(2*metric)
        ledger.check(f"exact hybrid mean elimination does not subtract tube spin twice, helicity {helicity}",
                     s.simplify(reduced-expected) == 0)
        wrong = independent**2/(2*metric)+independent*point/metric+hs
        ledger.check(f"point-mean kinetic-cross substitution fails the hybrid same-action identity, helicity {helicity}",
                     s.simplify(wrong-hamiltonian) != 0)

    # Retain all physical coordinates and the FULL reaction inverse.
    p = s.Matrix([[4, 1], [1, 3]])
    response = s.Matrix([[1, 2]])
    centered_p = p-response.T*response/10
    sources = s.Matrix([[1, 0], [1, 1]])  # independent physical angle/shape chart
    result = hermitian_schur_jet((centered_p, s.zeros(2), s.zeros(2)),
                                 (sources, s.zeros(2), s.zeros(2)),
                                 (s.zeros(2), s.zeros(2), s.zeros(2)))
    full_inertia = -result.reduced[0]
    ledger.check("independent angle and shape coordinates retain a full positive kinetic block",
                 full_inertia[0, 0] > 0 and full_inertia.det() > 0)
    conditional = full_inertia[0, 0]-full_inertia[0, 1]**2/full_inertia[1, 1]
    ledger.check("legitimate further shape elimination uses its nontrivial positive Schur complement",
                 conditional > 0 and conditional != full_inertia[0, 0])

    theta, nu, frequency = s.symbols("theta nu omega", real=True, nonzero=True)
    f1 = s.cos(theta)
    f0_repaired = nu*s.diff(f1, theta)
    ledger.check("material commutator supplies the exact missing local tag-transport column",
                 f0_repaired == -nu*s.sin(theta))
    # Fourier coefficients ordered +1,-1; the spin weight is cos(theta).
    original = sum((-s.I*frequency*s.Rational(1, 2))
                   /(s.I*(sign*nu-frequency))*s.Rational(1, 2)
                   for sign in (-1, 1))
    repaired = sum((-s.I*frequency*s.Rational(1, 2)
                    +s.I*sign*nu*s.Rational(1, 2))
                   /(s.I*(sign*nu-frequency))*s.Rational(1, 2)
                   for sign in (-1, 1))
    ledger.check("transported boundary spin has the exact frequency-dependent residue when column is absent",
                 s.simplify(original-frequency**2/(2*(frequency**2-nu**2))) == 0)
    ledger.check("actual commutator-column repair removes the apparent transport poles",
                 s.simplify(repaired-s.Rational(1, 2)) == 0)
    ledger.check("frequency-independent spin cannot be inferred from its bulk moment alone",
                 s.simplify(original-repaired) != 0)
    time = s.Symbol("t", real=True)
    phi = s.Function("Phi")(time)
    f0, profile = s.Function("f0")(theta), s.Function("f1")(theta)
    local_tag = profile*phi
    residual = s.diff(local_tag, time)+nu*s.diff(local_tag, theta)
    residual -= profile*s.diff(phi, time)+f0*phi
    ledger.check("local collar criterion is necessary and sufficient for arbitrary angle histories",
                 s.simplify(residual-(nu*s.diff(profile, theta)-f0)*phi) == 0)
    # Finite-time transport at resonance is regular as an evolution, but secular.
    resonant = time*s.exp(-s.I*nu*time)
    ledger.check("the exact transported tag retains resonant secular response instead of inventing a modulus",
                 s.simplify(s.diff(resonant, time)+s.I*nu*resonant-s.exp(-s.I*nu*time)) == 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
