"""Exact gradient-lift orthogonality and complete Schur/normal-form transfer."""

import sympy as s

from substrate_framework.euler_orbit import (
    hermitian_schur_jet,
    micropolar_kinetic_normal_form,
)
from substrate_framework.homogenization import sphere_second_moment
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0065-full-stationary-gradient-dominance")
    count = 11
    size = 2*count+5  # affine, dual responses, Q,S, three new raw cages
    form = s.zeros(size)
    moments = s.Matrix(count, 5, s.symbols("moment0:55"))
    c = s.symbols("c", nonzero=True, real=True)
    for i in range(count):
        form[i, count+i], form[count+i, i] = 1, -1
        for j in range(5):
            form[i, 2*count+j] = moments[i, j]
            form[2*count+j, i] = -moments[i, j]
    form[2*count, 2*count+1], form[2*count+1, 2*count] = c, -c
    change = s.zeros(size, 5)
    for j in range(5):
        change[2*count+j, j] = 1
        for i in range(count):
            change[count+i, j] = -moments[i, j]
    projected = s.simplify(change.T*form*change)
    target = s.zeros(5)
    target[0, 1], target[1, 0] = c, -c
    ledger.check("projected gradient cages preserve the exact physical Q/S KKS block",
                 projected == target)
    ledger.check("each new cage is KKS-orthogonal to every old and new column",
                 projected[:, 2:] == s.zeros(5, 3))
    ledger.check("all eleven affine moments of the new cages vanish",
                 (form*change)[:count, 2:] == s.zeros(count, 3))
    # The future base may retain nonzero rotation moments by adding dual responses.
    retained_moment = s.symbols("retained_rotation", real=True)
    changed_base = change.copy()
    changed_base[count+3, 1] += retained_moment
    ledger.check("gradient orthogonality survives a rotation-retaining base partner",
                 s.simplify(changed_base.T*form*changed_base) == target)

    wave, bond, direction = s.symbols("wave d direction", real=True)
    difference = 2*s.I*s.sin(wave*bond*direction/2)
    ledger.check("geometric centered neighbor tie changes no uniform angle",
                 difference.subs(wave, 0) == 0)
    ledger.check("gradient amplitude has the exact first bond factor",
                 s.diff(difference, wave).subs(wave, 0) == s.I*bond*direction)
    ledger.check("centered bond introduces no spurious second-order amplitude",
                 s.diff(difference, wave, 2).subs(wave, 0) == 0)
    moment2 = sphere_second_moment()
    q = s.Matrix(s.symbols("q0:3", real=True))
    n = s.Matrix(s.symbols("n0:3", real=True))
    averaged = (q.T*moment2*q)[0]*n.dot(n)
    ledger.check("three bond directions and Haar axes control every gradient polarization",
                 s.simplify(averaged-q.dot(q)*n.dot(n)/3) == 0)

    # All coefficients below are a finite noncommuting operator compression;
    # they test algebra, not an independent-cell physical approximation.
    p = (s.Matrix([[3, s.I], [-s.I, 2]]),
         s.Matrix([[1, 2*s.I], [-2*s.I, 0]]), s.Matrix([[2, 1], [1, 3]]))
    nbase = (s.Matrix([[1, s.I], [2, 1]]),
             s.Matrix([[s.I, 1], [0, 2*s.I]]), s.Matrix([[2, 0], [s.I, 1]]))
    hbase = (s.Matrix([[5, 1], [1, 6]]),
             s.Matrix([[0, s.I], [-s.I, 1]]), s.Matrix([[4, 1], [1, 3]]))
    l1, l2 = s.Matrix([[2*s.I, 1], [1, s.I]]), s.Matrix([[1, 2*s.I], [3, 1]])
    carrier = s.symbols("gradient_carrier", positive=True)
    at, al = s.symbols("A_t A_l", positive=True)
    leading = s.diag(at, al)
    hcross = s.Matrix([[1, s.I], [-s.I, 2]])
    nnew = (nbase[0], nbase[1]+l1, nbase[2]+l2)
    hnew = (hbase[0], hbase[1], hbase[2]+carrier*leading+hcross)
    old = hermitian_schur_jet(p, nbase, hbase)
    new = hermitian_schur_jet(p, nnew, hnew)
    ledger.check("complete momentum inverse jets stay exactly unchanged",
                 old.inverse_momentum == new.inverse_momentum)
    ledger.check("gradient-only lift changes no zeroth-order locking block",
                 old.reduced[0] == new.reduced[0])
    r0, r1, _ = old.inverse_momentum
    correction = (
        nbase[0].conjugate().T*r0*l2+l2.conjugate().T*r0*nbase[0]
        +nbase[1].conjugate().T*r0*l1+l1.conjugate().T*r0*nbase[1]
        +l1.conjugate().T*r0*l1
        +nbase[0].conjugate().T*r1*l1+l1.conjugate().T*r1*nbase[0]
    )
    delta2 = s.simplify(new.reduced[2]-old.reduced[2])
    ledger.check("full noncommuting Schur difference retains every mixed jet",
                 s.simplify(delta2-carrier*leading-hcross+correction) == s.zeros(2))
    ledger.check("only the positive retained gradient term grows with carrier",
                 delta2.diff(carrier) == leading)
    ledger.check("first Schur jet remains bounded independently of gradient carrier",
                 new.reduced[1].diff(carrier) == s.zeros(2))

    djet = (s.Matrix([[2, 1], [0, 1]]),
            s.Matrix([[s.I, 0], [1, 2*s.I]]), s.Matrix([[1, 0], [0, 2]]))
    zero = (s.zeros(2), s.zeros(2), s.zeros(2))
    kinetic = hermitian_schur_jet(p, djet, zero)
    ledger.check("full kinetic jets include nonzero D1,D2 and are carrier independent",
                 all(item.diff(carrier) == s.zeros(2) for item in kinetic.reduced))

    rho, inertia, alpha = s.symbols("rho j alpha", positive=True)
    cb = s.symbols("C_base", real=True)
    normal = micropolar_kinetic_normal_form(
        rho, inertia, alpha, cb+carrier*at,
        inertia/4, 0, -inertia/2, wave,
    )
    ledger.check("relative kinetic correction is retained before curvature positivity",
                 s.simplify(normal.transverse_curvature-cb-carrier*at
                            +alpha*inertia/rho) == 0)
    ledger.check("normal-form curvature preserves the positive carrier slope",
                 s.diff(normal.transverse_curvature, carrier) == at)
    ledger.check("normalizing field map is unchanged by the gradient-only cage",
                 normal.field_map.diff(carrier) == s.zeros(2))
    bmix, mu, mp = s.symbols("mixed mass_U mass_Phi", real=True)
    generalized = micropolar_kinetic_normal_form(
        rho, inertia, alpha, cb+carrier*at, mu, mp, bmix, wave,
    )
    ledger.check("the same dominance works with general fixed physical kinetic contrast",
                 s.diff(generalized.transverse_curvature, carrier) == at)

    intensity, area, length, lam, remainder, margin = s.symbols(
        "nu A_star d lambda_abs R_total margin", positive=True
    )
    sufficient = 3*lam*(remainder+margin)/(rho*intensity*length**2*area)
    lower = rho*intensity*length**2*area*sufficient/(3*lam)-remainder
    ledger.check("finite analytic threshold gives a strict common transverse/longitudinal bound",
                 s.simplify(lower-margin) == 0)
    print("Analytic oracle: full-gradient-dominance.md; full stationary reaction space retained.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
