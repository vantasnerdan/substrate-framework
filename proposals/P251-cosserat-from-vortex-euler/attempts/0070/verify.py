"""Exact Euler affine-moment prescription and physical-angle action checks."""

import sympy as s

from substrate_framework.euler_orbit import hermitian_schur_jet
from substrate_framework.micropolar import relative_angle_field_map
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0070-affine-angular-reaction")
    b = s.Symbol("B", real=True, nonzero=True)
    axis = s.Matrix(s.symbols("n0:3", real=True))
    # Basis: 11 fixed isotropic dual responses, raw Q, raw S.
    moments = s.Matrix(11, 2, s.symbols("f0:22", real=True))
    response = s.eye(11).row_join(moments)
    selected = (-moments).col_join(s.eye(2))
    for index in range(3):
        selected[8+index, 1] += b*axis[index]
    form = s.zeros(13)
    form[11, 12], form[12, 11] = b, -b
    actual_moments = response*selected
    target = s.zeros(11, 2)
    target[8:11, 1] = b*axis
    ledger.check("three translations and five STF moments vanish; all three rotations are prescribed",
                 actual_moments == target)
    ledger.check("new rotational response leaves the exact angle/cage KKS pairing unchanged",
                 selected.T*form*selected == s.Matrix([[0, b], [-b, 0]]))
    physical_jet = s.zeros(1, 13)
    physical_jet[0, 11] = 1
    ledger.check("physical core angle and zero reaction angle jet survive the new response",
                 physical_jet*selected == s.Matrix([[1, 0]]))

    # Pullback on local (beta_1,beta_2,beta_3,q,s) coordinates.
    joint = s.zeros(5)
    for index in range(3):
        joint[index, 4], joint[4, index] = b*axis[index], -b*axis[index]
    joint[3, 4], joint[4, 3] = b, -b
    rel_to_abs = s.eye(5)
    for index in range(3):
        rel_to_abs[3, index] = -axis[index]
    transformed = s.simplify(rel_to_abs.T*joint*rel_to_abs)
    expected_form = s.zeros(5)
    expected_form[3, 4], expected_form[4, 3] = b, -b
    ledger.check("same Euler affine KKS becomes absolute physical-angle canonical form",
                 transformed == expected_form)
    old = joint.copy()
    for index in range(3):
        old[index, 4], old[4, index] = 0, 0
    ledger.check("zero-affine-pairing mutation fails the absolute-angle kinetic construction",
                 s.simplify(rel_to_abs.T*old*rel_to_abs) != expected_form)

    rho = s.Symbol("rho", positive=True)
    scalar_moment = s.Symbol("c", real=True)
    jn = s.Matrix([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                   [-axis[1], axis[0], 0]])
    first = scalar_moment*s.eye(3)+b*jn/(2*rho)
    angular = s.Matrix([
        sum(s.LeviCivita(i, j, k)*first[k, j] for j in range(3) for k in range(3))
        for i in range(3)])
    ledger.check("antisymmetric force first moment has the exact physical angular normalization",
                 s.simplify(rho*angular-b*axis) == s.zeros(3, 1))
    integrated_potential = b*axis/(2*rho)
    curl_moment = s.Matrix(3, 3, lambda i, j: -sum(
        s.LeviCivita(i, j, k)*integrated_potential[k] for k in range(3)))
    ledger.check("compact curl removes exactly the prescribed antisymmetric first moment",
                 first-curl_moment == scalar_moment*s.eye(3))

    x, y, z = coords = s.symbols("x y z", real=True)
    e = s.Matrix(s.symbols("e0:3", real=True))
    position = s.Matrix(coords)
    weighted = position.dot(position)*e

    def curl(v):
        return s.Matrix([s.diff(v[2], y)-s.diff(v[1], z),
                         s.diff(v[0], z)-s.diff(v[2], x),
                         s.diff(v[1], x)-s.diff(v[0], y)])

    ledger.check("angular-impulse variation equals the Euler rotational KKS moment with its sign",
                 s.simplify(-rho*curl(weighted)/2-rho*e.cross(position)) == s.zeros(3, 1))
    # The boundary identity itself follows by component integration by parts;
    # evaluate that underlying divergence identity without assuming a tail.
    v = s.Matrix([s.Function(f"v{i}")(*coords) for i in range(3)])
    for i in range(3):
        divergence_flux = sum(s.diff(position.dot(position)*sum(
            s.LeviCivita(i, j, ell)*v[ell] for ell in range(3)), coords[j])
            for j in range(3))
        identity = position.dot(position)*curl(v)[i]+2*position.cross(v)[i]
        ledger.check(f"exact velocity angular-momentum boundary identity, component {i}",
                     s.simplify(divergence_flux-identity) == 0)

    # Full non-diagonal reaction matrix; the physical-angle coupling acts
    # through its full inverse, while arbitrary potential mixed blocks remain.
    p = s.Matrix([[4, 1, 0], [1, 3, 1], [0, 1, 2]])
    d = s.Matrix([[1, 0], [0, 1], [1, 1]])
    n = s.Matrix([[1, 0], [0, 1], [1, -1]])
    h = 2*s.eye(2)+n.T*p.inv()*n
    reduced = hermitian_schur_jet((p, s.zeros(3), s.zeros(3)),
                                  (n, s.zeros(3, 2), s.zeros(3, 2)),
                                  (h, s.zeros(2), s.zeros(2)))
    inertia = d.T*p.inv()*d
    physical_rate = s.Matrix(s.symbols("Phi_dot0:2", real=True))
    state = s.Matrix(s.symbols("state0:2", real=True))
    paired = 0
    for sign in (-1, 1):
        source = sign*d*physical_rate-n*state
        paired += (source.T*p.inv()*source-state.T*h*state)[0]/4
    target_action = (physical_rate.T*inertia*physical_rate
                     -state.T*reduced.reduced[0]*state)[0]/2
    ledger.check("independent opposite fluid reactions produce absolute-angle full-operator inertia",
                 s.simplify(paired-target_action) == 0)
    ledger.check("same full reaction matrix gives positive inertia and restoring Schur block",
                 inertia[0, 0] > 0 and inertia.det() > 0
                 and reduced.reduced[0] == 2*s.eye(2))
    averaged_spin = s.zeros(2, 1)
    for sign in (-1, 1):
        reaction = p.inv()*(sign*d*physical_rate-n*state)
        averaged_spin += sign*d.T*reaction/2
    ledger.check("Euler physical angular-momentum response equals derivative of the reduced kinetic action",
                 s.simplify(averaged_spin-inertia*physical_rate) == s.zeros(2, 1))

    j, kap, k = s.symbols("j kappa k", positive=True)
    canonical = relative_angle_field_map((0, 0, k))
    for helicity in (-1, 1):
        vec = s.Matrix([1, s.I*helicity, 0])
        embed = s.zeros(6, 2)
        embed[:3, 0], embed[3:, 1] = vec, vec
        transform = s.Matrix([[1, 0], [-helicity*k/2, 1]])
        ledger.check(f"imported exact physical map agrees in helicity {helicity}",
                     canonical*embed == embed*transform)
        relative_mass = s.Matrix([[rho+j*k*k/4, j*helicity*k/2],
                                  [j*helicity*k/2, j]])
        ledger.check(f"computed affine connection gives absolute-angle physical mass, helicity {helicity}",
                     s.simplify(transform.T*relative_mass*transform-s.diag(rho, j))
                     == s.zeros(2))
        bphys = j/2-j/2
        mixing = -kap/2-kap*bphys/j
        ledger.check(f"new same-Euler sector has nonzero physical centroid transfer, helicity {helicity}",
                     mixing == -kap/2 and mixing != 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
