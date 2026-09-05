"""Exact physical observation coordinates and geometric spin connection."""

import sympy as s

from substrate_framework.euler_phase import moving_phase_pullback
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0119-physical-angle-spin")
    t = s.Symbol("t", real=True)
    beta, r = s.symbols("beta r", nonzero=True, real=True)
    a, b = s.Function("a")(t), s.Function("b")(t)
    J = s.Matrix([[0, -1], [1, 0]])
    omega = -beta * J
    h = -beta * r * s.eye(2)
    c = s.Matrix([[a, b]])
    d = c.diff(t) + c * (r * J)
    transform = c.col_join(d)
    wronskian = s.factor(transform.det())
    ledger.check("physical angle/rate determinant includes the moving observation",
                 s.simplify(wronskian - (a*b.diff(t)-b*a.diff(t)-r*(a*a+b*b))) == 0)
    embedding = transform.inv()
    result = moving_phase_pullback(omega, h, embedding, embedding.diff(t))
    mass = beta / wronskian
    ledger.check("physical rate is the first transformed evolution row",
                 s.simplify(result.generator[0, 0]) == 0
                 and s.simplify(result.generator[0, 1]-1) == 0)
    ledger.check("actual symplectic pullback supplies physical coordinate mass",
                 s.simplify(result.symplectic + mass*J) == s.zeros(2))
    ledger.check("scalar equation retains the Wronskian rate",
                 s.simplify(result.generator[1, 1]-wronskian.diff(t)/wronskian) == 0)
    stiffness = s.simplify(result.hamiltonian[0, 0])
    ledger.check("moving action has the required mass connection cross",
                 s.simplify(result.hamiltonian[0, 1]-mass.diff(t)/2) == 0
                 and s.simplify(result.hamiltonian[1, 1]-mass) == 0)
    ledger.check("complete potential yields the scalar physical equation",
                 s.simplify(result.generator[1, 0]+stiffness/mass) == 0)
    spin1, spin2 = s.symbols("spin1 spin2", real=True)
    spin = s.Matrix([[spin1, spin2]])
    actual_rows = s.simplify(spin * embedding)
    ledger.check("physical spin coefficient is a determinant not canonical momentum",
                 s.simplify(actual_rows[1]-c.col_join(spin).det()/wronskian) == 0
                 and s.simplify(actual_rows[0]-spin.col_join(d).det()/wronskian) == 0)
    gamma = s.Symbol("gamma", real=True)
    rotating = s.Matrix([[s.cos(gamma*t), s.sin(gamma*t)]])
    rotating_w = s.trigsimp(rotating.col_join(rotating.diff(t)+rotating*r*J).det())
    ledger.check("a winding change can change physical scalar mass sign",
                 rotating_w == gamma-r)

    # Actual instantaneous material geometry, u=U e3 and k.u=1.
    # Phi=kappa*u.cross(a)/U^2, Igeo=positive tag mass*c_t*U^2/N^2.
    # Covariance Ck=c_t*u is the periodic particle covariance from0114.
    U, kx, ky, kap, lam = s.symbols("U kx ky kap lam", nonzero=True, real=True)
    k = s.Matrix([kx, ky, 1/U])
    u = s.Matrix([0, 0, U])
    ax, ay, az = s.symbols("ax ay az", real=True)
    accel = s.Matrix([ax, ay, az])
    f1, f2, fd1, fd2 = s.symbols("f1 f2 fd1 fd2", real=True)
    phi = s.Matrix([f1, f2, 0])
    # Differentiating u.Phi=0 fixes the longitudinal physical rate.
    phid = s.Matrix([fd1, fd2, -(ax*f1+ay*f2)/U])
    amplitude = (-u.cross(phi)+u*k.dot(u.cross(phi))) / kap
    ledger.check("actual transverse-angle inverse satisfies Kelvin plane",
                 s.simplify(k.dot(amplitude)) == 0
                 and s.simplify(kap*u.cross(amplitude)/U**2-phi) == s.zeros(3, 1))
    # Differentiate Phi = kap*u cross a/U^2. kap is the fixed phase/carrier scale.
    adx, ady = s.symbols("adx ady", real=True)
    # Tangent of k.a=0 uses kdot=-A^T k, so no a-dot tangent is assumed here.
    ad = s.Matrix([adx, ady, s.Symbol("adz", real=True)])
    phid_from_a = kap*(accel.cross(amplitude)+u.cross(ad))/U**2 - 2*az*phi/U
    spin_from_moments = kap*(amplitude.cross(accel)+u.cross(ad))/U**2
    connection = 2*u*accel.dot(phi)/U**2 - 2*k.dot(u.cross(phi))*accel.cross(u)/U**2
    ledger.check("full material spin retains the geometric angle connection",
                 s.simplify(spin_from_moments-phid_from_a-connection) == s.zeros(3, 1))
    ex, ey = s.symbols("ex ey", real=True)
    transverse_mark = s.Matrix([ex, ey, 0])
    axis = transverse_mark-2*U*(kx*ex+ky*ey)*s.Matrix([0, 0, 1])
    spin_over_inertia = phid+connection
    scalar = transverse_mark.dot(phi)
    scalar_rate = ex*fd1+ey*fd2  # Bishop mark derivative is tangent to u.
    scalar_connection = -2*(kx*ax+ky*ay)
    ledger.check("oblique physical axis has positive geometric spin-rate coefficient",
                 s.simplify(axis.dot(spin_over_inertia)-scalar_rate-scalar_connection*scalar) == 0)

    # Derive the two-chart signed determinant directly from the Euler amplitudes.
    A11, A12, A13, A22, A23 = s.symbols("A11 A12 A13 A22 A23", real=True)
    A = s.Matrix([[A11, A12, A13], [A12+lam*U, A22, A23],
                  [A13, A23, -A11-A22]])
    a1, a2 = s.symbols("a1 a2", real=True)
    aa = s.Matrix([a1, a2, -U*(kx*a1+ky*a2)])
    bb = -lam*k.cross(aa)/k.dot(k)
    raw_spin = aa.cross(A*u)+u.cross(A*aa+bb)
    angle_vector = kap*u.cross(aa)/U**2
    axes = (s.Matrix([1, 0, -2*U*kx]), s.Matrix([0, 1, -2*U*ky]))
    determinants = []
    for physical_axis in axes:
        row_c = s.Matrix([[s.diff(physical_axis.dot(angle_vector), var) for var in (a1, a2)]])
        row_s = s.Matrix([[s.diff(physical_axis.dot(raw_spin), var) for var in (a1, a2)]])
        determinants.append(s.factor(row_c.col_join(row_s).det()))
    summed = s.factor(sum(determinants))
    ledger.check("two fixed geometric charts have a sign-definite determinant sum",
                 s.simplify(summed-kap*lam/(U*k.dot(k))) == 0)
    ledger.check("harmonic limit removes the strict two-chart gyro margin",
                 s.simplify(summed.subs(lam, 0)) == 0)

    # Mechanical versus canonical torque: retain both masses and the connection.
    q = s.Function("q")(t)
    M, inertia, K, chi = (s.Function(name)(t) for name in ("M", "I", "K", "chi"))
    mechanical = inertia*q.diff(t)+chi*q
    acceleration = -(M.diff(t)*q.diff(t)+K*q)/M
    torque = s.diff(mechanical, t).subs(q.diff(t, 2), acceleration)
    expected = (inertia.diff(t)-inertia*M.diff(t)/M+chi)*q.diff(t)+(chi.diff(t)-inertia*K/M)*q
    ledger.check("physical torque retains variable inertia and connection current",
                 s.simplify(torque-expected) == 0)
    ledger.check("scalar inertia equality is an additional measured matching condition",
                 s.simplify(mechanical-M*q.diff(t)-(inertia-M)*q.diff(t)-chi*q) == 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
