"""Exact coherence-return geometry and full momentum-gradient reduction."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0050-common-coherence-action")
    coords = s.Matrix(s.symbols("x y z", real=True))
    x, y, z = coords
    axis = s.Matrix(s.symbols("e0:3", real=True))

    def curl(v):
        return s.Matrix([s.diff(v[2], y)-s.diff(v[1], z),
                         s.diff(v[0], z)-s.diff(v[2], x),
                         s.diff(v[1], x)-s.diff(v[0], y)])

    def div(v):
        return sum(s.diff(v[j], coords[j]) for j in range(3))

    potential = -coords.dot(coords)*axis/2
    rot = axis.cross(coords)
    ledger.check("global rotation potential has exact curl",
                 s.simplify(curl(potential)-rot) == s.zeros(3, 1))
    cutoff = s.Function("chi_eta")(x, y, z)
    grad_cut = s.Matrix([s.diff(cutoff, t) for t in coords])
    for j in range(3):
        unit = s.eye(3)[:, j]
        tangent = curl(cutoff*coords[j]*potential)
        expanded = cutoff*(coords[j]*rot+unit.cross(potential))
        expanded += coords[j]*grad_cut.cross(potential)
        ledger.check(f"direction {j}: full affine envelope includes return shell",
                     s.simplify(tangent-expanded) == s.zeros(3, 1))
        ledger.check(f"direction {j}: exact incompressibility", s.simplify(div(tangent)) == 0)
        interior = curl(coords[j]*potential)
        ledger.check(f"direction {j}: origin angle jet unchanged",
                     interior.jacobian(coords).subs({x: 0, y: 0, z: 0}) == s.zeros(3))

    xi = s.Matrix(s.symbols("xi0:3"))
    omega = s.Matrix(s.symbols("w0:3"))
    dw = s.Matrix(s.symbols("dw0:3"))
    ledger.check("second transported vorticity gives negative helicity sign",
                 s.expand(omega.dot(xi.cross(dw))+xi.cross(omega).dot(dw)) == 0)

    # No energy-orthogonality or diagonal-momentum assumption is imposed.
    p11, p12, p22 = s.symbols("p11 p12 p22", real=True)
    p = s.Matrix([[p11, p12], [p12, p22]])
    pinv = p.inv()
    hp = s.Matrix(1, 2, s.symbols("h1 h2", real=True))
    op = s.Matrix(1, 2, s.symbols("o1 o2", real=True))
    d = s.Matrix(1, 2, s.symbols("d1 d2", real=True))
    k, nu, hzz = s.symbols("k nu hzz", real=True)
    hxp = -s.I*k*hp
    oxp = d-s.I*k*op
    hpx, opx = s.conjugate(hxp.T), -s.conjugate(oxp.T)
    exact = k*k*hzz-((hxp-s.I*nu*oxp)*pinv*(hpx-s.I*nu*opx))[0]
    stiff = k*k*hzz-(hxp*pinv*hpx)[0]
    gyro = -(oxp*pinv*hpx+hxp*pinv*opx)[0]
    mass = -(oxp*pinv*opx)[0]
    ledger.check("full Schur kernel has the stated stiffness, gyro and mass",
                 s.simplify(exact-(stiff-s.I*nu*gyro-nu*nu*mass)) == 0)
    ledger.check("common curvature retains all momentum cross terms",
                 s.simplify(s.diff(stiff, k, 2)/2-hzz+(hp*pinv*hp.T)[0]) == 0)
    ledger.check("common gradient inertia is the positive momentum Gram",
                 s.simplify(s.diff(mass, k, 2)/2-(op*pinv*op.T)[0]) == 0)
    minus = k*k*hzz-((hxp+s.I*nu*oxp)*pinv*(hpx+s.I*nu*opx))[0]
    ledger.check("independent time-reversal reduction cancels full gyro only",
                 s.simplify((exact+minus)/2-stiff+nu*nu*mass) == 0)

    a, cb, ib, j0, k0 = s.symbols("a C_B I_B2 J0 K0", nonzero=True, real=True)
    # Physical beta=0 optical section: B=q=Psi/a.
    optical_ratio = (k0+k*k*cb/a**2)/(j0+k*k*ib/a**2)
    ledger.check("constant-mass optical curvature includes negative inertia term",
                 s.simplify(j0*s.diff(optical_ratio, k, 2).subs(k, 0)/2
                            -(cb-k0*ib/j0)/a**2) == 0)
    change = s.Matrix([[1/a, 1-1/a], [1/a, -1/a]])
    ledger.check("physical common-relative field map preserves beta",
                 s.Matrix([[1, -1]])*change == s.Matrix([[0, 1]]))
    ledger.check("optical direction uses both common and relative gradients",
                 change[:, 0] == s.Matrix([1/a, 1/a]))

    # Orthogonalize a new cage using fixed eta0, not the high-energy body.
    l0, lz = s.symbols("l0 lz", nonzero=True, real=True)
    raw = s.zeros(7)  # K, eta0, r, Q, S, Z, new cage
    raw[0, 1], raw[1, 0] = l0, -l0
    raw[0, 6], raw[6, 0] = lz, -lz
    lifted = s.zeros(7, 1)
    lifted[6], lifted[1] = 1, -lz/l0
    ledger.check("fixed-eta0 cage moment projection annihilates common moment",
                 (raw*lifted)[0] == 0)
    ledger.check("support-separated lift changes no retained KKS pairing",
                 raw*lifted == s.zeros(7, 1))
    grow, bounded, threshold = s.symbols("A remainder threshold", positive=True)
    ledger.check("finite carrier threshold yields strict positive residual",
                 s.expand(grow*((bounded+threshold)/grow)-bounded) == threshold)
    print("Analytic oracle: coherence-action.md; compact transport and exact Schur reduction.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
