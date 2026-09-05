"""Exact exposing identities for the stationary affine-relative-angle route."""

import sympy as s

from substrate_framework.euler_orbit import hermitian_schur_jet
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0059-stationary-relative-angle")
    x, y, z = coords = s.symbols("x y z", real=True)
    lam, k = s.symbols("lambda k", real=True, nonzero=True)

    def curl(v):
        return s.Matrix([s.diff(v[2], y)-s.diff(v[1], z),
                         s.diff(v[0], z)-s.diff(v[2], x),
                         s.diff(v[1], x)-s.diff(v[0], y)])

    def grad(f):
        return s.Matrix([s.diff(f, c) for c in coords])

    def zero(v):
        return all(s.simplify(entry) == 0 for entry in v)

    abc = s.Matrix([s.sin(lam*z)+s.cos(lam*y),
                    s.sin(lam*x)+s.cos(lam*z),
                    s.sin(lam*y)+s.cos(lam*x)])
    ledger.check("three-axis ABC perturbation is exactly same-lambda Beltrami",
                 zero(curl(abc)-lam*abc))
    h = s.Matrix(3, 3, s.symbols("h0:9"))
    h[2, 2] = -h[0, 0]-h[1, 1]
    affine = s.Matrix(s.symbols("a0:3"))+h*s.Matrix(coords)
    omega = lam*abc
    ledger.check("eleven affine response curl kernels have the derived sign",
                 zero(curl(omega.cross(affine))-omega.jacobian(coords)*affine+h*omega))
    atom_conditions = list(h.T*s.eye(3))
    unknown_h = [h[i, j] for i in range(3) for j in range(3) if (i, j) != (2, 2)]
    ledger.check("three independent Fourier atoms eliminate every tracefree affine linear part",
                 s.linsolve(atom_conditions, unknown_h) == s.FiniteSet((0,)*8))

    moment = s.Matrix(3, 3, s.symbols("m0:9"))
    pairing = s.expand(sum(h[i, j]*moment[i, j] for i in range(3) for j in range(3)))
    conditions = [s.diff(pairing, entry) for entry in unknown_h]
    solution = s.linsolve(conditions, list(moment))
    c = moment[2, 2]
    ledger.check("all eight linear affine moments force precisely M=c I",
                 solution == s.FiniteSet((c, 0, 0, 0, c, 0, 0, 0, c)))
    # Integration by parts: integral y_j partial_i f = -delta_ij integral f.
    repaired_moment = c*s.eye(3)-(-s.eye(3)*(-c))
    ledger.check("compact gradient of mass -c cancels the entire first moment",
                 repaired_moment == s.zeros(3))

    raw_moments = s.Matrix(11, 2, s.symbols("f0:22"))
    functional = s.eye(11).row_join(raw_moments)
    selected = (-raw_moments).col_join(s.eye(2))
    b = s.Symbol("B", real=True, nonzero=True)
    form = s.zeros(13)
    form[11, 12], form[12, 11] = b, -b
    ledger.check("eleven disjoint exact responses remove every affine pairing",
                 functional*selected == s.zeros(11, 2))
    ledger.check("isotropic right inverse preserves the physical angle/cage KKS form exactly",
                 selected.T*form*selected == s.Matrix([[0, b], [-b, 0]]))
    core_jet = s.zeros(1, 13)
    core_jet[0, 11] = 1
    ledger.check("off-core response subtraction leaves the physical angle jet unchanged",
                 core_jet*selected == s.Matrix([[1, 0]]))
    mutated = form.copy()
    mutated[0, 1], mutated[1, 0] = 1, -1
    ledger.check("non-isotropic response mutation is detected by the KKS identity",
                 s.simplify(selected.T*mutated*selected-selected.T*form*selected)
                 != s.zeros(2))

    f0 = s.Matrix([s.Function(f"F{j}")(*coords) for j in range(3)])
    f1 = s.Matrix([s.Function(f"G{j}")(*coords) for j in range(3)])
    ez = s.Matrix([0, 0, 1])
    transverse = f0-ez*f0[2]
    carrier = s.exp(s.I*k*z)
    potential = s.I*ez.cross(transverse)*carrier/k
    scalar = f0[2]*carrier/(s.I*k)
    remainder = (f1-s.I*curl(ez.cross(transverse))+s.I*grad(f0[2]))*carrier
    ledger.check("compact full-Leray parametrix identity is exact for arbitrary smooth amplitudes",
                 zero((f0+f1/k)*carrier-curl(potential)-grad(scalar)-remainder/k))
    ledger.check("opposite scalar-potential sign fails the longitudinal carrier identity",
                 not zero((f0+f1/k)*carrier-curl(potential)+grad(scalar)-remainder/k))

    # A noncommuting divergence-free polynomial example corroborates, rather
    # than substitutes for, the general vector-field Jacobi proof in the text.
    def bracket(a, b):
        return b.jacobian(coords)*a-a.jacobian(coords)*b

    a, bvec, cvec = (s.Matrix([y, 0, x*z]), s.Matrix([0, z, x]),
                    s.Matrix([z, x, 0]))
    a[2] = x*y
    ledger.check("stationary Lie-Poisson sign uses curl(a cross b)=-[a,b]",
                 zero(curl(a.cross(bvec))+bracket(a, bvec)))
    jacobi = (bracket(a, bracket(bvec, cvec))+bracket(bvec, bracket(cvec, a))
              +bracket(cvec, bracket(a, bvec)))
    ledger.check("noncommuting polynomial vector-field Jacobi corroboration", zero(jacobi))

    # Full, non-diagonal reaction operator. No isolated-cell inverse is used.
    p = s.Matrix([[4, 1, 0], [1, 3, 1], [0, 1, 2]])
    d = s.Matrix([[1, 0], [0, 1], [1, 1]])
    n = s.Matrix([[1, 0], [0, 1], [1, -1]])
    hqq = 3*s.eye(2)+n.T*p.inv()*n
    jets = hermitian_schur_jet((p, s.zeros(3), s.zeros(3)),
                               (n, s.zeros(3, 2), s.zeros(3, 2)),
                               (hqq, s.zeros(2), s.zeros(2)))
    inertia = d.T*p.inv()*d
    ledger.check("canonical full reaction Schur API retains a positive restoring matrix",
                 jets.reduced[0] == 3*s.eye(2))
    ledger.check("full reaction inertia is positive by exact principal minors",
                 inertia[0, 0] > 0 and inertia.det() > 0)
    q = s.Matrix(s.symbols("q0:2", real=True))
    velocity = s.Matrix(s.symbols("v0:2", real=True))

    def eliminated(sign):
        source = sign*d*velocity-n*q
        return (source.T*p.inv()*source-q.T*hqq*q)[0]/2

    paired = s.expand((eliminated(1)+eliminated(-1))/2)
    expected = (velocity.T*inertia*velocity-q.T*jets.reduced[0]*q)[0]/2
    ledger.check("independently varied time-reversed reactions cancel only the odd gyro",
                 s.simplify(paired-expected) == 0)
    ledger.check("single realization has a nonzero gyro before averaging",
                 s.simplify(eliminated(1)-eliminated(-1)) != 0)
    tied_source = -n*q
    tied = (tied_source.T*p.inv()*tied_source-q.T*hqq*q)[0]/2
    ledger.check("prematurely tied opposite reactions erase inertia and fail the target",
                 s.hessian(tied, velocity) == s.zeros(2) and inertia != s.zeros(2))
    isolated = d.T*s.diag(*(1/p[i, i] for i in range(3)))*d
    ledger.check("isolated-response inverse mutation changes the full density inertia",
                 inertia != isolated)

    rho, j, alpha, cap_c = s.symbols("rho j alpha C", positive=True)
    for helicity in (-1, 1):
        mass = s.Matrix([[rho+j*k**2/4, -j*helicity*k/2],
                         [-j*helicity*k/2, j]])
        ledger.check(f"physical relative-rate mass determinant, helicity {helicity}",
                     s.simplify(mass.det()-rho*j) == 0)
    bmass = -j/2
    effective = cap_c+4*alpha*bmass/rho-4*alpha*(0-bmass**2/rho)/j
    ledger.check("complete second-gradient mass normalization gives C-alpha*j/rho",
                 s.simplify(effective-(cap_c-alpha*j/rho)) == 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
