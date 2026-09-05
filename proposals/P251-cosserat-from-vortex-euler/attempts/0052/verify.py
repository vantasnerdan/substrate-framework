"""Exact centering, boundary-current and fixed-response Gram identities."""

import sympy as s

from substrate_framework.euler_orbit import reduce_euler_rotor_block
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0052-material-spin-joining")
    m1, m2 = s.symbols("m1 m2", positive=True)
    v = s.Matrix(s.symbols("V0:3", real=True))
    w = s.Matrix(s.symbols("w0:3", real=True))
    # Two parcels already exhibit the exact weighted mean constraint.
    w2 = -m1*w/m2
    kinetic = (m1*(v+w).dot(v+w)+m2*(v+w2).dot(v+w2))/2
    centered = (m1+m2)*v.dot(v)/2+(m1*w.dot(w)+m2*w2.dot(w2))/2
    ledger.check("mass-centroid decomposition removes the entire kinetic cross",
                 s.expand(kinetic-centered) == 0)
    ledger.check("canonical translation momentum is physical mass times centroid velocity",
                 s.simplify(s.Matrix([s.diff(kinetic, vi) for vi in v])-(m1+m2)*v) == s.zeros(3, 1))

    x, y, z, time = s.symbols("x y z t", real=True)
    coords = (x, y, z)
    r = s.Matrix(coords)
    u = s.Matrix([s.Function(f"u{i}")(*coords) for i in range(3)])

    def curl(f):
        return s.Matrix([s.diff(f[2], y)-s.diff(f[1], z),
                         s.diff(f[0], z)-s.diff(f[2], x),
                         s.diff(f[1], x)-s.diff(f[0], y)])

    def div_tensor(tensor):
        return s.Matrix([sum(s.diff(tensor[i, j], coords[j]) for j in range(3)) for i in range(3)])

    def ax(tensor):
        return s.Matrix([sum(s.LeviCivita(i, j, k)*tensor[j, k]
                             for j in range(3) for k in range(3)) for i in range(3)])

    ledger.check("finite-cell angular impulse retains its exact boundary curl term",
                 s.simplify(curl(r.dot(r)*u)-2*r.cross(u)-r.dot(r)*curl(u)) == s.zeros(3, 1))
    boundary = s.Matrix([s.Function(f"B{i}")(x, y, z, time) for i in range(3)])
    improvement = s.Matrix(3, 3, lambda i, j: sum(s.LeviCivita(i, j, k)*boundary[k]
                                                for k in range(3))/2)
    ledger.check("antisymmetric improvement has exactly the boundary spin as axial vector",
                 ax(improvement) == boundary)
    ledger.check("canonical momentum correction is half the curl of boundary spin",
                 s.simplify(div_tensor(improvement)-curl(boundary)/2) == s.zeros(3, 1))
    ledger.check("boundary-spin momentum improvement preserves mass continuity",
                 s.simplify(sum(s.diff(div_tensor(improvement)[i], coords[i]) for i in range(3))) == 0)
    angular_flux = s.Matrix.hstack(*(r.cross(improvement[:, j]) for j in range(3)))
    ledger.check("total angular momentum changes only by the explicit boundary flux",
                 s.simplify(r.cross(div_tensor(improvement))-boundary-div_tensor(angular_flux)) == s.zeros(3, 1))
    ledger.check("improved spin balance retains the stress-axial source exactly",
                 s.simplify(-boundary.diff(time)+ax(improvement.diff(time))) == s.zeros(3, 1))
    ledger.check("time derivative commutes with the force-current improvement",
                 s.simplify(div_tensor(improvement).diff(time)-div_tensor(improvement.diff(time))) == s.zeros(3, 1))

    # Fixed-response pseudoinverse on a genuinely rank-deficient range.
    response = s.Matrix([[1, 2], [2, 4], [0, 0]])
    gram = response*response.T
    coefficients = s.Matrix(s.symbols("a0:2"))
    mean = response*coefficients
    inverse = gram.pinv()
    ledger.check("rank-deficient response has an exact Moore-Penrose inverse",
                 gram*inverse*gram == gram and inverse*gram*inverse == inverse)
    ledger.check("fixed-response projection cancels every actual mean component",
                 s.simplify(mean-gram*inverse*mean) == s.zeros(3, 1))
    ledger.check("rank-three invertibility is not assumed for centering", gram.rank() == 1)
    l_parallel = response.T*s.Matrix([1, 0, 0])
    gram_cross = response*l_parallel
    residual_parallel = (l_parallel.T*l_parallel)[0]-(gram_cross.T*inverse*gram_cross)[0]
    ledger.check("a nonzero common functional can be lost under centering",
                 l_parallel != s.zeros(2, 1) and residual_parallel == 0)
    extension = s.Matrix([1, 0, 1])
    response_extended = response.row_join(s.zeros(3, 1))
    cross_extended = response_extended*extension
    residual = (extension.T*extension)[0]-(cross_extended.T*inverse*cross_extended)[0]
    ledger.check("independent angular response leaves a strictly positive residual Gram norm",
                 residual.is_positive)

    # Call the shared Euler-orbit API instead of duplicating its reduction.
    hessian = s.Matrix([[5, 1, 1], [1, 4, 1], [1, 1, 3]])
    b, c = s.symbols("b c", positive=True)
    reduced = reduce_euler_rotor_block(hessian, b, c)
    p = hessian.extract([0, 2], [0, 2])
    coupling = hessian.extract([0, 2], [1])
    rates = s.Matrix(s.symbols("Bdot qdot", real=True))
    q = s.Symbol("q", real=True)
    d = s.diag(b, c)
    plus = p.inv()*(d*rates-coupling*q)
    minus = p.inv()*(-d*rates-coupling*q)
    tr, ts = s.symbols("surface_r surface_s", real=True)
    surface_row = s.Matrix([[tr, ts]])
    paired_surface = (surface_row*plus-surface_row*minus)[0]/2
    ledger.check("paired surface spin retains its complete Euler-derived rate response",
                 s.simplify(paired_surface-(surface_row*p.inv()*d*rates)[0]) == 0)
    paired_canonical = b*(plus[0]-minus[0])/2
    ledger.check("canonical bulk spin equals the shared API's common-angle momentum",
                 s.simplify(paired_canonical-(reduced.kinetic*rates)[0]) == 0)
    ledger.check("surface spin generally differs from canonical bulk spin",
                 s.diff(paired_surface, tr) != 0 and s.diff(paired_surface, ts) != 0)
    print("Exact action/current bridge established; selected EPS centered-response rank remains explicit.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
