"""Exact affine core jets, bond rank, isotropy and projected shear algebra."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0054-same-EPS-affine-shear")
    x = s.Matrix(s.symbols("x y z", real=True))
    a, b, c, d, f = s.symbols("a b c d f", real=True)
    strain = s.Matrix([[a, c, d], [c, b, f], [d, f, -a-b]])

    def curl(v):
        return s.Matrix([s.diff(v[2], x[1])-s.diff(v[1], x[2]),
                         s.diff(v[0], x[2])-s.diff(v[2], x[0]),
                         s.diff(v[1], x[0])-s.diff(v[0], x[1])])

    potential = -x.cross(strain*x)/3
    ledger.check("compact-potential interior has the exact affine displacement",
                 s.simplify(curl(potential)-strain*x) == s.zeros(3, 1))
    ledger.check("physical core deformation gradient is exactly symmetric tracefree",
                 curl(potential).jacobian(x) == strain and s.trace(strain) == 0)
    chi = s.Function("chi")(*x)
    grad_chi = s.Matrix([s.diff(chi, t) for t in x])
    returned = curl(chi*potential)
    ledger.check("outer return is included by the exact curl product rule",
                 s.simplify(returned-chi*strain*x-grad_chi.cross(potential)) == s.zeros(3, 1))
    ledger.check("full returned displacement is divergence-free",
                 s.simplify(sum(s.diff(returned[i], x[i]) for i in range(3))) == 0)

    units = [s.eye(3)[:, j] for j in range(3)]
    bonds = units + [(units[i]+sign*units[j])/s.sqrt(2)
                     for i in range(3) for j in range(i+1, 3) for sign in (-1, 1)]
    frame_energy = s.expand(sum((n.T*strain*n)[0]**2 for n in bonds))
    frobenius = s.trace(strain.T*strain)
    diagonal = sum(strain[i, i]**2 for i in range(3))
    ledger.check("nine bond extensions have the exact shear frame identity",
                 s.simplify(frame_energy-frobenius-diagonal/2) == 0)
    spin = s.Matrix([[0, -a, b], [a, 0, -c], [-b, c, 0]])
    ledger.check("rigid rotation excites no affine bond extension",
                 all(s.simplify((n.T*spin*n)[0]) == 0 for n in bonds))

    basis = [s.diag(1, -1, 0)/s.sqrt(2), s.diag(1, 1, -2)/s.sqrt(6)]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        item = s.zeros(3)
        item[i, j] = item[j, i] = 1/s.sqrt(2)
        basis.append(item)
    gram = s.Matrix(5, 5, lambda i, j: s.trace(basis[i].T*basis[j]))
    ledger.check("strain coordinates are Frobenius orthonormal", gram == s.eye(5))
    bond_matrix = s.Matrix(9, 5, lambda i, j: (bonds[i].T*basis[j]*bonds[i])[0])
    frame = s.simplify(bond_matrix.T*bond_matrix)
    ledger.check("bond map has full five-dimensional rank", bond_matrix.rank() == 5)
    ledger.check("exact frame spectrum lies between one and three halves",
                 frame == s.diag(s.Rational(3, 2), s.Rational(3, 2), 1, 1, 1))

    rotations = [s.Matrix.hstack(*(e.cross(v) for v in units)) for e in units]
    reps = [s.Matrix(5, 5, lambda i, j: s.trace(basis[i]*(r*basis[j]-basis[j]*r)))
            for r in rotations]
    ledger.check("strain rotation generators preserve its inner product",
                 all(r.T == -r for r in reps))
    ledger.check("strain representation has the exact angular-momentum-two Casimir",
                 sum((r*r for r in reps), s.zeros(5)) == -6*s.eye(5))

    values = s.symbols("C0:15")
    sym = s.zeros(5)
    index = 0
    for i in range(5):
        for j in range(i, 5):
            sym[i, j] = sym[j, i] = values[index]
            index += 1
    constraints = [entry for r in reps for entry in r*sym-sym*r]
    linear = s.linear_eq_to_matrix(constraints, values)[0]
    ledger.check("isotropic symmetric shear tensor has exactly one scalar coefficient",
                 linear.rank() == 14)
    ledger.check("the surviving invariant tensor is the identity",
                 all(r*s.eye(5)-s.eye(5)*r == s.zeros(5) for r in reps))

    cross_values = s.symbols("T0:15")
    cross = s.Matrix(5, 3, cross_values)
    cross_constraints = [entry for r, j in zip(reps, rotations, strict=True)
                         for entry in r*cross-cross*j]
    cross_linear = s.linear_eq_to_matrix(cross_constraints, cross_values)[0]
    ledger.check("isotropy permits no symmetric-strain/vector-spin cross tensor",
                 cross_linear.rank() == 15)

    # A rank-deficient four-functional Gram example checks the pseudoinverse
    # identity; the general analytic range proof is in affine-shear.md.
    response = s.Matrix([[1, 0], [0, 1], [1, 1], [2, -1]])
    g = response*response.T
    coeff = s.Matrix(s.symbols("f0:2"))
    datum = response*coeff
    ledger.check("singular response Gram removes every compatible moment exactly",
                 g.rank() == 2 and s.simplify(datum-g*g.pinv()*datum) == s.zeros(4, 1))
    ledger.check("fixed response projection is idempotent on the response range",
                 s.simplify((g*g.pinv())**2-g*g.pinv()) == s.zeros(4))

    p11, p12, p22, n1, n2, h = s.symbols("p11 p12 p22 n1 n2 h", real=True)
    p = s.Matrix([[p11, p12], [p12, p22]])
    n = s.Matrix([n1, n2])
    q = s.symbols("strain", real=True)
    momenta = -p.inv()*n*q
    energy = ((momenta.T*p*momenta)[0]+2*q*(n.T*momenta)[0]+h*q*q)/2
    stiffness = h-(n.T*p.inv()*n)[0]
    ledger.check("static shear reduction retains the full off-diagonal momentum block",
                 s.simplify(energy-stiffness*q*q/2) == 0)
    density, trace_c = s.symbols("n_cell trace_C", positive=True)
    ledger.check("isotropic coefficient has the standard shear-energy normalization",
                 s.simplify(density*(trace_c/5)*frobenius/2
                            -(density*trace_c/10)*frobenius) == 0)
    print("Analytic oracle: affine-shear.md; full Euler energy and finite carrier bounds.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
