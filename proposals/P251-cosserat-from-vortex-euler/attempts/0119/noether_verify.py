"""Exact current normalization and the genuine rigid-rotation alternative."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0119-Noether-current")
    t = s.Symbol("t", real=True)
    beta = s.Symbol("beta", nonzero=True, real=True)
    c1, c2, s1, s2 = s.symbols("c1 c2 s1 s2", real=True)
    omega = s.Matrix([[0, beta], [-beta, 0]])
    c, spin = s.Matrix([[c1, c2]]), s.Matrix([[s1, s2]])
    bracket = (c*(-omega.inv())*spin.T)[0]
    ledger.check("mechanical spin generates unit scalar angle only at the measured determinant",
                 s.simplify(bracket-c.col_join(spin).det()/beta) == 0)
    q, v, M, K = (s.Function(name)(t) for name in ("q", "v", "M", "K"))
    symmetric = -M*(q*v.diff(t)-v*q.diff(t))/2-K*q*q/2-M.diff(t)*q*v/2-M*v*v/2
    canonical = M*v*q.diff(t)-M*v*v/2-K*q*q/2
    ledger.check("the complete moving one-form has its exact boundary derivative",
                 s.simplify(symmetric+s.diff(M*q*v/2, t)-canonical) == 0)
    size, inertia, connection = (s.Function(name)(t) for name in ("size", "I", "connection"))
    angle = size*q
    normalized_spin = (inertia*angle.diff(t)+connection*angle)/size
    ledger.check("normalizing a physical axis retains inertia and changes only its connection",
                 s.simplify(normalized_spin-inertia*q.diff(t)
                            -(connection+inertia*size.diff(t)/size)*q) == 0)

    x, y, z = s.symbols("x y z", real=True)
    coords = (x, y, z)
    rho = s.Symbol("rho", positive=True)
    # Independent actual smooth Euler fixture; covariance supplies the general law.
    u = s.Matrix([s.sin(z)+2*s.cos(y), 3*s.sin(x)+s.cos(z),
                  2*s.sin(y)+3*s.cos(x)])
    pressure = -rho*u.dot(u)/2
    axis = s.Matrix([1, 2, 3])
    rotation = axis.cross(s.Matrix(coords))
    du = axis.cross(u)-u.jacobian(coords)*rotation
    dp = -(s.Matrix([s.diff(pressure, coordinate) for coordinate in coords]).dot(rotation))
    residual = du.jacobian(coords)*u+u.jacobian(coords)*du
    residual += s.Matrix([s.diff(dp, coordinate)/rho for coordinate in coords])
    ledger.check("a genuine full spatial rotation is a neutral linearized Euler direction",
                 all(s.trigsimp(value) == 0 for value in residual))
    ledger.check("the same rotation direction respects incompressibility",
                 s.trigsimp(sum(s.diff(du[i], coords[i]) for i in range(3))) == 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
