"""Exact joint residual (0241 equations 4-8); no Euler supplier is assumed."""
import sympy as s

from substrate_framework.micropolar import (
    MicropolarCoefficients,
    micropolar_fourier_stiffness,
)
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger('P251-0241-actual-joint-residual')
    rho, j, nu, a, ct, cl = s.symbols('rho j nu a cT cL', positive=True)
    k, t = s.symbols('k t', real=True)
    alpha, beta = j*nu**2/4, j/(2*rho)
    curl = s.I*s.Matrix([[0, -k, 0], [k, 0, 0], [0, 0, 0]])
    transform = s.BlockMatrix([[s.eye(3), -beta*curl], [curl/2, s.eye(3)]]).as_explicit()
    idx = [0, 1, 3, 4, 5]
    transform = transform.extract(idx, idx)
    mass = s.diag(rho, rho, j, j, j)
    freq = s.diag(a*k**2, a*k**2, nu**2+ct*k**2, nu**2+ct*k**2, nu**2+cl*k**2)
    gt, gl = j*(ct-alpha/rho), j*cl
    coeff = MicropolarCoefficients(0, rho*a, alpha, gl/2, 0, gt)
    stiff = s.Matrix(micropolar_fourier_stiffness([0, 0, k], coeff)).extract(idx, idx)
    defect = s.simplify(stiff*transform-mass*transform*freq)
    checks.check('Eq6: all physical joint mismatch coefficients below degree3 vanish',
                 defect.applyfunc(lambda x: s.series(x, k, 0, 3).removeO()) == s.zeros(5))
    expected = s.zeros(5)
    expected[:2, 2:] = beta*(rho*ct-rho*a-alpha)*k**2*curl[:2, :]
    expected[2:, :2] = j*(ct-a-alpha/rho)*k**2*curl[:, :2]/2
    checks.check('Eq6: exact cubic mismatch including both off-diagonal rows', s.simplify(defect-expected) == s.zeros(5))
    z = s.Matrix([s.Function(f'z{i}')(t) for i in range(5)])
    e = s.Matrix([s.Function(f'e{i}')(t) for i in range(5)])
    physical = transform*z+e
    residual = mass*physical.diff(t, 2)+stiff*physical
    rhs = mass*transform*(z.diff(t, 2)+freq*z)+defect*z+mass*e.diff(t, 2)+stiff*e
    checks.check('Eq5: actual history and physical observation errors both remain',
                 s.simplify(residual-rhs) == s.zeros(5, 1))
    mutated = transform.copy()
    mutated[:2, 2:] *= 2
    wrong = s.simplify(stiff*mutated-mass*mutated*freq)
    checks.check('wrong optical current normalization exposes a linear coupling defect',
                 wrong.applyfunc(lambda x: s.expand(x).coeff(k, 1)) != s.zeros(5))
    q = s.Function('q')(t)
    U = s.Function('U')(t)
    # One arbitrary tensor component of Q; the identity is componentwise.
    Q = q*s.diff(U, t)
    x = s.symbols('x:3', real=True)
    velocity = s.Matrix([s.Function(f'u{i}')(t, *x) for i in range(3)])
    tensor = s.Matrix(3, 3, lambda i, jj: q*sum(
        s.LeviCivita(i, jj, kk)*velocity[kk] for kk in range(3)))
    spin = s.Matrix([s.Function(f'S{i}')(t, *x) for i in range(3)])
    mu = s.Matrix(3, 3, lambda i, jj: s.Function(f'm{i}{jj}')(t, *x))
    div = lambda m: s.Matrix([sum(s.diff(m[i, jj], x[jj]) for jj in range(3))
                             for i in range(3)])
    improved_spin = spin-div(tensor)
    improved_mu = mu-tensor.diff(t)
    checks.check('Eq7: derived superpotential leaves the full angular divergence unchanged',
                 s.simplify(improved_spin.diff(t)-div(improved_mu)
                            -spin.diff(t)+div(mu)) == s.zeros(3, 1))
    checks.check('omitting q-dot flux leaves a nonzero angular divergence defect',
                 s.simplify(improved_spin.diff(t)-div(mu-q*tensor.diff(t)/q
                                      +s.diff(q,t)*tensor/q)
                            -spin.diff(t)+div(mu)) != s.zeros(3, 1))
    # Use an independent integration variable, preserving the endpoint.
    r = s.Symbol('r', real=True)
    primitive = q*U-q.subs(t, 0)*U.subs(t, 0)-s.Integral(
        s.diff(q.subs(t, r), r)*U.subs(t, r), (r, 0, t))
    checks.check('Eq8: complete accumulated current differentiates to the physical flux',
                 s.simplify(s.diff(primitive, t)-Q) == 0)
    checks.check('Eq8: initial accumulated correction vanishes without deleting initial G',
                 s.simplify(primitive.subs(t, 0)) == 0)
    checks.check('omitting q-dot memory leaves an exposed nonzero current source',
                 s.simplify(s.diff(q*U, t)-Q) == s.diff(q, t)*U
                 and s.diff(q, t)*U != 0)
    return checks.finish()


if __name__ == '__main__':
    raise SystemExit(main())
