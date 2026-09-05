"""Exact candidate solve: all products retained, no Fourier truncation."""

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0134-cell-ansatz")
    trig = [ef.trig(axis, kind=kind) for axis in range(3) for kind in ("sin", "cos")]
    coefficients = sp.symbols("w0:18", real=True)
    corrector = tuple(ef.add(*(ef.scale(trig[j], coefficients[6 * i + j]) for j in range(6))) for i in range(3))
    cases = [(1, (0, 0, 1), (1, 0, 0)), (0, (0, 0, 1), (1, 0, 0)),
             (1, (1, 1, 0), (1, -1, 0)), (1, (1, 1, 1), (1, -1, 0))]
    for amplitude, n, shift in cases:
        label = f"C={amplitude}, n={n}, U={shift}"
        u = (ef.add(ef.trig(2, kind="sin"), ef.scale(ef.trig(1), amplitude)),
             ef.add(ef.trig(0, kind="sin"), ef.trig(2)),
             ef.add(ef.scale(ef.trig(1, kind="sin"), amplitude), ef.trig(0)))
        ledger.check(f"{label}: exact stationary Beltrami field", not ef.divergence(u) and all(not ef.add(ef.curl(u)[j], ef.scale(u[j], -1)) for j in range(3)))
        tangent = tuple(ef.add(*(ef.scale(ef.derivative(component, j), -shift[j]) for j in range(3))) for component in u)
        energy = ef.scale(ef.add(*(ef.mul(component, component) for component in u)), sp.Rational(1, 2))
        pressure0 = ef.add(*(ef.scale(ef.derivative(energy, j), shift[j]) for j in range(3)))
        u_n = ef.add(*(ef.scale(u[j], n[j]) for j in range(3)))
        rhs = tuple(ef.scale(ef.add(ef.mul(u_n, tangent[j]), ef.scale(pressure0, n[j])), -1) for j in range(3))
        b = tuple(ef.add(ef.transport(u, corrector)[j], ef.transport(corrector, u)[j], ef.scale(rhs[j], -1)) for j in range(3))
        residual = ef.leray(b)
        divergence = ef.add(ef.divergence(corrector), *(ef.scale(tangent[j], n[j]) for j in range(3)))
        equations = [sp.re(v).expand() for field in (*residual, divergence) for v in field.values()]
        equations += [sp.im(v).expand() for field in (*residual, divergence) for v in field.values()]
        equations = [v for v in equations if v != 0]
        matrix, source = sp.linear_eq_to_matrix(equations, coefficients)
        solutions = sp.linsolve((matrix, source), coefficients)
        print(f"{label}: equations={matrix.rows}, variables={matrix.cols}, solution={solutions}")
        if solutions is sp.EmptySet:
            ledger.check(f"{label}: inconsistency is exposed by exact augmented rank", matrix.rank() < matrix.row_join(source).rank())
            print(f"{label}: route_scope=FIRST_SHELL_ONLY; broader cell space remains active")
        else:
            solution = next(iter(solutions))
            substitutions = dict(zip(coefficients, solution))
            ledger.check(f"{label}: candidate solves EVERY exact residual coefficient", all(sp.simplify(e.subs(substitutions)) == 0 for e in equations))
            solved = tuple({wave: sp.expand(value.subs(substitutions)) for wave, value in component.items()} for component in corrector)
            w_n = ef.add(*(ef.scale(solved[j], n[j]) for j in range(3)))
            stress = sp.Matrix([sp.simplify(ef.add(ef.mul(u[i], w_n), ef.mul(solved[i], u_n)).get(ef.ZERO, 0)) for i in range(3)])
            nv = sp.Matrix(n)
            projection = sp.eye(3) - nv * nv.T / nv.dot(nv)
            print(f"{label}: actual negative stress response={list(sp.simplify(-projection * stress))}")
    raise SystemExit(ledger.finish())


if __name__ == "__main__":
    main()
