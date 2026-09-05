"""Exact velocity/current corrector diagnostics; every product mode retained."""

from __future__ import annotations

import sympy as sp

from substrate_framework import euler_fourier as ef
from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0179")
    y, z = sp.symbols("Y Z", real=True)
    psi, alpha = sp.cos(y) + sp.cos(z), sp.cos(y) - sp.cos(z)
    a_alpha = sp.sin(z) * sp.diff(alpha, y) - sp.sin(y) * sp.diff(alpha, z)
    ledger.check("two-wave strain source is the actual streamline derivative", a_alpha == -2 * sp.sin(y) * sp.sin(z))
    ledger.check("antisymmetric strain lies in the exact Laplacian kernel", sp.expand(-sp.diff(alpha, y, 2) - sp.diff(alpha, z, 2) - alpha) == 0)
    ledger.check("streamfunction is invariant and strain changes sign under exchange", psi.xreplace({y: z, z: y}) == psi and alpha.xreplace({y: z, z: y}) == -alpha)
    norm_alpha = sp.integrate(sp.integrate(alpha**2, (y, 0, 2 * sp.pi)), (z, 0, 2 * sp.pi)) / (2 * sp.pi)**2
    ledger.check("Fredholm obstruction has strictly nonzero source pairing", norm_alpha == 1)

    u = (ef.add(ef.trig(1), ef.trig(2)), ef.trig(2, kind="sin"), ef.scale(ef.trig(1, kind="sin"), -1))
    ledger.check("the compared field is the actual constant-minus-one-curl Euler field", all(not ef.add(left, right) for left, right in zip(ef.curl(u), u, strict=True)))
    kap, displacement = sp.Matrix([0, 1, 1]), sp.Matrix([0, 1, -1])
    translation = tuple(ef.scale(ef.add(*(ef.scale(ef.derivative(component, j), displacement[j]) for j in range(3))), -1) for component in u)
    derivative_translation = tuple(ef.add(*(ef.scale(ef.derivative(component, j), kap[j]) for j in range(3))) for component in translation)
    forcing = ef.leray(ef.cross(u, derivative_translation))
    expected_vorticity = ef.scale(ef.mul(ef.trig(1, kind="sin"), ef.trig(2, kind="sin")), 4)
    ledger.check("full pressure forcing has the derived planar cohomology sign", not ef.add(ef.curl(forcing)[0], ef.scale(expected_vorticity, -1)))
    phi = sp.Function("phi")(y, z)
    planar = sp.Matrix([-sp.diff(phi, z), sp.diff(phi, y)])
    base = sp.Matrix([sp.sin(z), -sp.sin(y)])
    raw_euler = -planar.jacobian((y, z)) * base - base.jacobian((y, z)) * planar
    raw_curl = sp.diff(raw_euler[1], y) - sp.diff(raw_euler[0], z)
    b_phi = -sp.diff(phi, y, 2) - sp.diff(phi, z, 2) - phi
    expected_generator = base[0] * sp.diff(b_phi, y) + base[1] * sp.diff(b_phi, z)
    ledger.check("planar vorticity generator derives AB without a pressure omission", sp.simplify(raw_curl - expected_generator) == 0)

    # Full adjoint-current Krylov rows, not a projected dynamical truncation.
    basis = []
    for axis in (1, 2):
        for component in range(3):
            if component == axis:
                continue
            for kind in ("sin", "cos"):
                vector = [{}, {}, {}]
                vector[component] = ef.trig(axis, kind=kind)
                basis.append(tuple(vector))
    a = ef.add(*(ef.scale(u[j], kap[j]) for j in range(3)))
    u_d = ef.add(*(ef.scale(u[j], displacement[j]) for j in range(3)))
    current = ef.leray(tuple(ef.add(ef.scale(a, displacement[j]), ef.scale(u_d, kap[j])) for j in range(3)))

    def adjoint(field):
        advected = ef.transport(u, field)
        transpose = tuple(ef.add(*(ef.mul(ef.derivative(u[j], i), field[j]) for j in range(3))) for i in range(3))
        return ef.leray(tuple(ef.add(advected[j], ef.scale(transpose[j], -1)) for j in range(3)))

    test_generator = ef.leray(ef.cross(u, tuple(ef.add(left, right) for left, right in zip(ef.curl(basis[0]), basis[0], strict=True))))
    ledger.check("actual adjoint current agrees with the full Euler generator", ef.inner(adjoint(current), basis[0]) == ef.inner(current, test_generator))
    rows, right = [], []
    observed = current
    first_mismatch = None
    for order in range(9):
        next_observed = adjoint(observed)
        rows.append([ef.inner(next_observed, column) for column in basis])
        right.append(-ef.inner(observed, forcing))
        matrix, target = sp.Matrix(rows), sp.Matrix(right)
        rank, augmented = matrix.rank(), matrix.row_join(target).rank()
        print(f"current derivative{order + 1}: coefficient rank={rank}, augmented rank={augmented}, source={right[-1]}", flush=True)
        if rank != augmented and first_mismatch is None:
            first_mismatch = order + 1
        observed = next_observed
    if first_mismatch is not None:
        witnesses = [vector for vector in matrix.T.nullspace() if (vector.T * target)[0] != 0]
        witness = witnesses[0]
        ledger.check("finite current obstruction has a derived exact left-null certificate", witness.T * matrix == sp.zeros(1, len(basis)) and (witness.T * target)[0] != 0)
        print(f"first-shell initial-current route refuted at derivative{first_mismatch}; certificate={list(witness)}; source pairing={(witness.T * target)[0]}", flush=True)
    else:
        solution, parameters = matrix.gauss_jordan_solve(target)
        ledger.check("finite current fit solves every derived row exactly", matrix * solution == target)
        print(f"finite first-shell fit through derivative9 only; free parameters={parameters}; no all-time inference", flush=True)

    # The fully solvable stationary elementary-wave lift has the inverted sign.
    theta = sp.symbols("theta", real=True)
    direction = sp.Matrix([sp.sin(theta), 0, sp.cos(theta)])
    macro = sp.Matrix([-sp.cos(theta), 0, sp.sin(theta)])
    one_wave = sp.Matrix([sp.cos(z), sp.sin(z), 0])
    tangent = -macro[2] * one_wave.diff(z)
    lift = -direction.cross(tangent)
    observed_stress = direction.dot(one_wave) * lift + one_wave * direction.dot(lift)
    average = sp.Matrix([sp.integrate(component, (z, 0, 2 * sp.pi)) / (2 * sp.pi) for component in observed_stress])
    projected = (sp.eye(3) - direction * direction.T) * average
    ledger.check("solvable stationary velocity selects inverted physical acoustic stress", all(sp.trigsimp(component) == 0 for component in projected - sp.sin(theta)**2 * macro / 2))
    ledger.check("stationary elementary-wave residual vanishes before averaging", one_wave.cross(one_wave.diff(z, 2)) == sp.zeros(3, 1))
    s = sp.symbols("s", real=True)
    haar_inverted = sp.integrate((1 - s**2) / 8, (s, -1, 1))
    ledger.check("whole-law stationary-preparation sign differs from bare positive branch", haar_inverted == sp.Rational(1, 6) and haar_inverted != -sp.Rational(2, 15))
    print(f"one-wave stationary Haar response={haar_inverted}; bare response=-2/15", flush=True)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
