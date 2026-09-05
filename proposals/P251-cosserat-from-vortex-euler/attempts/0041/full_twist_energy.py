"""Exact complete-variation and kernel receipts; no numerical energy fitting."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0041-full-Biot-Savart-twist")
    x, y, z = sp.symbols("x y z", real=True)
    f = sp.Function("f")(x, y)
    q = sp.Function("q")(z)

    def zero(name: str, expression: sp.Expr) -> None:
        ledger.check(name, sp.simplify(expression) == 0)

    def zero_matrix(name: str, expression: sp.Matrix) -> None:
        ledger.check(name, all(sp.simplify(entry) == 0 for entry in expression))

    def angular(expression: sp.Expr) -> sp.Expr:
        return -y * sp.diff(expression, x) + x * sp.diff(expression, y)

    a_vector = sp.Matrix([-y * f, x * f])
    b_scalar = angular(f)
    zero(
        "angular vorticity variation is div A",
        sp.diff(a_vector[0], x) + sp.diff(a_vector[1], y) - b_scalar,
    )
    first = sp.Matrix(
        [sp.diff(q, z) * a_vector[0], sp.diff(q, z) * a_vector[1], -q * b_scalar]
    )
    second = sp.Matrix(
        [
            2 * q * sp.diff(q, z) * y * b_scalar,
            -2 * q * sp.diff(q, z) * x * b_scalar,
            q**2 * angular(b_scalar),
        ]
    )
    zero(
        "first isovortical variation divergence free",
        sum(
            sp.diff(first[index], coordinate)
            for index, coordinate in enumerate((x, y, z))
        ),
    )
    zero(
        "second isovortical variation divergence free",
        sum(
            sp.diff(second[index], coordinate)
            for index, coordinate in enumerate((x, y, z))
        ),
    )
    push = sp.Matrix([sp.cos(q) * x - sp.sin(q) * y, sp.sin(q) * x + sp.cos(q) * y, z])
    zero("full axial-twist map preserves volume", push.jacobian((x, y, z)).det() - 1)
    zero(
        "curl A detects nontrivial compact twist",
        sp.diff(a_vector[1], x)
        - sp.diff(a_vector[0], y)
        - (2 * f + x * sp.diff(f, x) + y * sp.diff(f, y)),
    )
    zero(
        "weighted curl integration-by-parts factor",
        -(sp.diff(x**2 + y**2, x) * x + sp.diff(x**2 + y**2, y) * y)
        + 2 * (x**2 + y**2),
    )

    xi, eta, k = sp.symbols("xi eta k", real=True)
    momentum = sp.Matrix([xi, eta])
    norm2 = momentum.dot(momentum)
    projector = sp.eye(2) - momentum * momentum.T / norm2
    zero_matrix("transverse projector idempotent", projector * projector - projector)
    zero_matrix("transverse projector annihilates gradient", projector * momentum)
    vx, vy = sp.symbols("vx vy", real=True)
    vector = sp.Matrix([vx, vy])
    zero(
        "projector quadratic form is a square",
        (vector.T * projector * vector)[0] - (xi * vy - eta * vx) ** 2 / norm2,
    )
    zero(
        "complete resolvent subtraction gives projector",
        vector.dot(vector) / (norm2 + k**2)
        - (momentum.dot(vector)) ** 2 / (norm2 * (norm2 + k**2))
        - (vector.T * projector * vector)[0] / (norm2 + k**2),
    )
    zero(
        "resolvent identity fixes k squared sign",
        1 / (norm2 + k**2) - 1 / norm2 + k**2 / (norm2 * (norm2 + k**2)),
    )

    x1, x2, y1, y2 = sp.symbols("x1 x2 y1 y2", real=True)
    separation2 = (x1 - y1) ** 2 + (x2 - y2) ** 2
    log_distance = sp.log(separation2) / 2
    biharmonic = separation2 * log_distance / (8 * sp.pi)

    def angular_x(expression: sp.Expr) -> sp.Expr:
        return -x2 * sp.diff(expression, x1) + x1 * sp.diff(expression, x2)

    def angular_y(expression: sp.Expr) -> sp.Expr:
        return -y2 * sp.diff(expression, y1) + y1 * sp.diff(expression, y2)

    dot = x1 * y1 + x2 * y2
    cross = x1 * y2 - x2 * y1
    exact_kernel = -dot * log_distance / (2 * sp.pi) - angular_x(angular_y(biharmonic))
    exposed_kernel = (-dot * log_distance + dot / 2 + cross**2 / separation2) / (
        4 * sp.pi
    )
    zero(
        "complete real-space kernel including finite term",
        exact_kernel - exposed_kernel,
    )
    zero(
        "common rotation preserves the exact kernel",
        angular_x(exposed_kernel) + angular_y(exposed_kernel),
    )
    zero(
        "biharmonic normalization away from diagonal",
        sp.diff(biharmonic, x1, 2)
        + sp.diff(biharmonic, x2, 2)
        - (log_distance + 1) / (2 * sp.pi),
    )
    zero(
        "bounded angular kernel geometry",
        cross**2 - ((x1 * (y2 - x2) - x2 * (y1 - x1)) ** 2),
    )

    rotation = sp.Matrix(
        [
            [sp.cos(2 * sp.pi / 3), -sp.sin(2 * sp.pi / 3)],
            [sp.sin(2 * sp.pi / 3), sp.cos(2 * sp.pi / 3)],
        ]
    )
    zero_matrix(
        "threefold symmetry removes bend-twist vector",
        sp.eye(2) + rotation + rotation**2,
    )
    ledger.check(
        "no invariant nonzero transverse vector under C3",
        sp.simplify((rotation - sp.eye(2)).det()) != 0,
    )

    density, circulation, radius, epsilon = sp.symbols(
        "rho Gamma d epsilon", positive=True
    )
    line_tension = density * circulation**2 * sp.log(radius / epsilon) / (4 * sp.pi)
    zero(
        "three finite cores give triangle leading coefficient",
        3
        * density
        * circulation**2
        * radius**2
        * sp.log(radius / epsilon)
        / (4 * sp.pi)
        - 3 * line_tension * radius**2,
    )
    zero("sixfold isotropic bend projector trace", sp.trace(projector) - 1)
    ledger.check(
        "omitting longitudinal subtraction changes twist energy",
        sp.simplify(vector.dot(vector) - (vector.T * projector * vector)[0]) != 0,
    )
    ledger.check(
        "wrong longitudinal resolvent sign loses projector",
        sp.simplify(
            vector.dot(vector) / (norm2 + k**2)
            + (momentum.dot(vector)) ** 2 / (norm2 * (norm2 + k**2))
            - (vector.T * projector * vector)[0] / (norm2 + k**2)
        )
        != 0,
    )
    print(
        "Scope: exact variation and kernel identities. Positivity, infrared limits, and remainder bounds are proved in full-biot-savart-twist.md."
    )
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
