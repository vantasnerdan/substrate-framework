#!/usr/bin/env python3
"""Exact algebraic regressions for P253/0015; no numerical approximation."""

from __future__ import annotations

import json

import sympy as sp


rho, z, c = sp.symbols("rho z c", real=True)


def bracket(a: sp.Expr, b: sp.Expr) -> sp.Expr:
    return sp.diff(a, rho) * sp.diff(b, z) - sp.diff(a, z) * sp.diff(b, rho)


# Translating-frame generator: {q,Psi-c rho}={q,Psi}+c q_z.
q = rho**2 * z + sp.sin(z) + rho
psi = rho * z**2 + sp.cos(rho)
frame_identity = sp.simplify(
    bracket(q, psi - c * rho) - bracket(q, psi) - c * sp.diff(q, z)
)
assert frame_identity == 0


# Poisson integration-by-parts identity used in the first variation.  Its
# pointwise difference is a divergence, whose integral vanishes under the
# compact-support/periodic hypotheses in construction.md.
chi = rho**3 + rho * z + z**2
cyclic_difference = (
    psi * bracket(q, chi)
    - q * bracket(chi, psi)
    - sp.diff(q * psi * sp.diff(chi, z), rho)
    + sp.diff(q * psi * sp.diff(chi, rho), z)
)
assert sp.simplify(cyclic_difference) == 0


# Opposite high-time-frequency signs in the loop-action Hessian.
t, period = sp.symbols("t period", positive=True)
k = sp.symbols("k", integer=True, positive=True)
n = 2 * sp.pi * k / period
J = sp.Matrix([[0, -1], [1, 0]])
e = sp.Matrix([1, 0])
f = sp.Matrix([0, -1])


def omega(u: sp.Matrix, v: sp.Matrix) -> sp.Expr:
    return (u.T * J * v)[0]


eta_plus = e * sp.cos(n * t) + f * sp.sin(n * t)
eta_minus = e * sp.cos(n * t) - f * sp.sin(n * t)
integrand_plus = sp.trigsimp(omega(eta_plus, sp.diff(eta_plus, t)))
integrand_minus = sp.trigsimp(omega(eta_minus, sp.diff(eta_minus, t)))
assert sp.simplify(integrand_plus - n) == 0
assert sp.simplify(integrand_minus + n) == 0
omega_integral_plus = 2 * sp.pi * k
omega_integral_minus = -2 * sp.pi * k
action_principal_plus = -omega_integral_plus
action_principal_minus = -omega_integral_minus


# With i_X Omega=dH and -dTheta=Omega, variation of
# int Theta(qdot)-H is -Omega(eta,qdot)-dH(eta).
g1, g2, e1, e2 = sp.symbols("g1 g2 e1 e2", real=True)
gradient = sp.Matrix([g1, g2])
test_tangent = sp.Matrix([e1, e2])
hamiltonian_vector = J * gradient
assert sp.simplify(omega(hamiltonian_vector, test_tangent) - gradient.dot(test_tangent)) == 0
action_variation_on_flow = sp.simplify(
    -omega(test_tangent, hamiltonian_vector) - gradient.dot(test_tangent)
)
assert action_variation_on_flow == 0


# Period-twist monodromy in the (phase, amplitude) basis.
Tprime = sp.symbols("Tprime", real=True, nonzero=True)
M = sp.Matrix([[1, -Tprime], [0, 1]])
N = M - sp.eye(2)
assert M.T * J * M == J
assert sp.expand(M.charpoly().as_expr()) == sp.Symbol("lambda") ** 2 - 2 * sp.Symbol("lambda") + 1
assert N != sp.zeros(2)
assert N**2 == sp.zeros(2)


print(
    json.dumps(
        {
            "frame_generator_identity": "pass",
            "periodic_poisson_cyclic_identity": "pass",
            "local_action_sign_convention": "i_X Omega=dH, -dTheta=Omega",
            "action_variation_on_hamiltonian_flow": "pass",
            "loop_action_principal_high_frequency_signs": [
                str(action_principal_plus),
                str(action_principal_minus),
            ],
            "period_twist_monodromy_symplectic": "pass",
            "period_twist_jordan_nilpotency": "N != 0 and N^2 = 0",
        },
        indent=2,
        sort_keys=True,
    )
)
