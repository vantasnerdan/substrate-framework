"""Exact exposing oracle for P253/0079 finite-dimensional claims."""

import importlib.util
from pathlib import Path

import sympy as sp


MODULE_PATH = Path(__file__).with_name("schwinger_hopf_blocks.py")
SPEC = importlib.util.spec_from_file_location("p253_0079_blocks", MODULE_PATH)
assert SPEC and SPEC.loader
blocks = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(blocks)


J = blocks.standard_complex_structure()
assert sp.simplify(J * J + sp.eye(4)) == sp.zeros(4)

# A generic real operator is decomposed into commuting and anticommuting rows.
entries = sp.symbols("v0:16", real=True)
V = sp.Matrix(4, 4, entries)
VC, VA = blocks.complex_linear_split(V, J)
assert sp.simplify(VC * J - J * VC) == sp.zeros(4)
assert sp.simplify(VA * J + J * VA) == sp.zeros(4)
assert sp.simplify(VC + VA - V) == sp.zeros(4)

# A complex-linear generator realifies to a matrix commuting with J.
a, b, c, d = sp.symbols("a b c d", real=True)
H = sp.Matrix([[a, c - sp.I * d], [c + sp.I * d, b]])
C = -sp.I * H
C_real = blocks.realify_complex(C)
assert sp.simplify(C_real * J - J * C_real) == sp.zeros(4)

# Exact Pauli axes: diagonal, cosine and sine responses are non-collinear.
sigma_x = sp.Matrix([[0, 1], [1, 0]])
sigma_y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
sigma_z = sp.Matrix([[1, 0], [0, -1]])
assert blocks.pauli_noncollinearity(sigma_x, sigma_y) == sp.Matrix([0, 0, 1])
assert blocks.pauli_noncollinearity(sigma_z, sigma_x) == sp.Matrix([0, 1, 0])

# Collinear controls are explicitly rejected.
scale = sp.symbols("scale", real=True)
assert blocks.pauli_noncollinearity(sigma_x, scale * sigma_x) == sp.zeros(3, 1)

# The corrected two-sided gate-time bounds distinguish outbound leakage from
# Q-to-P return and have the exact zero-error and leading limits.
eps_a, eps_out, kappa_back = sp.symbols(
    "eps_A eps_out kappa_back", nonnegative=True
)
beta_bound, p_gain, q_gain, action_change = blocks.coupled_gate_bounds(
    eps_a, eps_out, kappa_back
)
assert beta_bound.subs(eps_a, 0) == 0
assert p_gain.subs({eps_a: 0, kappa_back: 0}) == 1
assert q_gain.subs({eps_a: 0, eps_out: 0, kappa_back: 0}) == 0
assert action_change.subs({eps_a: 0, kappa_back: 0}) == 0
assert sp.diff(beta_bound, eps_a).subs(eps_a, 0) == 1
assert sp.diff(action_change, eps_a).subs(
    {eps_a: 0, kappa_back: 0}
) == 2
assert sp.diff(action_change, kappa_back).subs(
    {eps_a: 0, kappa_back: 0}
) == 2
assert sp.diff(q_gain, eps_out).subs(
    {eps_a: 0, eps_out: 0, kappa_back: 0}
) == 1
# One-way leakage alone cannot perturb the P-action bound.
assert sp.diff(action_change, eps_out) == 0

# On a rational massive ray, the O(N^-2) delta shift cancels the O(N^-1)
# curvature mismatch against the Theta(N) derivative at leading C1 order.
k_star, p, large_n, delta_c, f_prime = sp.symbols(
    "k_star p N Delta_c F_prime", nonzero=True
)
delta_0, delta_shift = blocks.rational_ray_crossing_shift(
    k_star, p, large_n, delta_c, f_prime
)
assert delta_0 == k_star / (large_n * p)
linearized_mismatch = (
    large_n * p * f_prime * delta_shift + delta_0 * delta_c
)
assert sp.simplify(linearized_mismatch) == 0

# The stronger existence argument needs only C0 convergence: at the IVT
# endpoints the first physical wave number is exactly k_star +/- sqrt(e).
error_envelope = sp.symbols("e_bar", positive=True)
delta_center, delta_minus, delta_plus = blocks.massive_bracketing_interval(
    k_star, p, large_n, error_envelope
)
assert sp.simplify(large_n * p * delta_center - k_star) == 0
assert sp.simplify(
    large_n * p * delta_minus - (k_star - sp.sqrt(error_envelope))
) == 0
assert sp.simplify(
    large_n * p * delta_plus - (k_star + sp.sqrt(error_envelope))
) == 0
assert sp.simplify(error_envelope / sp.sqrt(error_envelope)) == sp.sqrt(
    error_envelope
)

print("PASS V_C commutes with J and V_A anticommutes with J")
print("PASS realified complex-linear Hamiltonian generator commutes with J")
print("PASS sigma_x/sigma_y and sigma_z/sigma_x Pauli axes are non-collinear")
print("PASS collinear Pauli controls are rejected")
print("PASS two-sided gate bounds separate outbound leakage and return feedback")
print("PASS rational-ray O(N^-2) shift cancels the leading curvature mismatch")
print("PASS C0 rational-ray bracket has error o(h) and width 2h/(Np)")
