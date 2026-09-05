"""Exact identities for the 0253 logarithmic flat-edge construction.

This is symbolic exposure only.  It neither solves the finite-radius free-boundary
problem nor supplies numerical evidence.
"""

import sympy as sp


d, a, rho, R, x, lam = sp.symbols(
    "d a rho R x lam", positive=True, finite=True
)

phi_d = sp.exp(-1 / d)
edge_ratio = sp.simplify(sp.diff(phi_d, d, 2) / phi_d)
assert sp.simplify(edge_ratio - (d**-4 - 2 * d**-3)) == 0

# For rho < a, d=a-rho and Phi(rho)=exp[-1/(a-rho)].
phi_rho = sp.exp(-1 / (a - rho))
laplace_2_ratio = sp.simplify(
    (sp.diff(phi_rho, rho, 2) + sp.diff(phi_rho, rho) / rho) / phi_rho
)
expected_laplace_2_ratio = (
    (a - rho) ** -4
    - 2 * (a - rho) ** -3
    - 1 / (rho * (a - rho) ** 2)
)
assert sp.simplify(laplace_2_ratio - expected_laplace_2_ratio) == 0

# If r=R+x and rho=sqrt(x^2+z^2), a circular cross-section obeys
# Delta*Phi=Phi''+R Phi'/(r rho).  Compare it with the straight radial
# Laplacian Phi''+Phi'/rho.
phi, phi_prime, laplace_2 = sp.symbols("phi phi_prime laplace_2")
r = R + x
laplace_star = laplace_2 - x * phi_prime / (r * rho)
g = lam**2 * phi + laplace_2
finite_R_residual = sp.factor(
    -laplace_star - lam**2 * phi + (r**2 / R**2) * g
)
expected_residual = sp.factor(
    x * phi_prime / (r * rho) + (2 * x / R + x**2 / R**2) * g
)
assert sp.simplify(finite_R_residual - expected_residual) == 0

# The reciprocal term x/(R+x) cannot be represented on a nontrivial circular
# level as A(rho)+B(rho)(R+x)^2.  Its third x derivative is nonzero.
reciprocal_third_derivative = sp.factor(sp.diff(x / (R + x), x, 3))
assert reciprocal_third_derivative != 0

# With Phi=exp(-1/T), the singular logarithmic equation becomes a degenerate
# but coefficient-regular equation in the defining function T.
T = sp.symbols("T", positive=True)
G_T = sp.simplify(
    T**4
    * (
        lam**2
        + T**-4
        - 2 * T**-3
        - T**-2 / (a - T)
    )
)
assert sp.simplify(
    G_T - (1 - 2 * T - T**2 / (a - T) + lam**2 * T**4)
) == 0
assert sp.limit(G_T, T, 0, dir="+") == 1

print("edge_ratio =", edge_ratio)
print("radial_laplacian_ratio =", laplace_2_ratio)
print("finite_R_residual =", finite_R_residual)
print("d3_dx3[x/(R+x)] =", reciprocal_third_derivative)
print("regularized_G(T) =", G_T)
print("regularized_G(0+) =", sp.limit(G_T, T, 0, dir="+"))
print("EXACT symbolic checks: PASS")
