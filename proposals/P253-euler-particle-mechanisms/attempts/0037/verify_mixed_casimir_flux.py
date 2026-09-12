"""Exact symbolic checks for the local mixed-Casimir conservation law."""

import sympy as sp


r, z = sp.symbols("r z", positive=True)
psi = sp.Function("Psi")(r, z)
xi = sp.Function("xi")(r, z)
zeta = sp.Function("zeta")(r, z)
D = sp.Function("D")
G = sp.Function("G")


def bracket(f, g):
    return (sp.diff(f, r) * sp.diff(g, z)
            - sp.diff(f, z) * sp.diff(g, r)) / r


# Use the material derivative directly, with G'(xi)=xi D(xi).
xi_t = -bracket(psi, xi)
zeta_t = -bracket(psi, zeta) + 2 * xi * sp.diff(xi, z) / r**4
q_t = zeta_t * D(xi) + zeta * sp.diff(D(xi), xi) * xi_t
transport_residual = sp.expand(q_t + bracket(psi, zeta * D(xi)))
expected = 2 * xi * D(xi) * sp.diff(xi, z) / r**4
assert sp.simplify(transport_residual - expected) == 0

# The meridional bracket is an axial flux plus a radial boundary flux.
q = sp.Function("q")(r, z)
divergence_residual = (
    r * bracket(psi, q)
    - sp.diff(sp.diff(psi, r) * q, z)
    + sp.diff(sp.diff(psi, z) * q, r)
)
assert sp.simplify(divergence_residual) == 0

# For a regular-label cutoff D(L(r)), the density is a radial divergence plus
# a bounded interior energy pairing and the axial second derivative.
stream = sp.Function("stream")(r, z)
label = sp.Function("L")(r)
weighted = D(label)
laplace_star = sp.diff(stream, r, 2) - sp.diff(stream, r) / r + sp.diff(stream, z, 2)
regular_label_residual = (
    -weighted * laplace_star / r
    + sp.diff(weighted * sp.diff(stream, r) / r, r)
    - sp.diff(weighted, r) * sp.diff(stream, r) / r
    + weighted * sp.diff(stream, z, 2) / r
)
assert sp.simplify(regular_label_residual) == 0

# Chain rule identifying the source as an axial derivative of 2G/r^4.
chain_residual = (
    2 * sp.diff(G(xi), z) / r**4
    - 2 * sp.Subs(sp.diff(G(sp.Symbol("s")), sp.Symbol("s")),
                  sp.Symbol("s"), xi).doit() * sp.diff(xi, z) / r**4
)
chain_residual = chain_residual.subs(
    sp.Subs(sp.diff(G(sp.Symbol("s")), sp.Symbol("s")),
            sp.Symbol("s"), xi),
    xi * D(xi),
)
assert sp.simplify(chain_residual) == 0

print("4 exact mixed-Casimir flux checks passed")
