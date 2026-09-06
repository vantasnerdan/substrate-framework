"""Exact distributional exposing test; not a fluid simulation."""

import sympy as sp

z = sp.Symbol("z", real=True)
kappa = sp.Symbol("kappa", positive=True)
psi = sp.Heaviside(-z)
residual = sp.diff(psi, z, 2)*kappa**2/2
phi = z*sp.exp(-z*z)
pairing = sp.integrate(residual*phi, (z, -sp.oo, sp.oo))
restricted_pairing = sp.integrate(residual*sp.exp(-z*z), (z, -sp.oo, sp.oo))
assert residual == -kappa**2*sp.DiracDelta(z, 1)/2
assert pairing == kappa**2/2
assert restricted_pairing == 0
print("Schrodinger residual:", residual)
print("Unrestricted exposing test:", pairing)
print("Zero-normal-derivative test:", restricted_pairing)
print("PASS: full distribution residual is nonzero; restricted tests miss it.")
