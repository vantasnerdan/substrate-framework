#!/usr/bin/env python3
"""Materialized source-equivalent form of the original P253/0010 one-liner.

The historical output is preserved separately from the Codex session
transcript.  This file was not rerun while repairing the provenance receipt.
"""

import sympy as s


xi, e, r, z, a = s.symbols("xi e r z a", nonzero=True)
f = xi**2 / (2 * r**2) + a * xi / e - e * z * xi
xs = s.solve(s.diff(f, xi), xi)[0]
print("xi_star =", s.factor(xs))
print("reduced =", s.factor(f.subs(xi, xs)))

eta, chi = s.symbols("eta chi")
q = chi**2 / r**2 - 2 * e * eta * chi
print(
    "square_identity_residual =",
    s.expand(q - ((chi / r - e * r * eta) ** 2 - e**2 * r**2 * eta**2)),
)
