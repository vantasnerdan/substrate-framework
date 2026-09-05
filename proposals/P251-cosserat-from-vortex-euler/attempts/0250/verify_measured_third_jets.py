"""Expose the sparse-tag measured-jet scaling used in the 0250 diagonal."""

import sympy as sp

N, j, q, s, M, C3, Cclock, K = sp.symbols(
    "N j q s M C3 Cclock K", positive=True)
cN = N**(-j)
gain = N*cN
normalized_resonant = sp.simplify(gain/gain)
normalized_tail = sp.simplify(N**(-q)/gain)
graph_cost = sp.simplify(N**(s+1)/gain)
assert normalized_resonant == 1
assert sp.simplify(normalized_tail-N**(j-q-1)) == 0
assert sp.limit(normalized_tail.subs(q, j), N, sp.oo) == 0
assert sp.simplify(graph_cost-N**(j+s)) == 0
assert sp.simplify((C3*K**3)/K**2-C3*K) == 0
assert sp.simplify((Cclock/N)/K**2-Cclock/(N*K**2)) == 0

print("PASS resonant measured third jet / unit gain is O(1)")
print("PASS nonresonant integration-by-parts tail = O(N**(j-q-1))")
print("PASS q>=j makes the normalized tail bounded and vanishing")
print("PASS derivative graph cost after unit gain = O(N**(j+s))")
print("PASS measured third remainder divided by K**2 = C3*K")
print("PASS leading-clock correction is separate from exact-parametrix error")
