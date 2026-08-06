"""Independent recomputation of C-SIG-001 S-value. No repo imports."""
import math

h = 1.0
xi = 0.0
S = 0.0
species = []
# kink-antikink: mass 8, multiplicity 2 (kink + antikink)
for name, m, mult in [("kink", 8.0, 2)]:
    contrib = mult * (1/6 - xi) * m**2 * math.log((8*math.pi/m)**2) / (2*math.pi)
    species.append((name, m, mult, contrib))
    S += contrib
# breather tower E_n = 16 sin(n h/16), n=1..25, n*h < 8*pi
for n in range(1, 26):
    E = 16.0 * math.sin(n * h / 16.0)
    assert n * h < 8 * math.pi
    contrib = (1/6 - xi) * E**2 * math.log((8*math.pi/E)**2) / (2*math.pi)
    species.append((f"breather n={n}", E, 1, contrib))
    S += contrib

q = 1.0 / S
print(f"S = {S:.12f}")
print(f"q = {q:.14f}")
print(f"expected S = 120.058077992, diff = {S - 120.058077992:.3e}")
print(f"expected q = 0.00832930209048, diff = {q - 0.00832930209048:.3e}")
for s in species[:5]:
    print(s)
