#!/usr/bin/env python3
"""Short exact checks for P253/0029; no floating-point or production numerics."""

from fractions import Fraction as Q
import sympy as sp


def d(j: int) -> Q:
    assert abs(j) >= 2
    return -Q(j, 4 * (j * j - 1) * abs(j))


checks = []

# Exact d_j parity and sharp high-mode bounds.
checks.append(all(d(-j) == -d(j) for j in range(2, 100)))
checks.append(max(abs(d(j)) for j in range(-99, 100) if abs(j) >= 2) == Q(1, 12))
checks.append(max(abs(d(j) / j) for j in range(-99, 100) if abs(j) >= 2) == Q(1, 24))

# Worst-case autonomous gaps under epsilon^2 |c3| <= 1/2.
eta_gap = Q(1, 2)
same_gap = Q(1, 3) - eta_gap * Q(1, 6)
opposite_gap = 4 * Q(1, 3) - 2 * Q(13, 24) - eta_gap * Q(1, 6)
checks.extend([same_gap == Q(1, 4), opposite_gap == Q(1, 6)])

# Positive principal metric under epsilon^2 |c3| <= 3/4.
eta_metric = Q(3, 4)
metric_lower = Q(1, 3) - Q(13, 48) - eta_metric * Q(1, 24)
metric_upper = Q(2, 3) + eta_metric * Q(1, 24)
checks.extend([metric_lower == Q(1, 32), metric_upper == Q(67, 96)])
checks.append(metric_upper / metric_lower == Q(67, 3))

# Homological sign in e^{-X} L e^X: R + i Delta X = 0 for X=iR/Delta.
R, Delta = sp.symbols("R Delta", nonzero=True, real=True)
X = sp.I * R / Delta
checks.append(sp.simplify(R + sp.I * Delta * X) == 0)

# Resonant derivative elimination: alpha' ell is replaced by -beta times
# the spatial frequency difference.
alpha, ell, c1, c2, c3, n, s, eps, e = sp.symbols(
    "alpha ell c1 c2 c3 n s eps e", nonzero=True, real=True
)
beta, c1p, c2p, c3p = sp.symbols("beta c1p c2p c3p", real=True)
raw = beta * alpha * ell + c1p * n + c2p * s + eps**2 * c3p * e
resonance_sub = {alpha * ell: -(c1 * n + c2 * s + eps**2 * c3 * e)}
reduced = n * (c1p - beta * c1) + s * (c2p - beta * c2) + eps**2 * e * (c3p - beta * c3)
checks.append(sp.simplify(raw.subs(resonance_sub) - reduced) == 0)

# A=J L is symplectic, and a positive symmetric L is its invariant metric.
a, b, c = sp.symbols("a b c", real=True)
J = sp.Matrix([[0, 1], [-1, 0]])
Lmat = sp.Matrix([[a, b], [b, c]])
A = J * Lmat
checks.append(sp.simplify(A.T * J + J * A) == sp.zeros(2))
checks.append(sp.simplify(A.T * Lmat + Lmat * A) == sp.zeros(2))

# Exact contour KKS/Hessian normalization gives the principal generator.
jj, ee, mm = sp.symbols("jj ee mm", positive=True, real=True)
W = -(ee**4 / jj) * J
LH = (ee**2 * mm / jj) * sp.eye(2)
A_from_kks = -W.inv() * LH
checks.append(sp.simplify(A_from_kks + ee**-2 * mm * J) == sp.zeros(2))

# A common positive physical factor cancels from i_X Omega=dH.
k, Om, H = sp.symbols("k Om H", nonzero=True)
checks.append(sp.simplify((k * H) / (k * Om) - H / Om) == 0)

assert all(checks)
print(f"exact checks passed: {len(checks)}")
print(f"autonomous gaps: same={same_gap}, opposite={opposite_gap}")
print(f"principal metric bounds: [{metric_lower}, {metric_upper}], ratio={metric_upper/metric_lower}")
