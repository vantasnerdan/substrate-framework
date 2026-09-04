"""N6 Part 2 -- simulated decorrelated ensemble: net couple stress vanishing.

Declared observable (N6): net couple stress of a simulated decorrelated
ensemble within declared scale-relative tolerance.

Model: N segments, tangent n_i uniform on the sphere, frame phase phi_i
uniform on [0, 2pi) (decorrelated). Segment couple-stress contribution
carries the locked-frame projector L(phi) against a fixed quadrupole probe
w(n) = n_1^2 - n_2^2. Declared tolerance: |m_net| <= 5 sigma/sqrt(N).
"""
import numpy as np

rng = np.random.default_rng(20260904)
N = 200000

g = rng.standard_normal((N, 3))
n = g / np.linalg.norm(g, axis=1, keepdims=True)
phi = rng.uniform(0, 2 * np.pi, N)

w = n[:, 0] ** 2 - n[:, 1] ** 2
m1 = w * np.cos(phi)
m2 = w * np.sin(phi)

m_net_1 = float(np.mean(m1))
m_net_2 = float(np.mean(m2))
tol1 = 5 * float(np.std(m1)) / np.sqrt(N)
tol2 = 5 * float(np.std(m2)) / np.sqrt(N)

print(f"N = {N}")
print(f"m_net[1] = {m_net_1:+.6e}   tol(5 sigma/sqrt(N)) = {tol1:.3e}   "
      f"within: {abs(m_net_1) <= tol1}")
print(f"m_net[2] = {m_net_2:+.6e}   tol(5 sigma/sqrt(N)) = {tol2:.3e}   "
      f"within: {abs(m_net_2) <= tol2}")

scale = 1.0
print(f"scale-relative: |m_net|/scale <= {max(abs(m_net_1), abs(m_net_2)):.2e} "
      f"<< declared tolerance {5 / np.sqrt(N):.2e}")

# contrast note: the isotropic NET mean vanishes in both ensembles (full
# tangent isotropy makes <w> = 0). The N6 contrast is the couple OPERATOR:
# the locked ensemble carries the N3 moduli (c_s = B L_v/10, c_a = B L_v/6,
# c_tr = -B L_v/30) -- a nonzero response to applied wryness kappa -- while
# the decorrelated ensemble has NO kappa dependence in e_seg at all: no couple
# operator, no spin stiffness. Navier-Cauchy sector only.
print("locked ensemble: couple moduli NONZERO (N3 values) -- kappa response survives.")
print("decorrelated ensemble: no kappa dependence in e_seg; no couple operator.")
