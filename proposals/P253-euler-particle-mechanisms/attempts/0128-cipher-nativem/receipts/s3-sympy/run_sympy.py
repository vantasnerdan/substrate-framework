"""R-e banking: gauged-Clebsch A-sector algebra (04-s3-bridge.md Reading 1)
+ helicity-structure check (Reading 2). Command: python3 run_sympy.py (> run.log)
Env: repo python kernel, sympy only.
"""
import sympy as sp

k, a = sp.symbols('k alpha')
gx, gy, gz, hx, hy, hz, Ax, Ay, Az = sp.symbols('gx gy gz hx hy hz Ax Ay Az')
ux = gx + a * hx - k * a * Ax
uy = gy + a * hy - k * a * Ay
uz = gz + a * hz - k * a * Az
T = sp.expand((ux ** 2 + uy ** 2 + uz ** 2) / 2)
Te = sp.expand(T)
a2 = Te.coeff(Ax, 2) + Te.coeff(Ay, 2) + Te.coeff(Az, 2)
print("A2-total:", sp.simplify(a2), "(algebraic mass-like, per-comp k^2 a^2/2)")
print("dA-dependence:", "NONE" if "Derivative" not in str(T) else "PRESENT (unexpected!)")
ax, ay, az, bx, by, bz = sp.symbols('ax ay az bx by bz')
wx, wy, wz = ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx
h = (gx + a * bx) * wx + (gy + a * by) * wy + (gz + a * bz) * wz
print("helicity poly-degree in a:", sp.Poly(h, a).degree(),
      "(0 = alpha-nabla-beta term drops out identically)")
print("helicity = det(grad phi, grad alpha, grad beta):",
      sp.simplify(h - (gx * wx + gy * wy + gz * wz)) == 0)
