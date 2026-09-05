"""Independent real-phase normalization audit of the annular supplier."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0203-independent-annular-review")
    phase = s.symbols("phase", real=True)
    f, aa, bb, shear, vort = s.symbols("F A B Wprime Z", real=True)
    real = s.Matrix([f*s.cos(phase), -aa*s.sin(phase), -bb*s.sin(phase)])
    imag = s.Matrix([f*s.sin(phase), aa*s.cos(phase), bb*s.cos(phase)])
    density = s.simplify(s.Matrix([0, -shear, vort]).dot(real.cross(imag)))
    checks.check("real phase pair retains both actual vorticity components",
                 s.simplify(density-f*(shear*bb+vort*aa)) == 0)
    checks.check("full angular action measure is two pi, not pi",
                 s.simplify(s.integrate(density, (phase, 0, 2*s.pi))
                            -2*s.pi*f*(shear*bb+vort*aa)) == 0)
    r, m, k, doppler, p, angular = s.symbols("r m k sigma P Omega", positive=True)
    reduced = density.subs({aa: m*p/(r*doppler**2)-2*angular*f/doppler,
                            bb: k*p/doppler**2})
    expected = f*((m*vort/r+k*shear)*p/doppler**2
                  -2*angular*vort*f/doppler)
    checks.check("physical cross product gives the quoted complete KKS integrand",
                 s.simplify(reduced-expected) == 0)
    beta = s.symbols("beta", negative=True)
    omega, clock = s.symbols("omega c", positive=True)
    unit = s.Matrix([[0, 1], [-1, 0]])
    embedding = s.Matrix([[0, 1/(omega*clock)], [1/clock, 0]])
    mass = -beta/(omega*clock**2)
    checks.check("measured angle and physical rate inherit positive full phase",
                 s.simplify(embedding.T*(beta*unit)*embedding-mass*unit)
                 == s.zeros(2))
    checks.check("same embedding inherits the entire positive laboratory energy",
                 s.simplify(embedding.T*(-beta*omega*s.eye(2))*embedding
                            -mass*s.diag(omega**2, 1)) == s.zeros(2))
    checks.check("clock mass is strictly positive without a unit spin assumption",
                 mass.is_positive is True)
    n0, n1, n2 = s.symbols("n0 n1 n2")
    matrix = s.Matrix([[1, n, n*(n-1)] for n in (n0, n1, n2)])
    checks.check("fixed tag derivative rows have the full Vandermonde independence",
                 s.factor(matrix.det()-(n1-n0)*(n2-n0)*(n2-n1)) == 0)
    checks.check("admissible profile trial has a strict analytic negative margin",
                 s.simplify((angular/2)**2-2*angular**2)
                 == -s.Rational(7, 4)*angular**2)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
