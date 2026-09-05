"""Exact pressure, scale and physical-mean identities for the same array."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0141-same-array-pressure-jets")
    x, y = s.symbols("x y", real=True)
    psi = s.Function("psi")(x, y)
    bernoulli = s.Function("B")(psi)
    energy = s.Symbol("C", positive=True)
    vertical = s.sqrt(2*(energy-bernoulli))
    planar = s.Matrix([-s.diff(psi, y), s.diff(psi, x)])
    checks.check("Bernoulli lift axial speed is transported by the same planar flow",
                 s.simplify(planar.dot(s.Matrix([s.diff(vertical, x), s.diff(vertical, y)]))) == 0)
    horizontal_curl = s.Matrix([s.diff(vertical, y), -s.diff(vertical, x)])
    factor = s.diff(bernoulli, psi)/vertical
    checks.check("same-field Bernoulli lift gives the correct horizontal force-free curl",
                 (horizontal_curl-factor*planar).applyfunc(s.simplify) == s.zeros(2, 1))

    p1, p2, a1, a2, k = s.symbols("p1 p2 a1 a2 k", real=True)
    denominator = (p1+s.I*a1)**2+(p2+s.I*a2)**2+k*k
    real_part = s.re(s.expand(denominator))
    imag_part = s.im(s.expand(denominator))
    checks.check("weighted complete pressure denominator has the massive positive real part",
                 s.expand(real_part-(p1*p1+p2*p2+k*k-a1*a1-a2*a2)) == 0)
    checks.check("weighted pressure modulus retains the imaginary shift term",
                 s.expand(denominator*s.conjugate(denominator)-real_part**2-imag_part**2) == 0)
    checks.check("positive axial carrier is essential to the stated weight domain",
                 real_part.subs({p1: 0, p2: 0, k: 0}) == -a1*a1-a2*a2)

    q1, q2, r1, r2 = s.symbols("Q1 Q2 R1 R2", real=True)
    bloch_phase = s.exp(s.I*(q1*r1+q2*r2))
    for index, coordinate in enumerate((q1, q2)):
        vector_component = (r1, r2)[index]
        checks.check(f"actual intercell Bloch derivative {index} keeps its displacement moment",
                     s.simplify(s.diff(bloch_phase, coordinate)-s.I*vector_component*bloch_phase) == 0)
    checks.check("mixed transverse Bloch jet contains the full cell-displacement product",
                 s.simplify(s.diff(bloch_phase, q1, q2)+r1*r2*bloch_phase) == 0)
    physical_mean_weight = s.exp(-s.I*(q1*x+q2*y))
    checks.check("physical mean has an algebraic cell moment independent of hopping",
                 s.simplify(s.diff(physical_mean_weight, q1).subs({q1: 0, q2: 0})+s.I*x) == 0)

    # Exact thin-shell radial mass-border limit: Y=1 inside; the shell
    # equation (rY')'=-rQY carries integrated strength 2.
    radius = s.Symbol("a", positive=True)
    constant, log_coefficient = s.symbols("A B", real=True)
    radial = s.Symbol("r", positive=True)
    exterior = constant+log_coefficient*s.log(radial)
    solved = s.solve([exterior.subs(radial, radius)-1,
                      (radial*s.diff(exterior, radial)).subs(radial, radius)+2],
                     (constant, log_coefficient))
    checks.check("thin Rankine radial mass border has nonzero logarithmic response",
                 solved[log_coefficient] == -2)

    step, carrier = s.symbols("K k0", real=True)
    coefficients = s.symbols("e0:5", real=True)
    perturbation = sum(coefficients[n]*(k-carrier)**n for n in range(5))
    central = s.expand(perturbation.subs(k, carrier+step)
                       +perturbation.subs(k, carrier-step)-2*perturbation.subs(k, carrier))
    checks.check("symmetric physical carrier difference cancels onsite and first-jet errors",
                 central == 2*coefficients[2]*step**2+2*coefficients[4]*step**4)
    checks.check("its natural gradient error is controlled by the actual second carrier jet",
                 s.diff(central, step, 2).subs(step, 0)
                 == 2*s.diff(perturbation, k, 2).subs(k, carrier))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
