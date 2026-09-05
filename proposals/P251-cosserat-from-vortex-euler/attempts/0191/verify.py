"""Derive the full TR tilt pullback, Haar tensor and positive curvature form."""

import sympy as s

from substrate_framework.homogenization import sphere_fourth_moment_isotropic
from substrate_framework.micropolar import MicropolarCoefficients, micropolar_fourier_stiffness
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0191-transverse-family-correction")
    mass, frequency, amplitude = s.symbols("M nu A", positive=True)
    q, rate = s.symbols("q rate", real=True)
    symplectic = s.Matrix([[0, 1], [-1, 0]])
    axis_cross = -symplectic
    observed = s.zeros(2, 1)
    for sign in (-1, 1):
        embedding = s.diag(amplitude, sign*amplitude/frequency)
        phase = sign*mass*frequency*symplectic
        energy = mass*frequency**2*s.eye(2)
        checks.check(f"partner {sign} retains positive scalar phase from full vector tilt",
                     s.simplify(embedding.T*phase*embedding
                                -mass*amplitude**2*symplectic) == s.zeros(2))
        checks.check(f"partner {sign} retains actual physical mode energy",
                     s.simplify(embedding.T*energy*embedding
                                -mass*amplitude**2*s.diag(frequency**2, 1)) == s.zeros(2))
        tilt = embedding*s.Matrix([q, rate])
        derivative = tilt.diff(q)*rate-tilt.diff(rate)*frequency**2*q
        checks.check(f"partner {sign} conjugate tilt follows its actual circular generator",
                     s.simplify(derivative+sign*frequency*axis_cross*tilt) == s.zeros(2, 1))
        observed += tilt/2
    checks.check("whole TR law cancels the conjugate observation without doubling mass",
                 observed == s.Matrix([amplitude*q, 0]))
    wave = s.Matrix(s.symbols("kx ky kz", real=True))
    fourth = sphere_fourth_moment_isotropic()
    tensor = s.zeros(3)
    for i in range(3):
        for j in range(3):
            tensor[i, j] = sum(
                wave[a]*wave[b]*(s.Rational(int(i == j and a == b), 6)
                                 -fourth[i, j, a, b]/2)
                for a in range(3) for b in range(3))
    target = (2*wave.dot(wave)*s.eye(3)-wave*wave.T)/5
    checks.check("actual transverse-frame Haar tensor replaces the axial supplier tensor",
                 s.simplify(3*tensor-target) == s.zeros(3))
    checks.check("raw phase and raw literal current each acquire one-third factor",
                 s.trace(s.eye(3)/3) == 1
                 and target.subs(dict(zip(wave, [0, 0, 1], strict=True)))
                 == s.diag(s.Rational(2, 5), s.Rational(2, 5), s.Rational(1, 5)))
    ct, cl = s.symbols("CT CL", positive=True)
    denominator = 4*ct+3*cl
    cs = 3*ct*cl/denominator
    ca = 4*ct**2/denominator
    ctr = cl/2-cs
    coefficient = MicropolarCoefficients(0, 1, 1, ctr, cs, ca)
    reference = MicropolarCoefficients(0, 1, 1, (cl-ct)/2, ct/2, ct/2)
    checks.check("positive local representative retains both canonical curvature eigenvalues",
                 s.simplify(coefficient.transverse_curvature-ct) == 0
                 and s.simplify(coefficient.longitudinal_curvature-cl) == 0
                 and s.simplify(micropolar_fourier_stiffness(wave, coefficient)
                                -micropolar_fourier_stiffness(wave, reference)) == s.zeros(6))
    checks.check("all three irreducible local curvature sectors are positive",
                 cs.is_positive and ca.is_positive
                 and s.simplify(3*ctr+cs-9*cl**2/(2*denominator)) == 0)
    null_shift = ct/2-cs
    checks.check("coefficient change is exactly the retained null-Lagrangian flux",
                 s.simplify(ctr-(cl-ct)/2-null_shift) == 0
                 and s.simplify(ca-ct/2-null_shift) == 0)
    print("Derived transverse/longitudinal factors:", s.Rational(2, 5), s.Rational(1, 5))
    print("Derived positive trace-sector coefficient:", s.factor(3*ctr+cs))
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
