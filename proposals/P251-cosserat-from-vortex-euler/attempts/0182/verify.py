"""Derive nonunit-overlap current normalization from both actual forms."""

import sympy as s

from substrate_framework.micropolar import MicropolarCoefficients, micropolar_fourier_stiffness
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0182-nonunit-current-map")
    k, rho, j, eta, inertia, ma, j2, a, gap, bt, bl = s.symbols(
        "k rho j eta I ma j2 a gap BT BL", positive=True)

    def jet(matrix):
        return matrix.applyfunc(lambda value: s.series(value, k, 0, 3).removeO().expand())

    physical = s.Matrix([[1-inertia*k**2/(4*rho), -eta*j*k/(2*rho)], [k/2, 1]])
    improvement = s.Matrix([[1, (eta-1)*j*k/(2*rho)], [0, 1]])
    effective_i = inertia+(1-eta)*j
    composed = jet(improvement*physical)
    target = s.Matrix([[1-effective_i*k**2/(4*rho), -j*k/(2*rho)], [k/2, 1]])
    checks.check("actual current improvement preserves the angle and complete map",
                 s.simplify(composed-target) == s.zeros(2))
    mass0 = s.diag(rho+ma*k**2, j+j2*k**2)
    stiff0 = s.diag(rho*a*k**2, j*gap+(j*bt+gap*j2)*k**2)
    inverse_physical = jet(physical.inv())
    mass_physical = jet(inverse_physical.T*mass0*inverse_physical)
    checks.check("nonunit physical overlap has nonzero mixed mass",
                 s.simplify(mass_physical[0, 1]-(eta-1)*j*k/2) == 0)
    inverse = jet(composed.inv())
    mass = jet(inverse.T*mass0*inverse)
    stiff = jet(inverse.T*stiff0*inverse)
    expected_mass = s.diag(rho+(ma+effective_i/2-j/4)*k**2,
                          j+(j2-j**2/(4*rho))*k**2)
    checks.check("both actual gradient masses retain the overlap/current change",
                 s.simplify(mass-expected_mass) == s.zeros(2))
    normalize = s.eye(2)-s.diag(1/rho, 1/j)*(mass-s.diag(rho, j))/2
    normal_mass = jet(normalize.T*mass*normalize)
    normal_stiff = jet(normalize.T*stiff*normalize)
    checks.check("one derivative transformation normalizes both inertias",
                 s.simplify(normal_mass-s.diag(rho, j)) == s.zeros(2))
    ct = j*bt-j**2*gap/(4*rho)
    coefficients = MicropolarCoefficients(0, rho*a, j*gap/4, (j*bl-ct)/2, ct/2, ct/2)
    canonical = s.Matrix(micropolar_fourier_stiffness([0, 0, k], coefficients))
    helicity = s.Matrix([1, s.I, 0])/s.sqrt(2)
    embedding = s.zeros(6, 2)
    embedding[:3, 0], embedding[3:, 1] = helicity, helicity
    checks.check("full normalized potential agrees with the canonical Cosserat pencil",
                 s.simplify(normal_stiff-embedding.conjugate().T*canonical*embedding) == s.zeros(2))
    checks.check("overlap and tagged-gradient masses are retained only in physical reconstruction",
                 not normal_stiff.has(eta, inertia, ma, j2))
    optical = physical*s.Matrix([0, 1])
    reconstructed = improvement.inv()*composed*s.Matrix([0, 1])
    checks.check("physical optical transfer survives exactly when overlap is positive",
                 s.simplify(optical[0]+eta*j*k/(2*rho)) == 0
                 and s.simplify(reconstructed[0]-optical[0]) == 0
                 and optical[0].subs(eta, 0) == 0)
    checks.check("uniform rigid-frame response is unchanged by the curl improvement",
                 improvement.subs(k, 0) == s.eye(2))
    print("Derived C_T:", ct)
    print("Literal physical optical transfer:", optical[0])
    print("Stationary mode/current existence remains the actual0181 obligation.")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
