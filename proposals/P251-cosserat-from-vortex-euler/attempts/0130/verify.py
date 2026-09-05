"""Quasiperiodic Euler lift and nonempty physical response window."""

import sympy as s

from substrate_framework.euler_fourier import ZERO, leray
from substrate_framework.euler_observation import material_tag_fourier_dipole
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0130-finite-window-join")
    p = s.Matrix(s.symbols("p0:3", real=True))
    projector = s.eye(3)-p*p.T/p.dot(p)
    checks.check("full physical pressure projector is transverse and idempotent",
                 s.simplify(p.T*projector) == s.zeros(1, 3)
                 and s.simplify(projector**2-projector) == s.zeros(3))
    checks.check("pressure force has zero physical curl",
                 s.simplify(p.cross((s.eye(3)-projector)*s.Matrix([1, 2, 3])))
                 == s.zeros(3, 1))
    # A near resonance does not create a large Leray multiplier.
    tiny = s.Symbol("epsilon", positive=True)
    ray = s.Matrix([2, -3, 1])*tiny
    p_ray = s.simplify(projector.subs(dict(zip(p, ray)), simultaneous=True))
    checks.check("Leray norm is unchanged as a nonzero physical mode approaches zero",
                 p_ray.diff(tiny) == s.zeros(3) and p_ray.eigenvals() == {0: 1, 1: 2})
    checks.check("zero physical mode retains its actual harmonic velocity",
                 leray(({ZERO: 2}, {ZERO: -3}, {ZERO: 1}))
                 == ({ZERO: 2}, {ZERO: -3}, {ZERO: 1}))

    matrix = s.Matrix([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
    k = s.Matrix(s.symbols("k0:3", real=True))
    lift = matrix.col_join(k.T)
    label = s.Matrix([2, -1, 3, 5])
    physical_wave = lift.T*label
    checks.check("lifted Fourier differentiation equals the physical chain rule",
                 physical_wave == matrix.T*label[:3, 0]+5*k)
    coefficients = s.Matrix(s.symbols("u0:3", real=True))
    checks.check("four-dimensional transport divergence is physical incompressibility",
                 s.expand(label.dot(lift*coefficients)-physical_wave.dot(coefficients)) == 0)
    # Nonzero four-label can represent a constant PHYSICAL mode.
    resonant_k = -matrix.T*label[:3, 0]/5
    checks.check("resonant nonzero labels are not discarded as nonconstant modes",
                 physical_wave.subs(dict(zip(k, resonant_k)), simultaneous=True) == s.zeros(3, 1)
                 and label != s.zeros(4, 1))

    density, inertia, frequency, norm_k, c_k = s.symbols("rho j Omega k C_k", positive=True)
    axis = s.Matrix([1, 0, 0])
    wave = s.Matrix([0, 0, norm_k])
    spin_change = -inertia*frequency*axis
    delta_hybrid = -material_tag_fourier_dipole(wave, spin_change, s.zeros(3))
    checks.check("the ideal measured spin induces the transverse hybrid velocity",
                 s.simplify(delta_hybrid/density-s.I*inertia*frequency*wave.cross(axis)/(2*density))
                 == s.zeros(3, 1))
    relative_lower = s.Rational(1, 2)-s.Rational(1, 16)-s.Rational(1, 16)
    checks.check("linear response window leaves a strict physical signal margin",
                 relative_lower == s.Rational(3, 8) and relative_lower > 0)
    nonlinear_lower = relative_lower-s.Rational(1, 16)
    checks.check("small genuine nonlinear flows retain a nonzero normalized signal",
                 nonlinear_lower == s.Rational(5, 16) and nonlinear_lower > s.Rational(1, 4))
    bound = norm_k*inertia*frequency*(s.Rational(1, 2)-s.Rational(1, 16))-c_k*norm_k**2
    edge = inertia*frequency/(16*c_k)
    checks.check("declared wavenumber endpoint satisfies the unsimplified error inequality",
                 s.simplify(bound.subs(norm_k, edge)-relative_lower*edge*inertia*frequency) == 0)
    checks.check("deleting the actual spin deletes this coupling mechanism",
                 delta_hybrid.subs(inertia, 0) == s.zeros(3, 1))
    print("Scope: actual finite-window first-jet coupling, not autonomous second-gradient closure")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()
