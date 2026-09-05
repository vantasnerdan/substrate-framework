"""Exact whole-field parity, axial Haar and inherited optical action checks."""

import sympy as s

from substrate_framework.euler_displacement_preparation import transverse_pair_average
from substrate_framework.euler_optical_response import (
    OpticalModeJet,
    correlated_optical_preparation,
)
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0202-actual-axial-parity-interface")
    direction = s.symbols("t0:3", real=True)
    polarization = s.symbols("d0:3", real=True)
    t, d = s.Matrix(direction), s.Matrix(polarization)
    k = s.Matrix(s.symbols("K0:3", real=True))
    tensor = s.Matrix(3, 3, lambda i, j: transverse_pair_average(
        3*t[i]*t[j]*t.dot(k)**2*d.dot(d), direction, polarization))
    expected = (k.dot(k)*s.eye(3)+2*k*k.T)/5
    checks.check("actual axial direction gives the full fourth-moment tensor",
                 s.simplify(tensor-expected) == s.zeros(3))
    checks.check("axial longitudinal coefficient is three fifths",
                 s.simplify(tensor*k-s.Rational(3, 5)*k.dot(k)*k) == s.zeros(3, 1))
    transverse = k.cross(s.Matrix([1, 2, 3]))
    checks.check("axial transverse coefficient is one fifth",
                 s.simplify(tensor*transverse-k.dot(k)*transverse/5) == s.zeros(3, 1))
    mirror = s.diag(-1, 1, 1)
    carrier = s.Matrix([0, 0, 1])
    axis = mirror.det()*mirror*carrier
    phi = s.Matrix(s.symbols("Phi0:3", real=True))
    amp, delta, companion = s.symbols("A Delta Lambda", real=True)
    theta = amp*carrier.dot(phi)
    theta_mirror = amp*axis.dot(phi)
    checks.check("physical axial current adds under positive mirror pairing",
                 s.simplify((carrier*delta*theta+axis*delta*theta_mirror)/2
                            -carrier*delta*theta) == s.zeros(3, 1))
    checks.check("literal polar centroid companion cancels under the same pairing",
                 s.simplify(carrier*companion*theta
                            +mirror*carrier*companion*theta_mirror) == s.zeros(3, 1))
    checks.check("same carrier increment preserves cancellation of its derivatives",
                 (mirror*carrier).dot(k) == carrier.dot(k))
    checks.check("canceled polar mean retains its quadratic companion norm",
                 s.simplify((companion**2*theta**2+companion**2*theta_mirror**2)/2
                            -companion**2*theta**2) == 0)
    mass, nu = s.symbols("M nu", positive=True)
    unit = s.Matrix([[0, 1], [-1, 0]])
    reverse = s.diag(1, -1)
    phase = mass*unit
    energy = mass*s.diag(nu**2, 1)
    checks.check("physical rate reversal preserves rather than doubles phase",
                 (phase+reverse.T*(-phase)*reverse)/2 == phase)
    checks.check("the identical reversal preserves the full positive energy",
                 (energy+reverse.T*energy*reverse)/2 == energy)
    modes = (OpticalModeJet(2, 1, s.Rational(1, 3), 3, s.Rational(1, 7), 2),
             OpticalModeJet(2, 3, -1, 5, -2, s.Rational(3, 2)))
    preparation = correlated_optical_preparation(modes, (s.Rational(2, 5),
                                                         s.Rational(3, 5)), 7)
    checks.check("unchanged canonical preparation solves all four full-jet rows",
                 all(value == 0 for value in (
                     preparation.normalization_residual, preparation.variance_residual,
                     preparation.curvature_residual, preparation.energy_residual)))
    checks.check("positive probabilities retain positive inherited optical density",
                 preparation.J0.is_positive is True)
    j0, j2, curvature = s.symbols("J0 J2 B", positive=True)
    full_mass = j0*s.eye(3)+j2*tensor/2
    stiffness = nu**2*full_mass+j0*nu*curvature*tensor
    inverse_second_jet = s.eye(3)/j0-j2*tensor/(2*j0**2)
    error = s.expand(inverse_second_jet*stiffness
                     -nu**2*s.eye(3)-nu*curvature*tensor)
    epsilon = s.symbols("epsilon", real=True)
    scaled = error.subs(dict(zip(k, epsilon*k)))
    checks.check("gradient mass is retained in the second-jet action quotient",
                 all(s.diff(value, epsilon, order).subs(epsilon, 0) == 0
                     for value in scaled for order in range(4)))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
