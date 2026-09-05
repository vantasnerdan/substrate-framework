"""Independent real-field variations of the Fourier operator and local action."""

import pytest
import sympy as s

from substrate_framework.micropolar import (
    MicropolarCoefficients,
    isotropic_micropolar_energy,
    micropolar_fourier_stiffness,
    uniform_phase_average,
)


def coefficients():
    return MicropolarCoefficients(*s.symbols("lam mu alpha ctr cs ca", real=True))


def test_fourier_operator_from_real_sine_cosine_action_variation():
    c = coefficients()
    k = s.Matrix(s.symbols("k0:3", real=True))
    u = s.Matrix(s.symbols("u0:3", real=True))
    phi = s.Matrix(s.symbols("p0:3", real=True))
    # U=u*cos(k.x), Phi=phi*sin(k.x). The two quadratures average to 1/2.
    averaged = (isotropic_micropolar_energy(-u*k.T, phi, s.zeros(3), c)
                +isotropic_micropolar_energy(s.zeros(3), s.zeros(3, 1), phi*k.T, c))/2
    actual = s.hessian(averaged, list(u)+list(phi))
    phase = s.diag(1, 1, 1, -s.I, -s.I, -s.I)
    operator = micropolar_fourier_stiffness(k, c)
    expected = phase.conjugate().T*operator*phase/2
    assert s.simplify(actual-expected) == s.zeros(6)
    assert operator.conjugate().T == operator


def test_local_angular_reaction_and_rigid_frame_covariance():
    c = coefficients()
    h = s.Matrix(3, 3, s.symbols("h0:9", real=True))
    g = s.Matrix(3, 3, s.symbols("g0:9", real=True))
    phi = s.Matrix(s.symbols("p0:3", real=True))
    w = isotropic_micropolar_energy(h, phi, g, c)
    stress = s.Matrix(3, 3, lambda i, j: s.diff(w, h[i, j]))
    axial = s.Matrix([sum(s.LeviCivita(i, j, k)*stress[j, k]
                         for j in range(3) for k in range(3)) for i in range(3)])
    assert s.simplify(s.Matrix([s.diff(w, p) for p in phi])-axial) == s.zeros(3, 1)
    angle = s.Matrix(s.symbols("a0:3", real=True))
    rotation = s.Matrix([[0, -angle[2], angle[1]],
                         [angle[2], 0, -angle[0]], [-angle[1], angle[0], 0]])
    transformed = isotropic_micropolar_energy(h+rotation, phi+angle, g, c)
    assert s.expand(transformed-w) == 0


def test_both_transverse_helicities_and_longitudinal_spin():
    c = coefficients()
    k = s.symbols("k", positive=True)
    operator = micropolar_fourier_stiffness([0, 0, k], c)
    for helicity in (-1, 1):
        axis = s.Matrix([1, s.I*helicity, 0])/s.sqrt(2)
        columns = s.zeros(6, 2)
        columns[:3, 0], columns[3:, 1] = axis, axis
        projected = s.simplify(columns.conjugate().T*operator*columns)
        target = s.Matrix([[(c.shear+c.locking)*k*k, -2*c.locking*helicity*k],
                            [-2*c.locking*helicity*k,
                             c.transverse_curvature*k*k+4*c.locking]])
        assert s.simplify(projected-target) == s.zeros(2)
    assert s.expand(operator[5, 5]-4*c.locking
                    -2*(c.symmetric_curvature+c.trace_curvature)*k*k) == 0
    # Historical subtraction of grad-div rather than addition is exposed.
    wrong = 2*(c.skew_curvature-c.trace_curvature)*k*k+4*c.locking
    assert s.simplify(operator[5, 5]-wrong) != 0
    assert operator[:3, 5] == s.zeros(3, 1)


def test_structure_free_action_and_argument_shapes():
    k = s.symbols("k", real=True)
    fluid = MicropolarCoefficients(0, 0, 0, 0, 0, 0)
    assert micropolar_fourier_stiffness([0, 0, k], fluid) == s.zeros(6)
    with pytest.raises(ValueError, match="shape"):
        micropolar_fourier_stiffness([1, 2], fluid)


def test_product_phase_closure_retains_energy_but_not_conservative_stiffness():
    p, r, q1, q2 = s.symbols("p r q1 q2", real=True)
    stiffness = s.Symbol("K", positive=True)
    energy = stiffness*(1-s.cos(p-r+q1-q2))
    average = uniform_phase_average(energy, (p, r))
    assert average == stiffness
    assert s.hessian(average, [q1, q2]) == s.zeros(2)
    # Identical uniform one-point marginals, different joint law: locking survives.
    correlated = uniform_phase_average(energy.subs(r, p), (p,))
    assert s.diff(correlated, q1, 2).subs({q1: 0, q2: 0}) == stiffness
    assert s.simplify(correlated-stiffness*(1-s.cos(q1-q2))) == 0
    with pytest.raises(ValueError, match="distinct"):
        uniform_phase_average(energy, (p, p))
