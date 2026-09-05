"""Independent full-action variations, including the ensemble-order mutation."""

import pytest
import sympy as s

from substrate_framework.euler_orbit import (
    affine_cage_rotation_map,
    isotropic_axis_gradient,
    micropolar_kinetic_normal_form,
    reduce_euler_rotor_block,
)


def test_full_mixed_hessian_reduction_by_direct_variation():
    h = s.Matrix([[7, 2, 1], [2, 5, -1], [1, -1, 4]])
    cb, ci = s.symbols("cb ci", nonzero=True, real=True)
    r, q, p, bd, qd = s.symbols("r q p bd qd", real=True)
    coordinates = s.Matrix([r, q, p])
    action = cb*r*bd+ci*p*qd-(coordinates.T*h*coordinates)[0]/2
    momenta = s.solve([s.diff(action, r), s.diff(action, p)], [r, p])
    result = reduce_euler_rotor_block(h, cb, ci)
    velocity = s.Matrix([bd, qd])
    reconstructed = ((velocity.T*result.kinetic*velocity)[0]/2
                     +q*result.gyro.dot(velocity)-result.stiffness*q*q/2)
    assert s.simplify(action.subs(momenta)-reconstructed) == 0
    assert result.gyro[0] != 0


def test_time_reversed_ensemble_retains_independent_reactions():
    h = s.Matrix([[7, 2, 1], [2, 5, -1], [1, -1, 4]])
    q, bd, qd = s.symbols("q bd qd", real=True)
    rp, sp, rm, sm = s.symbols("rp sp rm sm", real=True)
    plus, minus = s.Matrix([rp, q, sp]), s.Matrix([rm, q, sm])
    action = ((2*rp*bd+3*sp*qd-(plus.T*h*plus)[0]/2)
              +(-2*rm*bd-3*sm*qd-(minus.T*h*minus)[0]/2))/2
    momenta = s.solve([s.diff(action, p) for p in (rp, sp, rm, sm)], (rp, sp, rm, sm))
    result = reduce_euler_rotor_block(h, 2, 3)
    velocity = s.Matrix([bd, qd])
    target = (velocity.T*result.kinetic*velocity)[0]/2-result.stiffness*q*q/2
    assert s.simplify(action.subs(momenta)-target) == 0
    reverse = reduce_euler_rotor_block(h, -2, -3)
    assert reverse.kinetic == result.kinetic
    assert reverse.stiffness == result.stiffness
    assert reverse.gyro == -result.gyro
    # The wrong tied-reaction ensemble loses all time-kinetic dependence.
    mutation = action.subs({rm: rp, sm: sp})
    assert s.diff(mutation, bd) == s.diff(mutation, qd) == 0


def test_general_physical_cage_map_and_frame_covariance():
    m = s.Matrix([[3, s.Rational(1, 2)], [s.Rational(1, 2), 2]])
    result = affine_cage_rotation_map(m, 5)
    psid, betad = s.symbols("psid betad", real=True)
    qd = (psid-betad)/result.factor
    velocity = s.Matrix([betad+qd, qd])
    actual = (velocity.T*m*velocity)[0]/2
    target = (result.spin_inertia*psid**2+result.cage_inertia*betad**2)/2
    assert s.simplify(actual-target) == 0
    assert result.spin_inertia > 0 and result.cage_inertia > 0
    assert result.stiffness == 5/result.factor**2


def test_singular_cage_and_invalid_orbit_data_are_exposed():
    with pytest.raises(ValueError, match="singular"):
        affine_cage_rotation_map([[1, -1], [-1, 2]], 1)
    with pytest.raises(ValueError, match="positive definite"):
        reduce_euler_rotor_block(s.diag(1, -1, 1), 1, 1)
    with pytest.raises(ValueError, match="nonzero"):
        reduce_euler_rotor_block(s.eye(3), 0, 1)
    with pytest.raises(ValueError, match="symmetric"):
        reduce_euler_rotor_block([[1, 1, 0], [0, 1, 0], [0, 0, 1]], 1, 1)
    with pytest.raises(ValueError, match="finite"):
        reduce_euler_rotor_block(s.diag(s.oo, 1, 1), 1, 1)
    with pytest.raises(ValueError, match="finite"):
        reduce_euler_rotor_block(s.eye(3), s.nan, 1)


def test_isotropic_tensor_from_independent_spherical_moments():
    cp, cl, density = s.symbols("cp cl nu", positive=True)
    result = isotropic_axis_gradient(s.diag(cp, cp, cl), [0, 0, 1], density)
    # C_rot=cp*I+(cl-cp)*n*n.T; exact uniform-sphere moments, not API formulas.
    def moment(indices):
        powers = [indices.count(i) for i in range(3)]
        if any(power % 2 for power in powers):
            return s.S.Zero
        numerator = s.prod(s.factorial2(power-1) for power in powers)
        return numerator/s.factorial2(sum(powers)+1)

    for i in range(3):
        for k in range(3):
            for j in range(3):
                for ell in range(3):
                    actual = density*(cp*int(j == ell)*moment([i, k])
                                      +(cl-cp)*moment([i, k, j, ell]))
                    target = (result.norm_coefficient*int(i == k)*int(j == ell)
                              +result.mixed_coefficient*(int(i == j)*int(k == ell)
                                                       +int(i == ell)*int(k == j)))
                    assert s.simplify(actual-target) == 0
    # Independently rotating the axis and C is the wrong physical ensemble.
    wrong = density*(2*cp+cl)/9
    assert s.simplify(result.norm_coefficient-wrong) != 0


def test_curvature_coercivity_does_not_require_positive_trace_modulus():
    result = isotropic_axis_gradient(s.diag(100, 100, 1), [0, 0, 1])
    assert result.trace_modulus < 0
    assert result.symmetric_modulus > 0 and result.skew_modulus > 0
    assert 3*result.trace_modulus+result.symmetric_modulus > 0
    g = s.Matrix(3, 3, s.symbols("g0:9"))
    sym, skew = (g+g.T)/2, (g-g.T)/2
    direct = (result.norm_coefficient*sum(value**2 for value in g)
              +result.mixed_coefficient*(s.trace(g)**2+s.trace(g*g)))/2
    decomposed = (result.trace_modulus*s.trace(g)**2
                  +result.symmetric_modulus*sum(value**2 for value in sym)
                  +result.skew_modulus*sum(value**2 for value in skew))
    assert s.expand(direct-decomposed) == 0
    with pytest.raises(ValueError, match="unit"):
        isotropic_axis_gradient(s.eye(3), [0, 0, 2])
    with pytest.raises(ValueError, match="positive"):
        isotropic_axis_gradient(s.eye(3), [0, 0, 1], 0)


@pytest.mark.parametrize("helicity", [-1, 1])
def test_same_field_map_normalizes_mass_and_keeps_potential_correction(helicity):
    rho, j, alpha = s.symbols("rho j alpha", positive=True)
    c, mu, mp, b, k, a = s.symbols("C mu mp b k A", real=True)
    result = micropolar_kinetic_normal_form(rho, j, alpha, c, mu, mp, b, k, helicity)
    mass = s.Matrix([[rho+mu*k*k, b*helicity*k], [b*helicity*k, j+mp*k*k]])
    potential = s.Matrix([[a*k*k, -2*alpha*helicity*k],
                          [-2*alpha*helicity*k, 4*alpha+c*k*k]])
    transform = result.field_map
    pulled_mass = transform.T*mass*transform
    pulled_potential = transform.T*potential*transform
    expected_mass = s.diag(rho, j)
    expected_potential = s.Matrix([[a*k*k, -2*alpha*helicity*k],
                                   [-2*alpha*helicity*k,
                                    4*alpha+result.transverse_curvature*k*k]])
    for actual, expected in ((pulled_mass, expected_mass),
                             (pulled_potential, expected_potential)):
        for value in actual-expected:
            for order in range(3):
                assert s.simplify(s.expand(value).coeff(k, order)) == 0
    # Omitting either correction changes the actual optical gradient coefficient.
    assert s.simplify(s.expand(pulled_potential[1, 1]).coeff(k, 2)-c) != 0
    wrong_mass = s.Matrix([[1-mu*k*k/(2*rho), b*helicity*k/rho],
                           [0, 1-(mp-b*b/rho)*k*k/(2*j)]])
    assert s.simplify(s.expand((wrong_mass.T*mass*wrong_mass)[0, 1]).coeff(k, 1)) != 0
    assert transform.subs(k, 0) == s.eye(2)


def test_normal_form_exposes_structure_free_and_helicity_boundaries():
    with pytest.raises(ValueError, match="spin inertia"):
        micropolar_kinetic_normal_form(1, 0, 1, 1, 1, 1, 1, 1)
    with pytest.raises(ValueError, match="helicity"):
        micropolar_kinetic_normal_form(1, 1, 1, 1, 1, 1, 1, 1, 0)
