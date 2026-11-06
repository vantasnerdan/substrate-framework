"""Exact polynomial and literal Fock/CAR oracle for full force channels.

Importing runs nothing. Execute explicitly with PYTHONPATH=src. This fixture
uses the real CanonicalForceFilter arrays/evaluator, not an ODE replay or a
material packet. References integrate exact SymPy polynomials on nonuniform
intervals and lift actual dense operators to occupation-bit Fock space.
No reference calls the channel implementation's quadrature or moment helpers.
This is not a physical material-action, passive-kernel or loss admission.
"""
from __future__ import annotations

from dataclasses import replace
import json

import numpy as np
import sympy as sp

from substrate_framework.physical_factor_preparation import (
    ActionFactors, FactorFunctional, FactorPreparation,
)
from canonical_force_filter import CanonicalForceFilter
from force_channel_representation import decompose_force_channels


ATOL = 3e-11
RTOL = 3e-12
SENSITIVITY_FLOOR = 1e-6


def _canonical(matrix):
    """Canonicalize every exact block BEFORE the next matrix composition."""
    return matrix.applyfunc(lambda value: sp.expand_complex(sp.expand(value)))


def _numeric(matrix):
    return np.asarray(matrix.tolist(), dtype=complex)


def _maximum(value):
    return float(np.max(np.abs(value), initial=0.0))


def _check(actual, expected, label):
    actual, expected = np.asarray(actual), np.asarray(expected)
    if actual.shape != expected.shape or not np.allclose(actual, expected, atol=ATOL, rtol=RTOL):
        raise AssertionError(label+': shape or numerical invariant failed')
    return _maximum(actual-expected)


def _sensitive(correct, wrong, label):
    difference = _maximum(np.asarray(correct)-np.asarray(wrong))
    if difference <= SENSITIVITY_FLOOR:
        raise AssertionError(label+': fixture does not expose the mutation')
    return difference


def _integrate_local(matrix, x):
    """Exact integral on [-1,1] by monomial antiderivatives, not quadrature."""
    matrix = _canonical(matrix)
    return matrix.applyfunc(lambda value: sp.Add(*(
        coefficient*sp.Rational(2, power[0]+1)
        for power, coefficient in sp.Poly(value, x).terms() if power[0] % 2 == 0)))


def _chebyshev_blocks(polynomial, x):
    """Exact triangular polynomial conversion, never a numerical fit."""
    remainder = _canonical(polynomial)
    blocks = [sp.zeros(*polynomial.shape) for _ in range(8)]
    for degree in range(7, -1, -1):
        tk = sp.Poly(sp.chebyshevt(degree, x), x)
        block = remainder.applyfunc(
            lambda value: sp.Poly(value, x).nth(degree)/tk.LC())
        blocks[degree] = _canonical(block)
        remainder = _canonical(remainder-blocks[degree]*tk.as_expr())
    if remainder != sp.zeros(*polynomial.shape):
        raise AssertionError('exact polynomial conversion lost a coefficient')
    return blocks


def polynomial_fixture():
    """Degree-seven, continuous, complex two-probe kernels; one null ionic mode."""
    x = sp.symbols('x', real=True)
    rat, imaginary = sp.Rational, sp.I
    times = (rat(0), rat(1, 7), rat(11, 10), rat(9, 5))
    exact, coefficient_arrays, boundary_arrays = {}, {}, {}
    for kind, channels, offset in (('ionic', 4, 1), ('field', 2, 8)):
        def endpoint(knot):
            return _canonical(sp.Matrix(channels, 2, lambda channel, probe:
                0 if kind == 'ionic' and channel == 3 else
                rat((channel+1)*(knot+2)+offset+probe, 11)
                + imaginary*rat((channel+2)*(probe+1)-knot+offset, 13)))

        endpoints = [endpoint(knot) for knot in range(len(times))]
        segments, coefficients = [], []
        for segment in range(len(times)-1):
            def amplitude(shift):
                return _canonical(sp.Matrix(channels, 2, lambda channel, probe:
                    0 if kind == 'ionic' and channel == 3 else
                    rat((channel+1)*(probe+2)+shift+offset, 17+segment)
                    + imaginary*rat((channel+2)*(shift+1)-probe+offset, 19+2*segment)))
            # The T5 bubble produces non-negligible degree14 products. Other
            # bubbles and unequal endpoints break parity and channel degeneracy.
            linear = _canonical(((1-x)*endpoints[segment]+(1+x)*endpoints[segment+1])/2)
            bubble = _canonical((1-x*x)*(amplitude(1)*sp.chebyshevt(5, x)
                                         + amplitude(3)*sp.chebyshevt(2, x)
                                         + amplitude(6)*(x+rat(segment+1, 5))))
            polynomial = _canonical(linear+bubble)
            segments.append(polynomial)
            coefficients.append(np.stack([_numeric(block) for block in
                                           _chebyshev_blocks(polynomial, x)]))
        exact[kind] = segments
        coefficient_arrays[kind] = np.stack(coefficients)
        boundary_arrays[kind] = np.stack([_numeric(value) for value in endpoints])
    numeric_times = np.array([float(value) for value in times])
    canonical = CanonicalForceFilter(
        final_time=float(times[-1]), initial_dual=np.zeros((12, 2), complex),
        segment_bounds=np.column_stack((numeric_times[:-1], numeric_times[1:])),
        ionic_chebyshev=coefficient_arrays['ionic'],
        field_chebyshev=coefficient_arrays['field'], boundary_times=numeric_times,
        ionic_boundary=boundary_arrays['ionic'], field_boundary=boundary_arrays['field'],
        metadata={'q_index': 0, 'coordinates': [0, 1], 'polynomial_degree': 7,
                  'scope': 'exact finite oracle, not material filter authority'},
    )
    return canonical, x, times, exact


def exact_gram(left_segments, right_segments, x, times):
    gram = sp.zeros(left_segments[0].rows, right_segments[0].rows)
    for segment, (left, right) in enumerate(zip(left_segments, right_segments)):
        product = _canonical(left*_canonical(right.H))
        contribution = _canonical(_integrate_local(product, x)
                                  *((times[segment+1]-times[segment])/2))
        gram = _canonical(gram+contribution)
    return gram


def _modal_at(kind, canonical, time):
    knot = int(np.searchsorted(canonical.boundary_times, time, side='left'))
    if knot < len(canonical.boundary_times) and canonical.boundary_times[knot] == time:
        return kind['modal_boundary'][knot]
    segment = knot-1
    left, right = canonical.segment_bounds[segment]
    local = 2*(time-left)/(right-left)-1
    return np.polynomial.chebyshev.chebval(local, kind['modal_chebyshev'][segment])


def polynomial_checks(decompose):
    canonical, x, times, exact = polynomial_fixture()
    result = decompose(canonical)
    errors, sensitivities, references = {}, {}, {}
    _check(result['segment_bounds'], canonical.segment_bounds, 'accepted physical intervals')
    _check(result['boundary_times'], canonical.boundary_times, 'accepted knot times')
    for name in ('ionic', 'field'):
        kind, channels = result[name], exact[name][0].rows
        expected = _numeric(exact_gram(exact[name], exact[name], x, times))
        references[name] = expected
        errors[name+'_exact_integrated_gram'] = _check(kind['gram'], expected, name+' exact Gram')
        basis, eigenvalues = kind['basis'], kind['eigenvalues']
        if basis.shape != (channels, channels) or eigenvalues.shape != (channels,):
            raise AssertionError(name+': full channel count must survive including null modes')
        if np.any(np.diff(eigenvalues) > 0):
            raise AssertionError(name+': eigenvalues are not descending')
        errors[name+'_unitarity'] = _check(basis.conj().T@basis, np.eye(channels), name+' unitary basis')
        errors[name+'_spectral_reconstruction'] = _check(
            (basis*eigenvalues)@basis.conj().T, expected, name+' complete spectrum')
        errors[name+'_all_coefficients'] = _check(
            np.einsum('ia,skap->skip', basis, kind['modal_chebyshev']),
            canonical.arrays[name+'_chebyshev'], name+' coefficient reconstruction')
        errors[name+'_all_boundaries'] = _check(
            np.einsum('ia,bap->bip', basis, kind['modal_boundary']),
            canonical.arrays[name+'_boundary'], name+' knot reconstruction')
        # Interior asymmetry exposes x -> -x, which a symmetric integrated
        # Gram alone CANNOT expose. References evaluate exact local blocks.
        for segment, polynomial in enumerate(exact[name]):
            for local in (sp.Rational(-3, 5), sp.Rational(2, 7)):
                physical_time = times[segment]+(times[segment+1]-times[segment])*(local+1)/2
                expected_value = _numeric(_canonical(polynomial.subs(x, local)))
                actual_modal = _modal_at(kind, canonical, float(physical_time))
                label = name+f'_interior_{segment}_{local}'
                errors[label] = _check(basis@actual_modal, expected_value, label)
        for knot, time in enumerate(canonical.boundary_times):
            errors[name+f'_boundary_evaluation_{knot}'] = _check(
                basis@_modal_at(kind, canonical, time),
                canonical.arrays[name+'_boundary'][knot], name+' exact knot evaluation')
        # Explicit independently computed bad quadratures/weights demonstrate
        # that degree7 and physical dt, rather than sampled channel norms, matter.
        nodes7, weights7 = np.polynomial.legendre.leggauss(7)
        wrong7 = np.zeros_like(expected)
        no_jacobian = np.zeros_like(expected)
        for segment, polynomial in enumerate(exact[name]):
            evaluator = sp.lambdify(x, polynomial, 'numpy')
            for node, weight in zip(nodes7, weights7):
                value = np.asarray(evaluator(node), complex)
                wrong7 += float((times[segment+1]-times[segment])/2)*weight*(value@value.conj().T)
            product = _canonical(polynomial*_canonical(polynomial.H))
            no_jacobian += _numeric(_integrate_local(product, x))
        sensitivities[name+'_seven_point_rule'] = _sensitive(expected, wrong7, name+' wrong quadrature')
        sensitivities[name+'_missing_dt_jacobian'] = _sensitive(expected, no_jacobian, name+' wrong time weights')
        sensitivities[name+'_reversed_signed_intervals'] = _sensitive(expected, -expected, name+' signed interval orientation')
        original = _numeric(_canonical(exact[name][1].subs(x, sp.Rational(2, 7))))
        reversed_local = _numeric(_canonical(exact[name][1].subs(x, -sp.Rational(2, 7))))
        sensitivities[name+'_reversed_local_coordinate'] = _sensitive(original, reversed_local, name+' reversed interval polynomial')
    # Changing the field unit scale must not rotate or reweight ionic modes;
    # the reciprocal test also checks that there is no hidden pooled metric.
    for changed, untouched in (('field', 'ionic'), ('ionic', 'field')):
        scaled = replace(canonical, **{
            changed+'_chebyshev': canonical.arrays[changed+'_chebyshev']*17,
            changed+'_boundary': canonical.arrays[changed+'_boundary']*17,
        })
        scaled_result = decompose(scaled)
        errors[changed+'_unit_scale_leaves_'+untouched+'_gram'] = _check(
            scaled_result[untouched]['gram'], result[untouched]['gram'], 'unit-separated Gram')
        errors[changed+'_unit_scale_leaves_'+untouched+'_eigenvalues'] = _check(
            scaled_result[untouched]['eigenvalues'], result[untouched]['eigenvalues'], 'unit-separated spectrum')
        errors[changed+'_gram_scale'] = _check(scaled_result[changed]['gram'],
                                              result[changed]['gram']*17**2, 'physical scale squared')
    # This off-kind integral is deliberately ONLY a unit-invalid mutation
    # witness: it is not added to either physical Gram or output force.
    forbidden_cross = _numeric(exact_gram(exact['ionic'], exact['field'], x, times))
    sensitivities['unlike_units_have_nonzero_forbidden_cross_gram'] = _sensitive(
        forbidden_cross, np.zeros_like(forbidden_cross), 'ionic/field mixing')
    empty = replace(canonical, field_chebyshev=np.empty((3, 8, 0, 2), complex),
                    field_boundary=np.empty((4, 0, 2), complex),
                    metadata={**canonical.metadata, 'q_index': 1})
    empty_result = decompose(empty)
    for key, shape in (('basis', (0, 0)), ('gram', (0, 0)), ('eigenvalues', (0,)),
                       ('modal_chebyshev', (3, 8, 0, 2)), ('modal_boundary', (4, 0, 2))):
        if empty_result['field'][key].shape != shape:
            raise AssertionError('empty field lost its explicit '+key+' axes')
    errors['empty_field_leaves_ionic_gram'] = _check(
        empty_result['ionic']['gram'], references['ionic'], 'non-Gamma ionic representation')
    zero = replace(canonical, final_time=0.0, segment_bounds=np.empty((0, 2)),
                   boundary_times=np.array([0.0]),
                   ionic_chebyshev=np.empty((0, 8, 4, 2), complex),
                   field_chebyshev=np.empty((0, 8, 2, 2), complex),
                   ionic_boundary=canonical.ionic_boundary[:1],
                   field_boundary=canonical.field_boundary[:1])
    zero_result = decompose(zero)
    for name, channels in (('ionic', 4), ('field', 2)):
        errors[name+'_zero_time_gram'] = _check(zero_result[name]['gram'],
                                               np.zeros((channels, channels)), 'zero-time integral')
        errors[name+'_zero_time_full_basis'] = _check(
            zero_result[name]['basis'].conj().T@zero_result[name]['basis'],
            np.eye(channels), 'zero-time full basis')
        errors[name+'_zero_time_boundary'] = _check(
            np.einsum('ia,bap->bip', zero_result[name]['basis'], zero_result[name]['modal_boundary']),
            zero.arrays[name+'_boundary'], 'zero-time boundary retained')
    return canonical, result, errors, sensitivities


def _annihilators(modes):
    """Actual occupation-bit CAR matrices, with explicit fermionic signs."""
    answer = []
    for mode in range(modes):
        matrix = sp.zeros(2**modes)
        for occupation in range(2**modes):
            if occupation & (1 << mode):
                sign = (-1)**((occupation & ((1 << mode)-1)).bit_count())
                matrix[occupation ^ (1 << mode), occupation] = sign
        answer.append(matrix)
    return answer


def _lift(operator, car):
    result = sp.zeros(car[0].rows)
    for row in range(operator.rows):
        for column in range(operator.cols):
            pair = _canonical(car[row].H*car[column])
            result = _canonical(result+operator[row, column]*pair)
    return result


def physical_car_checks(canonical, representation):
    """Literal dense physical forces -> independent adjoints -> public CAR API.

The complex nonorthogonal raw frame and its density coefficients are left
unchanged. A two-orbital Gaussian mixture has a third empty PW orbital;
forces connect to that orbital, so saved-band-only covariance is incorrect.
Both directed orders and ionic/field cross covariances remain present.
"""
    rat, imaginary = sp.Rational, sp.I
    orbitals = sp.Matrix([[rat(3, 5), 4*imaginary/5, 0],
                          [4*imaginary/5, rat(3, 5), 0], [0, 0, 1]])
    raw_map = sp.Matrix([[1, imaginary/4], [rat(1, 7), 1]])
    inverse = _canonical(raw_map.inv())
    population = sp.diag(rat(2, 5), rat(3, 7))
    psi_exact = _canonical(orbitals[:, :2]*raw_map)
    left_density = _canonical(inverse*population)
    density_exact = _canonical(left_density*inverse.H)
    # Occupations are specified in the full orbitals frame. The empty third
    # mode is represented literally, not by a Wick/dense covariance formula.
    fock = sp.zeros(8)
    for bits in range(4):
        fock[bits, bits] = sp.prod(population[j, j] if bits & (1 << j)
                                  else 1-population[j, j] for j in range(2))
    car = _annihilators(3)
    for row in range(3):
        for column in range(3):
            anticommutator = _canonical(car[row]*car[column].H+car[column].H*car[row])
            if anticommutator != (sp.eye(8) if row == column else sp.zeros(8)):
                raise AssertionError('literal Fock fixture fails CAR')
    physical_exact = {}
    for kind, channels, offset in (('ionic', 4, 1), ('field', 2, 7)):
        for channel in range(channels):
            index = channel+offset
            operator = sp.Matrix([
                [rat(index, 5), rat(index+1, 7)+imaginary*rat(index+2, 9), rat(2-index, 8)+imaginary/3],
                [rat(index+1, 7)-imaginary*rat(index+2, 9), -rat(index+1, 6), rat(index, 4)-2*imaginary/5],
                [rat(2-index, 8)-imaginary/3, rat(index, 4)+2*imaginary/5, rat(index+3, 11)],
            ])
            physical_exact[kind, channel] = _canonical(operator)
    exact_lifts, exact_means = {}, {}
    for key, operator in physical_exact.items():
        left = _canonical(orbitals.H*operator)
        orbital_operator = _canonical(left*orbitals)
        exact_lifts[key] = _lift(orbital_operator, car)
        exact_means[key] = sp.expand(sp.trace(_canonical(fock*exact_lifts[key])))
    exact_covariances = {}
    for left in physical_exact:
        first = _canonical(fock*exact_lifts[left])
        for right in physical_exact:
            product = _canonical(first*exact_lifts[right])
            exact_covariances[left, right] = sp.expand(
                sp.trace(product)-exact_means[left]*exact_means[right])
    psi, density = _numeric(psi_exact), _numeric(density_exact)
    operators, kinds = {}, {}
    for kind in ('ionic', 'field'):
        originals = np.stack([_numeric(physical_exact[kind, channel])
                              for channel in range(representation[kind]['basis'].shape[0])])
        basis = representation[kind]['basis']
        modes = np.einsum('ia,ijk->ajk', basis, originals)
        kinds[kind] = originals, modes
        for label, matrices in (('channel', originals), ('mode', modes)):
            for index, operator in enumerate(matrices):
                operators[kind, label, index] = operator
    # Physical modal matrices are generally NOT Hermitian; supply the actual
    # dense Hilbert adjoint action, never reuse action for a complex direction.
    probes = {}
    for key, operator in operators.items():
        identity = ':'.join(map(str, key))
        probes[identity] = ActionFactors(operator@psi, operator.conj().T@psi,
            identity, identity+':actual-adjoint',
            'Ry/bohr' if key[0] == 'ionic' else 'Ry/field-coordinate',
            'three-orbital-finite-PW', 'literal dense force fixture')
    functional = FactorFunctional.from_actions(psi, probes,
        factor_identity='unchanged-nonorthogonal-frame', pw_basis_identity='three-orbital-finite-PW')
    preparation = FactorPreparation(functional, density, 'literal original Gaussian Fock ensemble')
    fock_numeric = _numeric(fock)
    numeric_car = [_numeric(value) for value in car]
    orbital_numeric = _numeric(orbitals)
    lifted, means = {}, {}
    for key, operator in operators.items():
        rotated = orbital_numeric.conj().T@operator@orbital_numeric
        lift = sum((rotated[row, column]*(numeric_car[row].conj().T@numeric_car[column])
                    for row in range(3) for column in range(3)), np.zeros((8, 8), complex))
        lifted[key], means[key] = lift, np.trace(fock_numeric@lift)
    errors, sensitivities = {}, {}
    covariance = {}
    for left in operators:
        label_left = ':'.join(map(str, left))
        errors['mean_'+label_left] = _check(preparation.mean(label_left), means[left], 'actual physical mean')
        for right in operators:
            label_right = ':'.join(map(str, right))
            independent = np.trace(fock_numeric@lifted[left]@lifted[right])-means[left]*means[right]
            actual = preparation.ordered(label_left, label_right)
            covariance[left, right] = independent
            errors['CAR_'+label_left+'__'+label_right] = _check(actual, independent, 'actual ordered CAR covariance')
            if left[1] == right[1] == 'channel':
                exact = complex(exact_covariances[(left[0], left[2]), (right[0], right[2])])
                errors['exact_'+label_left+'__'+label_right] = _check(actual, exact, 'exact SymPy Fock covariance')
    sample_times = (float(sp.Rational(3, 7)), float(sp.Rational(13, 10)))
    # Check bilinear operators, not just orthonormality/Gram energy. Different
    # times and probes, ordered pairs and cross-kind blocks all participate.
    physical_values = [canonical.evaluate(time) for time in sample_times]
    worst_wrong_operator, worst_wrong_covariance, worst_missing_cross = 0.0, 0.0, 0.0
    for kind_index, kind in enumerate(('ionic', 'field')):
        originals, modes = kinds[kind]
        basis = representation[kind]['basis']
        wrong_modes = np.einsum('ia,ijk->ajk', basis.conj(), originals)
        for time_index, time in enumerate(sample_times):
            kernel = physical_values[time_index][kind_index]
            modal = _modal_at(representation[kind], canonical, time)
            for probe in range(2):
                direct = np.einsum('i,ijk->jk', kernel[:, probe], originals)
                transformed = np.einsum('a,ajk->jk', modal[:, probe], modes)
                errors[kind+f'_physical_operator_{time_index}_{probe}'] = _check(
                    transformed, direct, 'nonconjugating physical force reconstruction')
                wrong = np.einsum('a,ajk->jk', modal[:, probe], wrong_modes)
                worst_wrong_operator = max(worst_wrong_operator, _maximum(wrong-direct))
    for left_index, left_kind in enumerate(('ionic', 'field')):
        a = _modal_at(representation[left_kind], canonical, sample_times[0])
        left_kernel = physical_values[0][left_index]
        left_originals, left_modes = kinds[left_kind]
        left_basis = representation[left_kind]['basis']
        for right_index, right_kind in enumerate(('ionic', 'field')):
            b = _modal_at(representation[right_kind], canonical, sample_times[1])
            right_kernel = physical_values[1][right_index]
            right_originals, right_modes = kinds[right_kind]
            right_basis = representation[right_kind]['basis']
            original_cov = np.array([[covariance[(left_kind, 'channel', i), (right_kind, 'channel', j)]
                                     for j in range(len(right_originals))] for i in range(len(left_originals))])
            modal_cov = np.array([[covariance[(left_kind, 'mode', i), (right_kind, 'mode', j)]
                                  for j in range(len(right_modes))] for i in range(len(left_modes))])
            direct = left_kernel.T@original_cov@right_kernel
            transformed = a.T@modal_cov@b
            label = left_kind+'__'+right_kind
            errors[label+'_two_time_two_probe_ordered_covariance'] = _check(
                transformed, direct, 'ordered output covariance under full complex basis')
            # A conjugated physical direction induces the wrong (conjugated)
            # channel mixing, even though the coefficient basis is unitary.
            wrong_cov = left_basis.conj().T@original_cov@right_basis.conj()
            worst_wrong_covariance = max(worst_wrong_covariance, _maximum(a.T@wrong_cov@b-direct))
            if left_kind == right_kind:
                worst_missing_cross = max(worst_missing_cross,
                    _maximum(a.T@np.diag(np.diag(modal_cov))@b-direct))
            else:
                sensitivities[label+'_cross_kind_covariance_is_retained'] = _sensitive(
                    direct, np.zeros_like(direct), 'dropping unlike-kind covariance (not unit mixing)')
    sensitivities['conjugated_physical_mode_operator'] = _sensitive(
        worst_wrong_operator, 0.0, 'conjugated physical mode operator')
    sensitivities['conjugated_physical_mode_CAR_covariance'] = _sensitive(
        worst_wrong_covariance, 0.0, 'conjugated physical mode covariance')
    sensitivities['dropping_modal_cross_covariances'] = _sensitive(
        worst_missing_cross, 0.0, 'missing modal covariance cross terms')
    return errors, sensitivities, len(operators), len(covariance)


def run_oracle(decompose=decompose_force_channels):
    canonical, result, polynomial_errors, polynomial_sensitivities = polynomial_checks(decompose)
    car_errors, car_sensitivities, operators, directed_pairs = physical_car_checks(canonical, result)
    return {
        'schema': 'P254-full-force-channel-representation-oracle-v1',
        'passed': True,
        'exact_polynomial_max_abs_errors': polynomial_errors,
        'independent_physical_CAR_max_abs_errors': car_errors,
        'mutation_sensitivities': {**polynomial_sensitivities, **car_sensitivities},
        'literal_physical_operators_including_complex_modes': operators,
        'independent_directed_Fock_covariances': directed_pairs,
        'ode_solves': 0,
        'scope': 'full-mode force representation and finite physical-operator/CAR invariance; no material action, full passive kernel or physical loss claim',
    }


if __name__ == '__main__':
    print(json.dumps(run_oracle(), indent=2, sort_keys=True, allow_nan=False))
