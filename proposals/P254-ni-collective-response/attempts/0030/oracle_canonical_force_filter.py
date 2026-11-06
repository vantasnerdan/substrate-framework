"""Owning exact-block and independent forced-ODE oracle for the causal filter.

Importing executes nothing. Run this file explicitly; JSON is emitted only
if every exact identity, forward convolution and sensitivity check passes.
Small dense/global blocks here are oracle-only, never a material state bank.
The independent reference does NOT call the consumer tangent or transpose.
No passive two-time noise, loss or lifetime assertion follows from this test.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / '0029'))
from reduced_electronic_correlations import ActiveBandCorrelations
from oracle_reduced_electronic_correlations import fixture, GlobalLinearization
from canonical_force_filter import build_canonical_force_filter


ERROR_LIMIT = 2e-8
SENSITIVITY_FLOOR = 2e-5


def _canonical(matrix):
    """Expand each exact block BEFORE composing it with the next block."""
    return matrix.applyfunc(sp.expand)


def _numeric(matrix):
    return np.asarray(matrix.tolist(), dtype=complex)


def exact_fixture():
    """Eight-dimensional exact complex moving model, with nonunit masses.

    Every electronic matrix is complex Hermitian but the bilinear dual is
    not an adjoint. rho(t), h(t) both move; screening, ionic feedback and
    global field feedback are all nonzero. All entries are rational in t.
    """
    t, r, i = sp.symbols('t', real=True), sp.Rational, sp.I
    rho = sp.Matrix([[r(2, 5)+t/100, r(1, 11)+i*(r(1, 9)+t/70)],
                     [r(1, 11)-i*(r(1, 9)+t/70), r(3, 5)-t/100]])
    h = sp.Matrix([[r(2, 7)+t/5, r(1, 3)+i*(r(1, 8)+t/9)],
                   [r(1, 3)-i*(r(1, 8)+t/9), -r(1, 4)+t*t/13]])
    g = sp.Matrix([[r(1, 2), r(2, 5)+i/3], [r(2, 5)-i/3, -r(1, 6)]])
    field = sp.Matrix([[r(2, 9), -r(1, 4)+2*i/7],
                       [-r(1, 4)-2*i/7, r(1, 5)]])
    potential = sp.Matrix([[r(1, 7), r(1, 8)-i/6],
                           [r(1, 8)+i/6, r(2, 11)]])
    kernel = r(3, 10)
    j = sp.diag(sp.Matrix([[0, 1], [-1, 0]]), sp.Matrix([[0, 1], [-1, 0]]))
    mechanical = _canonical(j*sp.diag(r(7, 5), r(1, 5), r(11, 6), r(1, 3)))
    rho, h, g, field, potential = map(_canonical, (rho, h, g, field, potential))
    exact = sp.zeros(8)
    for column in range(8):
        basis = sp.eye(8)[:, column]
        x, sigma = basis[:4, :], sp.Matrix(2, 2, list(basis[4:, :]))
        density = sp.expand(sp.trace(potential*sigma))
        dh = _canonical(g*x[0] + potential*(kernel*density) + field*x[2])
        gradient = sp.Matrix([sp.expand(sp.trace(g*sigma)), 0,
                              sp.expand(sp.trace(field*sigma)), 0])
        dx = _canonical(mechanical*x + j*gradient)
        ds = _canonical(-i*(_canonical(h*sigma)-_canonical(sigma*h)
                            + _canonical(dh*rho)-_canonical(rho*dh)))
        exact[:, column] = dx.col_join(sp.Matrix(4, 1, list(ds)))
    exact = _canonical(exact)
    # Canonical insertion is plus p-prime and plus Pi-prime, not a gradient
    # insertion via J, nor a mass-normalized coordinate acceleration.
    insertion = sp.zeros(8, 2)
    insertion[1, 0], insertion[3, 1] = 1, 1
    z = sp.Matrix(sp.symbols('z0:8'))
    state = sp.Matrix(sp.symbols('y0:8'))
    pairing = sp.expand((z.T*_canonical(exact*state))[0]
                       - (_canonical(exact.T*z).T*state)[0])
    if pairing != 0 or insertion.T*z != sp.Matrix([z[1], z[3]]):
        raise AssertionError('exact canonical bilinear/insertion identity failed')
    if (exact*insertion)[0, 0] != r(1, 5):
        raise AssertionError('physical mass belongs in L, not in the force insertion')
    rho_at = sp.lambdify(t, rho, 'numpy')
    h_at = sp.lambdify(t, h, 'numpy')
    matrix_at = sp.lambdify(t, exact, 'numpy')
    block = SimpleNamespace(
        ionic=_numeric(g)[None, None], ionic_reverse=_numeric(g.T)[None, None],
        potential=_numeric(potential)[None, :, :, None],
        density=_numeric(potential.T)[None, None], kernel=np.array([[float(kernel)]], complex),
    )
    c = SimpleNamespace(n=1, nq=1, ns=1, m=2, d=1, f=1,
                        target=np.array([[0]]), minus=np.array([0]),
                        mechanical=SimpleNamespace(generators=(_numeric(mechanical),),
                                                   j=(_numeric(j),)),
                        field_modes=_numeric(field)[None, None],
                        field_reverse_modes=_numeric(field.T)[None, None],
                        q=lambda q: block)

    def baseline(time):
        return np.asarray(rho_at(time), complex)[None], np.asarray(h_at(time), complex)[None]

    def reference_matrix(time):
        return np.asarray(matrix_at(time), complex)

    consumer = ActiveBandCorrelations(c, baseline, rtol=2e-11, atol=2e-13)
    # Compare the actual public implementation with exact rational blocks,
    # not merely a formal identity of two reference expressions.
    probe = sp.Matrix(8, 2, [r(k+1, 13)+i*r(17-k, 19) for k in range(16)])
    errors = {}
    for time in (r(0), r(2, 7), r(7, 5)):
        reference = _numeric(_canonical(exact.subs(t, time).T*probe))
        actual = consumer.transpose(0, float(time), _numeric(probe))
        errors[f'exact_transpose_at_{time}'] = _error(actual, reference)
        forward = consumer.tangent(0, float(time), _numeric(probe))
        errors[f'exact_tangent_at_{time}'] = _error(
            forward, _numeric(_canonical(exact.subs(t, time)*probe)))
    if max(errors.values()) > 3e-14:
        raise AssertionError(f'actual implementation differs from exact blocks: {errors}')
    return consumer, reference_matrix, errors


def _error(actual, expected):
    return float(np.linalg.norm(np.asarray(actual)-np.asarray(expected)) /
                 max(1.0, np.linalg.norm(expected)))


def _solve(rhs, interval, initial):
    if interval[0] == interval[1]:
        return np.asarray(initial).copy()
    result = solve_ivp(rhs, interval, np.asarray(initial).ravel(), method='DOP853',
                       rtol=8e-13, atol=8e-15, t_eval=[interval[1]],
                       first_step=min(1e-4, interval[1]-interval[0]))
    if not result.success:
        raise AssertionError(f'independent forward reference failed: {result.message}')
    return result.y[:, -1].reshape(np.shape(initial))


def _insertion(consumer, q):
    d, f = consumer.c.d, consumer.c.f if q == 0 else 0
    matrix = np.zeros((consumer.state_size(q), d+f), complex)
    matrix[d:2*d, :d] = np.eye(d)
    matrix[2*d+f:2*d+2*f, d:] = np.eye(f)
    return matrix


def _force(time, d, f):
    ionic = np.array([(1.0+.17j*(k+1))*np.cos((.7+.19*k)*time)
                      + (.2-.13j)*np.sin((1.3+.11*k)*time) for k in range(d)])
    field = np.array([(.31+.43j*(k+1))*np.exp((.11+.47j)*time)
                      + .23*np.sin(1.7*time) for k in range(f)])
    return np.concatenate((ionic, field))


def _convolution(filter_, force, *, transform=None):
    """Independent Gauss-Legendre integration of ONLY the public filter."""
    nodes, weights = np.polynomial.legendre.leggauss(20)
    answer = np.zeros(filter_.initial_dual.shape[1], complex)
    for left, right in filter_.segment_bounds:
        half, middle = (right-left)/2.0, (right+left)/2.0
        for node, weight in zip(nodes, weights):
            time = middle+half*node
            ionic, field = filter_.evaluate(time)
            kernels = np.concatenate((ionic, field), axis=0)
            if transform is not None:
                kernels = transform(time, kernels)
            answer += half*weight*(kernels.T@force(time))
    return answer


def _forward_forced(reference_matrix, insertion, final_time, initial, force, coordinates):
    terminal = _solve(lambda time, y: reference_matrix(time)@y + insertion@force(time),
                      (0.0, final_time), initial)
    return terminal[np.asarray(coordinates)]


def _forward_impulse(reference_matrix, insertion, time, final_time, coordinates):
    shape = insertion.shape
    terminal = _solve(lambda t, y: (reference_matrix(t)@y.reshape(shape)).ravel(),
                      (float(time), final_time), insertion)
    return terminal[np.asarray(coordinates)].T


def _forward_impulse_history(reference_matrix, insertion, final_time, coordinates):
    """One independent tiny forward fundamental solution for all probe times.

    Phi(T,s) B = Phi(T,0) solve(Phi(s,0), B). This avoids reintegrating the
    expensive independent global linearization at every boundary/interior
    witness; every original witness and the direct forced-ODE checks remain.
    The full transition bank here is tiny oracle-only data, never material data.
    """
    size = len(insertion)
    result = solve_ivp(
        lambda time, value: (reference_matrix(time) @ value.reshape(size, size)).ravel(),
        (0.0, final_time), np.eye(size, dtype=complex).ravel(), method='DOP853',
        rtol=8e-13, atol=8e-15, first_step=min(1e-4, final_time), dense_output=True)
    if not result.success:
        raise AssertionError(f'independent forward fundamental solve failed: {result.message}')
    terminal = result.y[:, -1].reshape(size, size)
    selected = np.asarray(coordinates)

    def impulse(time):
        initial = result.sol(float(time)).reshape(size, size)
        return (terminal @ np.linalg.solve(initial, insertion))[selected].T

    return impulse


def check_case(consumer, reference_matrix, q, final_time, coordinates, builder):
    """Endpoint/interior transfer, full preparation and forced convolution."""
    filter_ = builder(consumer, q, final_time, coordinates)
    d, f = consumer.c.d, consumer.c.f if q == 0 else 0
    nprobes, size = len(coordinates), consumer.state_size(q)
    insertion = _insertion(consumer, q)
    expected_shapes = {
        'initial_dual': (size, nprobes),
        'segment_bounds': (len(filter_.segment_bounds), 2),
        'ionic_chebyshev': (len(filter_.segment_bounds), 8, d, nprobes),
        'field_chebyshev': (len(filter_.segment_bounds), 8, f, nprobes),
        'boundary_times': (len(filter_.segment_bounds)+1,),
        'ionic_boundary': (len(filter_.segment_bounds)+1, d, nprobes),
        'field_boundary': (len(filter_.segment_bounds)+1, f, nprobes),
    }
    for key, shape in expected_shapes.items():
        if getattr(filter_, key).shape != shape:
            raise AssertionError(f'{key} has wrong compact array shape')
    if not (filter_.boundary_times[0] == 0.0 and
            filter_.boundary_times[-1] == final_time and
            np.all(np.diff(filter_.boundary_times) > 0.0)):
        raise AssertionError('accepted boundaries omit zero/terminal or reverse time')
    if not (np.array_equal(filter_.segment_bounds[:, 0], filter_.boundary_times[:-1]) and
            np.array_equal(filter_.segment_bounds[:, 1], filter_.boundary_times[1:])):
        raise AssertionError('segment bounds and accepted endpoints disagree')
    impulse_reference = _forward_impulse_history(
        reference_matrix, insertion, final_time, coordinates)
    times = list(filter_.boundary_times)
    for left, right in filter_.segment_bounds:
        times += [left+.23*(right-left), left+.71*(right-left)]
    # Both sides of an exact internal knot must select the physical interval,
    # rather than an interval in the solver's descending-time order.
    for time in filter_.boundary_times[1:-1]:
        times += [np.nextafter(time, 0.0), np.nextafter(time, final_time)]
    errors, references = {}, []
    for time in times:
        ionic, field = filter_.evaluate(time)
        reference = impulse_reference(time)
        actual = np.concatenate((ionic, field), axis=0)
        references.append((time, reference))
        errors[f'impulse_{len(references)-1}'] = _error(actual, reference)
    rng = np.random.default_rng(300011+q)
    initial = .3*(rng.normal(size=size)+1j*rng.normal(size=size))
    electronic_initial = initial.copy()
    electronic_initial[:len(consumer.c.mechanical.generators[q])] = 0.0
    zero_force = lambda time: np.zeros(d+f, complex)
    for name, state in (('full_preparation', initial), ('electronic_preparation', electronic_initial)):
        reference = _forward_forced(reference_matrix, insertion, final_time, state,
                                    zero_force, coordinates)
        errors[name] = _error(filter_.initial_dual.T@state, reference)
    for name, active in (('ionic', np.r_[np.ones(d), np.zeros(f)]),
                         ('field', np.r_[np.zeros(d), np.ones(f)]),
                         ('combined', np.ones(d+f))):
        if name == 'field' and not f:
            continue
        force = lambda time, active=active: active*_force(time, d, f)
        convolution = _convolution(filter_, force)
        reference = _forward_forced(reference_matrix, insertion, final_time,
                                    np.zeros(size, complex), force, coordinates)
        errors[f'{name}_forced_convolution'] = _error(convolution, reference)
    force = lambda time: _force(time, d, f)
    forced_reference = _forward_forced(reference_matrix, insertion, final_time,
                                       initial, force, coordinates)
    errors['preparation_plus_forcing'] = _error(
        filter_.initial_dual.T@initial+_convolution(filter_, force), forced_reference)
    errors['zero_time_force_vs_full_initial_dual'] = _error(
        np.concatenate(filter_.evaluate(0.0), axis=0), insertion.T@filter_.initial_dual)
    errors['terminal_displacement_force_is_zero'] = _error(
        np.concatenate(filter_.evaluate(final_time), axis=0), np.zeros((d+f, nprobes)))
    for invalid in (-1e-12, final_time+1e-12, np.nan, np.inf):
        try:
            filter_.evaluate(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError('filter extrapolates or accepts nonfinite insertion time')
    if max(errors.values()) > ERROR_LIMIT:
        raise AssertionError(f'independent forward oracle failed: {errors}')
    return filter_, errors, references, initial


class _WrongTranspose(ActiveBandCorrelations):
    """Deliberate numerical mutation, not an alternate implementation."""

    def __init__(self, original, reference, mutation):
        super().__init__(original.c, original.baseline, rtol=original.rtol, atol=original.atol)
        self.reference, self.mutation = reference, mutation

    def _transpose(self, q, time, value, block):
        matrix = self.reference(time)
        if self.mutation == 'forward_instead_of_transpose':
            return matrix@value
        if self.mutation == 'adjoint_instead_of_transpose':
            return matrix.conj().T@value
        if self.mutation == 'wrong_backward_sign':
            return -matrix.T@value
        if self.mutation == 'minus_q_swap':
            other = int(self.c.minus[q])
            return super()._transpose(other, time, value, self.c.q(other))
        raise AssertionError('unknown deliberate mutation')


def sensitivity(consumer, reference_matrix, q, final_time, coordinates,
                filter_, references, initial, builder):
    """Require actual end-to-end counterexamples, not source-text checks."""
    d, f = consumer.c.d, consumer.c.f if q == 0 else 0
    force = lambda time: _force(time, d, f)
    correct = filter_.initial_dual.T@initial+_convolution(filter_, force)
    gaps = {}
    mutations = ['forward_instead_of_transpose', 'adjoint_instead_of_transpose', 'wrong_backward_sign']
    if q != int(consumer.c.minus[q]):
        mutations.append('minus_q_swap')
    for mutation in mutations:
        wrong = builder(_WrongTranspose(consumer, reference_matrix, mutation), q, final_time, coordinates)
        gaps[mutation] = _error(wrong.initial_dual.T@initial+_convolution(wrong, force), correct)
    correct_convolution = _convolution(filter_, force)
    gaps['minus_physical_force'] = _error(-correct_convolution, correct_convolution)
    # A mass factor here duplicates the mass already inside L. Use distinct
    # nonunit masses in the counterexample, never the fixture's random weights.
    masses = np.arange(2, 2+d, dtype=float)
    gaps['extra_inverse_mass'] = _error(
        _convolution(filter_, force, transform=lambda t, k: k/np.r_[masses, np.ones(f)][:, None]),
        correct_convolution)
    gaps['conjugated_force_kernel'] = _error(
        _convolution(filter_, force, transform=lambda t, k: k.conj()), correct_convolution)
    insertion = _insertion(consumer, q)
    wrong = np.zeros_like(insertion)
    wrong[:d, :d] = np.eye(d)
    if f:
        wrong[2*d:2*d+f, d:] = np.eye(f)
    direct_wrong_selector = _forward_forced(reference_matrix, wrong, final_time,
                                            np.zeros(len(initial), complex), force, coordinates)
    gaps['coordinate_instead_of_momentum_selector'] = _error(direct_wrong_selector, correct_convolution)
    if f:
        only_field = lambda time: np.r_[np.zeros(d), _force(time, d, f)[d:]]
        field_correct = _convolution(filter_, only_field)
        field_wrong = _forward_forced(reference_matrix, wrong, final_time,
                                     np.zeros(len(initial), complex), only_field, coordinates)
        gaps['A_instead_of_Pi_selector'] = _error(field_wrong, field_correct)
    reversed_errors, reflected_errors, missing_boundary_errors = [], [], []
    for time, reference in references:
        reversed_value = np.concatenate(filter_.evaluate(final_time-time), axis=0)
        reversed_errors.append(_error(reversed_value, reference))
        if time not in filter_.boundary_times:
            segment = np.searchsorted(filter_.boundary_times, time)-1
            left, right = filter_.segment_bounds[segment]
            reflected = np.concatenate(filter_.evaluate(left+right-time), axis=0)
            reflected_errors.append(_error(reflected, reference))
        elif time != final_time:
            missing_boundary_errors.append(_error(np.zeros_like(reference), reference))
    gaps['descending_time_used_as_ascending'] = max(reversed_errors)
    gaps['reversed_local_polynomial_coordinate'] = max(reflected_errors)
    gaps['zero_or_missing_boundary_values'] = max(missing_boundary_errors)
    # Gamma physical mechanical transfer can be real despite complex electronic
    # dynamics. Adjoint sensitivity there is supplied by the preserved full
    # preparation dual; the non-Gamma case must expose conjugated force kernels.
    if q == 0:
        gaps.pop('conjugated_force_kernel')
    if min(gaps.values()) <= SENSITIVITY_FLOOR:
        raise AssertionError(f'oracle fixture does not expose a required mutation: {gaps}')
    return gaps


def run_oracle(builder=build_canonical_force_filter):
    exact_consumer, exact_matrix, symbolic_errors = exact_fixture()
    gamma, gamma_errors, gamma_references, gamma_initial = check_case(
        exact_consumer, exact_matrix, 0, 1.8, (0,), builder)
    gamma_gaps = sensitivity(exact_consumer, exact_matrix, 0, 1.8, (0,),
                             gamma, gamma_references, gamma_initial, builder)
    c, baseline, _ = fixture(n=3)
    consumer = ActiveBandCorrelations(c, baseline, rtol=2e-11, atol=2e-13)
    global_model = GlobalLinearization(c, baseline, 1)
    directed, directed_errors, directed_references, directed_initial = check_case(
        consumer, global_model.matrix, 1, 1.3, (2, 0), builder)
    directed_gaps = sensitivity(consumer, global_model.matrix, 1, 1.3, (2, 0),
                                directed, directed_references, directed_initial, builder)
    zero = builder(exact_consumer, 0, 0.0, (0,))
    if not (zero.segment_bounds.shape == (0, 2) and zero.boundary_times.tolist() == [0.0] and
            np.array_equal(zero.initial_dual, exact_consumer.ionic_probes(0, (0,))) and
            all(np.array_equal(v, np.zeros((1, 1))) for v in zero.evaluate(0.0))):
        raise AssertionError('zero-terminal-time filter mishandles the degenerate interval')
    # A strict finite-domain short baseline catches an unbounded initial-step
    # probe. It is a domain behavior check in addition to the physical oracle.
    short_time = 1e-8
    original_baseline = exact_consumer.baseline

    def bounded(time):
        if not 0.0 <= time <= short_time:
            raise AssertionError('solver queried outside the actual baseline interval')
        return original_baseline(time)

    short_consumer = ActiveBandCorrelations(exact_consumer.c, bounded,
                                            rtol=exact_consumer.rtol, atol=exact_consumer.atol)
    short = builder(short_consumer, 0, short_time, (0,))
    short_reference = _forward_impulse(exact_matrix, _insertion(exact_consumer, 0),
                                       0.0, short_time, (0,))
    short_error = _error(np.concatenate(short.evaluate(0.0), axis=0), short_reference)
    if short_error > 1e-14:
        raise AssertionError('short finite-domain impulse oracle failed')
    return {
        'schema': 'P254-canonical-force-filter-oracle-v1',
        'passed': True,
        'exact_blocks': 'expanded exact rational complex moving eight-state canonical model',
        'exact_block_errors': symbolic_errors,
        'gamma_forward_errors': gamma_errors,
        'directed_complex_forward_errors': directed_errors,
        'gamma_mutation_gaps': gamma_gaps,
        'directed_complex_mutation_gaps': directed_gaps,
        'short_baseline_domain_error': short_error,
        'error_limit': ERROR_LIMIT,
        'sensitivity_floor': SENSITIVITY_FLOOR,
        'fixtures': {'gamma': {'state_size': 8, 'moving_rho_and_h': True, 'ionic_mass': 5,
                               'field_kinetic_inverse_mass': '1/3', 'time': 1.8},
                     'directed': {'mesh': [3, 1, 1], 'q': 1, 'minus_q': 2,
                                  'spins': 2, 'bands': 2, 'coordinates': [2, 0], 'time': 1.3}},
        'scope': 'conditional causal physical-force filter only; not passive noise/loss or lifetime',
    }


if __name__ == '__main__':
    print(json.dumps(run_oracle(), indent=2, allow_nan=False))
