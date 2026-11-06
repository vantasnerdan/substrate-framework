"""Independent moment/quadrature proof and one actual-spectrum performance pilot."""
import argparse
import json
from pathlib import Path
from time import perf_counter
from tempfile import TemporaryDirectory
import numpy as np
import sympy as sp
from scipy.integrate import quad_vec
from scipy.linalg import expm
from scipy.special import eval_legendre
from density_moments import fourier_legendre_moments, SMALL_SERIES_TERMS
from density_integration import integrate_density_residual
from filtered_quadratic_bath import _fourier_legendre_moments, _chebyshev_to_legendre
from canonical_density_filter import load_canonical_density_filter
from consume_joint_influence import authenticate_handoff
from integrate_saved_spectral import load_saved_spectral
from substrate_framework.streamed_material_state import authenticate, file_binding
from reduced_projection_cache import ProjectionCache


def symbolic_identity():
    x = sp.Symbol('x'); z = sp.Symbol('z', real=True, nonzero=True)
    checks = []
    for order in range(4):
        polynomial = sp.Poly(sp.legendre(order, x), x)
        direct = sum(c*sp.diff(2*sp.sin(z)/z, z, power[0])/sp.I**power[0]
                     for power, c in polynomial.terms())
        reference = 2*sp.I**order*sp.expand_func(sp.jn(order, z))
        checks.append(sp.trigsimp(sp.expand(direct-reference)) == 0)
    k = SMALL_SERIES_TERMS+1
    tail = 1/sp.factorial(2*k+1)/(1-sp.Rational(1, (2*k+2)*(2*k+3)))
    assert all(checks) and tail < sp.Rational(1, 10**21)
    return dict(fourier_identities=checks, small_series_absolute_tail_bound=str(tail))


def moment_oracle():
    z = np.array([-60+.001j, -np.pi, -1.00000001, -.99999999, -1e-20+1e-20j,
                  0, 1e-20, .3+.4j, 1, 1.00000001+.1j, np.pi+.00001j,
                  7.2+.4j, 15.999, 16.0001, 35+.001j, -2j, 15j, 17j, -17j, 2.3])
    half = .4; gaps = z/half
    actual = fourier_legendre_moments(gaps, half, 0., 15)
    reference = _fourier_legendre_moments(gaps, half, 0., 15)
    assert np.allclose(actual, reference, rtol=2e-12, atol=2e-12)
    quadrature = np.array([quad_vec(lambda t: eval_legendre(n, t/half)*np.exp(1j*gaps*t),
                           -half, half, epsabs=2e-12, epsrel=2e-12)[0] for n in range(15)])
    assert np.allclose(actual, quadrature, rtol=2e-10, atol=2e-10)
    mutations = dict(real_projected_gaps=float(np.linalg.norm(
        fourier_legendre_moments(gaps.real, half, 0., 15)-reference)),
        wrong_phase_sign=float(np.linalg.norm(fourier_legendre_moments(-gaps, half, 0., 15)-reference)),
        wrong_jacobian=float(np.linalg.norm(fourier_legendre_moments(gaps, half/2, 0., 15)-reference)))
    assert min(mutations.values()) > 1e-6
    return dict(scipy_scaled_defect=float(np.max(np.abs(actual-reference)/(1+np.abs(reference)))),
                literal_quadrature_scaled_defect=float(np.max(np.abs(actual-quadrature)/(1+np.abs(quadrature)))),
                mutations=mutations, full_complex_arguments=True, degree14=True)


def convolution_oracle():
    rng = np.random.default_rng(4211)
    random = lambda shape: rng.normal(size=shape)+1j*rng.normal(size=shape)
    v = np.eye(3)+.12*random((3, 3)); vi = np.linalg.inv(v)
    energies = np.array([-18+.01j, .2-.04j, 21+.03j]); h = v@np.diag(energies)@vi
    saved = dict(right_eigenvectors=v, inverse_right_eigenvectors=vi,
                 energy_gaps=energies[:, None]-energies[None, :],
                 chebyshev_to_legendre=_chebyshev_to_legendre(8))
    column = random((3, 1)); column /= np.linalg.norm(column); p = column@column.conj().T
    vertices = random((3, 3, 3)); adjoints = random((3, 3, 3))
    residual = np.array([a-p@a@p for a in vertices])
    residual_adjoint = np.array([a-p@a@p for a in adjoints])
    bounds = np.array([[0., .31], [.31, .84]])
    coefficients = .07*random((2, 15, 3, 2))/(np.arange(15)[None, :, None, None]+1)**2
    coefficients[:, 14, 0, :] += .2
    provider = lambda mode: (vertices[mode], adjoints[mode])
    actual = integrate_density_residual(saved, provider, range(3), p, bounds, coefficients, row_block=2)
    direct = np.zeros((2, 2, 3, 3), complex)
    for segment, (left, right) in enumerate(bounds):
        def integrand(t):
            values = np.polynomial.chebyshev.chebval(2*(t-left)/(right-left)-1, coefficients[segment])
            force = np.einsum('mp,mij->pij', values, residual)
            adjoint = np.einsum('mp,mij->pij', values.conj(), residual_adjoint)
            return np.array([expm(1j*h*t)@force@expm(-1j*h*t),
                             expm(1j*h.conj().T*t)@adjoint@expm(-1j*h.conj().T*t)])
        direct += quad_vec(integrand, left, right, epsabs=2e-12, epsrel=2e-12)[0]
    errors = dict(forward=float(np.linalg.norm(actual['operators']-direct[0])),
                  independent_adjoint=float(np.linalg.norm(actual['adjoint_operators']-direct[1])))
    assert max(errors.values()) < 3e-10
    truncated = coefficients.copy(); truncated[:, 14] = 0
    low = integrate_density_residual(saved, provider, range(3), p, bounds, truncated)
    manufactured = integrate_density_residual(saved, lambda i:(vertices[i], vertices[i].conj().T),
                                               range(3), p, bounds, coefficients)
    mutations = dict(degree14_omission=float(np.linalg.norm(low['operators']-direct[0])),
                     manufactured_adjoint=float(np.linalg.norm(manufactured['adjoint_operators']-direct[1])))
    assert min(mutations.values()) > 1e-8
    return dict(errors=errors, mutations=mutations, probes=2, complex_nonhermitian_generator=True)


def actual_frequency_pilot(current_projection, original_snapshot, filter_receipt):
    producer, *_ = authenticate_handoff(current_projection, original_source_snapshots=[original_snapshot])
    saved = load_saved_spectral(producer['a37_filtered_operators'])
    packet = json.loads(Path(filter_receipt).read_bytes())
    filt = load_canonical_density_filter(packet)
    widths = np.diff(filt.density_segment_bounds, axis=1)[:, 0]
    segment = int(np.argmax(widths)); left, right = filt.density_segment_bounds[segment]
    start = perf_counter()
    old = _fourier_legendre_moments(saved['energy_gaps'], (right-left)/2, (right+left)/2, 15)
    old_seconds = perf_counter()-start
    start = perf_counter()
    new = fourier_legendre_moments(saved['energy_gaps'], (right-left)/2, (right+left)/2, 15)
    new_seconds = perf_counter()-start
    difference = float(np.max(np.abs(new-old)/(1+np.abs(old))))
    assert np.allclose(new, old, rtol=2e-12, atol=2e-12)
    return dict(original_spectral_source=producer['a37_filtered_operators'],
                original_filter=file_binding(Path(filter_receipt)), selected_widest_actual_segment=segment,
                segment_bounds=[float(left), float(right)], literal_gap_shape=list(saved['energy_gaps'].shape),
                all15_order_scaled_defect=difference, old_seconds=old_seconds, new_seconds=new_seconds,
                benchmark_scope='ONE actual widest accepted interval/all literal complex1038x1038 gaps; no H/eigensolve/state/old full-run replay')


def rectangular_recovery_oracle():
    rng = np.random.default_rng(4212)
    random = lambda shape: rng.normal(size=shape)+1j*rng.normal(size=shape)
    def coordinates(energies):
        v = np.eye(len(energies))+.1*random((len(energies), len(energies)))
        vi = np.linalg.inv(v)
        return dict(eigenvalues=energies, right_eigenvectors=v, inverse_right_eigenvectors=vi,
            energy_gaps=energies[:, None]-energies[None, :],
            chebyshev_to_legendre=_chebyshev_to_legendre(8)), v@np.diag(energies)@vi
    source, h = coordinates(np.array([-12+.02j, 2-.03j, 20+.04j]))
    target, ht = coordinates(np.array([-11-.02j, 22+.05j]))
    cs, ct = random((3, 1)), random((2, 1))
    ps = cs@cs.conj().T/np.linalg.norm(cs)**2
    pt = ct@ct.conj().T/np.linalg.norm(ct)**2
    vertices, adjoints = random((2, 3, 2)), random((2, 2, 3))
    residual = np.array([a-ps@a@pt for a in vertices])
    residual_adjoint = np.array([a-pt@a@ps for a in adjoints])
    coefficients = .02*random((1, 15, 2, 2)); coefficients[:, 14] += .2
    bounds = [[0., .6]]
    def forward(t):
        values = np.polynomial.chebyshev.chebval(2*t/.6-1, coefficients[0])
        return expm(1j*h*t)@np.einsum('mp,mij->pij', values, residual)@expm(-1j*ht*t)
    def backward(t):
        values = np.polynomial.chebyshev.chebval(2*t/.6-1, coefficients[0]).conj()
        return expm(1j*ht.conj().T*t)@np.einsum('mp,mij->pij', values, residual_adjoint)@expm(-1j*h.conj().T*t)
    exact_forward = quad_vec(forward, 0, .6, epsabs=2e-12, epsrel=2e-12)[0]
    exact_backward = quad_vec(backward, 0, .6, epsabs=2e-12, epsrel=2e-12)[0]
    provider = lambda mode:(vertices[mode], adjoints[mode])
    with TemporaryDirectory(prefix='p254-a42-row-recovery-oracle-') as directory:
        specification = dict(fixture='literal-complex-rectangular-two-probe-v1')
        cache = ProjectionCache(directory, specification)
        def interrupt(event):
            if event['stage'] == 'density_spectral_rows' and event['status'] == 'complete':
                raise InterruptedError('injected after first durable numeric block')
        try:
            integrate_density_residual(source, provider, range(2), ps, bounds, coefficients,
                target_saved=target, target_projector=pt, checkpoint=cache, row_block=2, progress=interrupt)
        except InterruptedError:
            pass
        assert len(cache.manifest['items']) == 1 and cache.manifest['status'] == 'partial'
        resumed = ProjectionCache(directory, specification)
        actual = integrate_density_residual(source, provider, range(2), ps, bounds, coefficients,
            target_saved=target, target_projector=pt, checkpoint=resumed, row_block=2)
        errors = dict(forward=float(np.linalg.norm(actual['operators']-exact_forward)),
                      independent_adjoint=float(np.linalg.norm(actual['adjoint_operators']-exact_backward)))
        assert max(errors.values()) < 3e-10
        assert actual['diagnostics']['new_spectral_row_blocks'] == 1
        assert actual['diagnostics']['reused_spectral_row_blocks'] == 1
        def forbidden_provider(mode):
            raise AssertionError('completed recovery must not replay physical vertices')
        complete = integrate_density_residual(source, forbidden_provider, range(2), ps, bounds, coefficients,
            target_saved=target, target_projector=pt,
            checkpoint=ProjectionCache(directory, specification), row_block=2)
        assert np.array_equal(actual['operators'], complete['operators'])
        assert np.array_equal(actual['adjoint_operators'], complete['adjoint_operators'])
        assert complete['diagnostics']['new_spectral_row_blocks'] == 0
        assert complete['diagnostics']['reused_spectral_row_blocks'] == 2
    return dict(errors=errors, partial_blocks_reused=1, missing_blocks_computed=1,
                completed_blocks_reused=2, completed_provider_replays=0, source_target_sizes=[3, 2])


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--current-projection', type=Path, required=True)
    parser.add_argument('--original-snapshot', type=Path, required=True)
    parser.add_argument('--filter-receipt', type=Path, required=True)
    args = parser.parse_args(); start = perf_counter()
    result = dict(schema='P254-A42-new-batched-complex-moment-oracle-v1',
                  symbolic=symbolic_identity(), moments=moment_oracle(), convolution=convolution_oracle(),
                  rectangular_recovery=rectangular_recovery_oracle())
    result['actual_frequency_pilot'] = actual_frequency_pilot(args.current_projection, args.original_snapshot,
                                                              args.filter_receipt)
    result.update(passed=True, seconds=perf_counter()-start,
                  scope='Encoded exact same finite-density moment/convolution algorithm; not new physics or all-source/wholeE accuracy')
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
