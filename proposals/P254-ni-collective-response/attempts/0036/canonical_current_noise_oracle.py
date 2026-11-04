"""Exact small-Fock oracle for the initial canonical uniform-field force boundary.

Execution is explicit: python canonical_current_noise_oracle.py --output NEWPATH.
The independent scalar/state/operator construction finishes before the adjacent
consumer is imported. No material run, old qualification proof, or initialization
is replayed. The Fock basis helpers are reused from the existing A30 oracle;
none of that oracle's reference calculations or production interfaces is called.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import importlib.util
import itertools
import json
import os
from pathlib import Path
import sys
import traceback

import numpy as np
import sympy as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CONTRACT = HERE / 'implementation-contract-001.json'
CONTRACT_SHA256 = 'f6f8a3e9b182db6fcf5a502fdf4d762f71d2c382efd98faf61a89c12136fafdd'
HELPER = HERE.parent / '0030' / 'oracle_material_noise_projection.py'
HELPER_SHA256 = '2893c6b35cffd2244564efffec2af9af003f59f545ef8f87852631de289dab19'
ATOL = 2e-12
RTOL = 2e-11
PARTS = ('R', 'E')


def require(predicate, message):
    if not bool(predicate):
        raise ValueError(message)


def binding(path, expected_sha256=None, expected_bytes=None):
    path = Path(path).resolve()
    data = path.read_bytes()
    result = dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
    if expected_sha256 is not None:
        require(result['sha256'] == expected_sha256, f'source digest mismatch: {path}')
    if expected_bytes is not None:
        require(result['bytes'] == expected_bytes, f'source length mismatch: {path}')
    return result


def authenticate_sources():
    plan_binding = binding(CONTRACT, CONTRACT_SHA256)
    plan = json.loads(CONTRACT.read_text())
    sources = [binding(row['path'], row['sha256'], row['bytes'])
               for row in plan['settled_sources']]
    helper_binding = binding(HELPER, HELPER_SHA256, 15773)
    return dict(preregistration=plan_binding, settled_sources=sources,
                fock_helpers=helper_binding, oracle=binding(__file__)), plan


def fock_helpers():
    """Reuse only occupation-bit annihilators and the existing symbolic lift."""
    src = str(ROOT / 'src')
    if src not in sys.path:
        sys.path.insert(0, src)
    spec = importlib.util.spec_from_file_location('_a36_existing_fock_helpers', HELPER)
    require(spec is not None and spec.loader is not None, 'Fock helper loader unavailable')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.annihilators, module.symbolic_lift


def clean(value):
    if isinstance(value, sp.MatrixBase):
        return value.applyfunc(sp.simplify)
    return sp.simplify(value)


def zero(value):
    if isinstance(value, sp.MatrixBase):
        return all(sp.simplify(entry) == 0 for entry in value)
    return sp.simplify(value) == 0


def exact(value):
    if isinstance(value, sp.MatrixBase):
        return [[sp.sstr(value[i, j]) for j in range(value.cols)]
                for i in range(value.rows)]
    return sp.sstr(value)


def numeric(value):
    return np.asarray(value.evalf(17).tolist(), dtype=np.complex128)


def complex_rows(value):
    array = np.asarray(value, dtype=np.complex128)
    return {'real': array.real.tolist(), 'imag': array.imag.tolist()}


def comparison(actual, expected):
    actual = np.asarray(actual, dtype=np.complex128)
    expected = np.asarray(expected, dtype=np.complex128)
    if actual.shape != expected.shape:
        return dict(passes=False, actual_shape=list(actual.shape),
                    expected_shape=list(expected.shape), shape_mismatch=True)
    difference = actual-expected
    finite = bool(np.all(np.isfinite(actual)))
    return dict(passes=finite and bool(np.allclose(actual, expected, atol=ATOL, rtol=RTOL)),
                finite=finite, max_abs_difference=float(np.max(np.abs(difference))),
                difference=complex_rows(difference))


def gaussian_fock(unitary, occupations, car):
    """Explicit independent-orbital Slater mixture, not a Wick covariance."""
    require(zero(unitary.H*unitary-sp.eye(2)), 'fixture orbitals must be unitary')
    require(all(0 < p < 1 for p in occupations), 'explicit mixed occupations required')
    dimension = car[0].rows
    creators = [sum((unitary[i, a]*car[i].H for i in range(2)), sp.zeros(dimension))
                for a in range(2)]
    density = sp.zeros(dimension)
    probabilities = []
    for bits in itertools.product((0, 1), repeat=2):
        probability = sp.prod(p if bit else 1-p for p, bit in zip(occupations, bits))
        vector = sp.zeros(dimension, 1)
        vector[0] = 1
        for a in reversed(range(2)):
            if bits[a]:
                vector = creators[a]*vector
        require(zero((vector.H*vector)[0]-1), 'Slater vector normalization failed')
        density += probability*vector*vector.H
        probabilities.append(probability)
    density = clean(density)
    require(zero(sp.trace(density)-1) and zero(density-density.H), 'invalid Gaussian state')
    require(sum(probabilities) == 1 and all(p > 0 for p in probabilities),
            'Gaussian mixture positivity failed')
    return density


def mean(state, operators):
    return clean(sp.Matrix([sp.trace(state*operator) for operator in operators]))


def covariance(state, left, right=None):
    """Direct ordered AB trace on an actual Fock density matrix."""
    right = left if right is None else right
    ml, mr = mean(state, left), mean(state, right)
    return clean(sp.Matrix(len(left), len(right),
                           lambda i, j: sp.trace(state*left[i]*right[j])-ml[i]*mr[j]))


def fixture():
    """Two explicit sectors at the actual N64 weight, an unrenormalized partial."""
    r, i = sp.Rational, sp.I
    ncells = 64
    weights = (r(1, ncells), r(1, ncells))
    polarizations = sp.Matrix([[r(3, 5), r(4, 5), 0],
                               [-r(12, 25), r(9, 25), r(4, 5)]])
    require(zero(polarizations*polarizations.T-sp.eye(2)), 'nonorthonormal fixture')
    require(not zero(polarizations-sp.eye(3)[:2, :]), 'trivial polarization fixture')
    frames = (sp.Matrix([[r(3, 5), 4*i/5], [4*i/5, r(3, 5)]]),
              sp.Matrix([[r(5, 13), r(12, 13)], [-12*i/13, 5*i/13]]))
    occupations = ((r(1, 4), r(2, 3)), (r(2, 5), r(3, 4)))
    h0 = (sp.Matrix([[1, (1+i)/3], [(1-i)/3, -r(2, 5)]]),
          sp.Matrix([[-r(2, 3), (2-i)/5], [(2+i)/5, r(4, 7)]]))
    vertices = ((sp.Matrix([[1, 1+i], [1-i, -2]]),
                 sp.Matrix([[r(2, 3), 2-i], [2+i, r(1, 5)]]),
                 sp.Matrix([[-1, 1+2*i], [1-2*i, 3]])),
                (sp.Matrix([[2, -1+2*i], [-1-2*i, -r(1, 3)]]),
                 sp.Matrix([[-r(3, 5), 3+i], [3-i, 1]]),
                 sp.Matrix([[r(4, 7), 2-3*i], [2+3*i, -2]])))
    return ncells, weights, polarizations, frames, occupations, h0, vertices


def independent_reference():
    """Finish the exact differentiated scalar and Fock answer before any map import."""
    annihilators, lift = fock_helpers()
    n, weights, e, frames, occupations, h0, vertices = fixture()
    a = sp.symbols('a_x a_y a_z', real=True)
    q = sp.symbols('q_0 q_1', real=True)
    shift = (sp.Matrix(1, len(q), q)*e).tolist()[0]
    substitution = dict(zip(a, shift))
    car = [sp.Matrix(item) for item in annihilators(2)]
    global_car = [sp.Matrix(item) for item in annihilators(4)]
    states = [gaussian_fock(u, d, car) for u, d in zip(frames, occupations)]
    # Occupation-bit mode zero is the low bit: cell zero is the right tensor factor.
    state = sp.kronecker_product(states[1], states[0])
    require(zero(sp.trace(state)-1), 'global Fock state is not normalized')
    projector = sp.diag(1, 0)
    global_h = sp.zeros(16)
    global_parts_h = {part: sp.zeros(16) for part in PARTS}
    current_totals = [sp.zeros(16) for _ in a]
    currents_by_sector = []
    local_covariances = []
    local_blocks = []
    sector_records = []
    derived_maps = []
    for s in range(2):
        require(all(zero(l-l.H) for l in vertices[s]), 'non-Hermitian physical vertex')
        require(any(not zero(x*y-y*x) for x, y in itertools.combinations(vertices[s], 2)),
                'fixture must have noncommuting vertices')
        h = h0[s]+sum((axis*l for axis, l in zip(a, vertices[s])), sp.zeros(2))
        require(zero(h-h.H), 'one-body scalar must be Hermitian for real field')
        local_h0 = lift(h0[s], car)
        require(not zero(states[s]*local_h0-local_h0*states[s]),
                'fixture must be coherent and nonstationary, not equilibrium')
        rho = clean(frames[s]*sp.diag(*occupations[s])*frames[s].H)
        require(rho[0, 1] != 0 and not zero(projector*rho-rho*projector),
                'retained projector must not commute with prepared coherent density')
        currents = [clean(-h.diff(axis)/2) for axis in a]
        # Source-vector coefficients are solved from differentiated operators,
        # not set equal to the proposed 2*E map.
        force_one_body = [clean(-h.subs(substitution).diff(coordinate)) for coordinate in q]
        basis = sp.Matrix.hstack(*(sp.Matrix(list(j)) for j in currents))
        require(basis.rank() == 3, 'current basis cannot expose all source-vector errors')
        rows = []
        for force in force_one_body:
            coefficients, free = basis.gauss_jordan_solve(sp.Matrix(list(force)))
            require(free.rows == 0, 'force source-vector coefficients are not unique')
            require(zero(sum((c*j for c, j in zip(coefficients, currents)), sp.zeros(2))-force),
                    'differentiated scalar source-vector reconstruction failed')
            rows.append(list(coefficients))
        derived_maps.append(sp.Matrix(rows))
        local_j = [lift(j, car) for j in currents]
        parts = {'R': [clean(projector*j*projector) for j in currents]}
        parts['E'] = [clean(j-r) for j, r in zip(currents, parts['R'])]
        lifted_parts = {part: [lift(j, car) for j in values] for part, values in parts.items()}
        local_c = covariance(states[s], local_j)
        blocks = {p+t: covariance(states[s], lifted_parts[p], lifted_parts[t])
                  for p in PARTS for t in PARTS}
        require(zero(sum(blocks.values(), sp.zeros(3))-local_c),
                'local ordered split covariance reconstruction failed')
        local_covariances.append(local_c)
        local_blocks.append(blocks)
        global_sector_car = global_car[2*s:2*s+2]
        hq = h.subs(substitution)
        global_h += n*weights[s]*lift(hq, global_sector_car)
        retained_h = projector*hq*projector
        for part, scalar_part in (('R', retained_h), ('E', hq-retained_h)):
            global_parts_h[part] += n*weights[s]*lift(scalar_part, global_sector_car)
        sector_j = [lift(j, global_sector_car) for j in currents]
        currents_by_sector.append(sector_j)
        for axis in range(3):
            current_totals[axis] += n*weights[s]*sector_j[axis]
        sector_records.append(dict(index=s, weight=exact(weights[s]), frame=exact(frames[s]),
                                   occupations=[exact(p) for p in occupations[s]],
                                   gaussian_fock_state=exact(states[s]),
                                   scalar_h_ry=exact(h), current_ha=tuple(exact(j) for j in currents),
                                   nonstationary_commutator=exact(clean(states[s]*local_h0-local_h0*states[s])),
                                   raw_current_covariance=exact(local_c),
                                   ordered_blocks={key: exact(value) for key, value in blocks.items()}))
    require(zero(derived_maps[0]-derived_maps[1]), 'sector-dependent source-vector convention')
    direct_force = [clean(-global_h.diff(coordinate)) for coordinate in q]
    direct_parts = {part: [clean(-h.diff(coordinate)) for coordinate in q]
                    for part, h in global_parts_h.items()}
    force_mean = mean(state, direct_force)
    differentiated_mean_scalar = clean(sp.Matrix([-sp.diff(sp.trace(state*global_h), x) for x in q]))
    require(zero(force_mean-differentiated_mean_scalar), 'force sign disagrees with scalar derivative')
    require(not zero(force_mean), 'force means must expose a global sign error')
    direct_covariance = covariance(state, direct_force)
    per_cell = clean(sum((w*c for w, c in zip(weights, local_covariances)), sp.zeros(3)))
    block_per_cell = {key: clean(sum((weights[s]*local_blocks[s][key] for s in range(2)), sp.zeros(3)))
                      for key in ('RR', 'RE', 'ER', 'EE')}
    direct_blocks = {p+t: covariance(state, direct_parts[p], direct_parts[t])
                     for p in PARTS for t in PARTS}
    require(zero(covariance(state, currents_by_sector[0], currents_by_sector[1])),
            'global sector state must actually be independent')
    require(zero(sum(direct_blocks.values(), sp.zeros(2))-direct_covariance),
            'global direct covariance split reconstruction failed')
    require(any(sp.im(value) != 0 for value in direct_covariance),
            'ordered imaginary force covariance is not exposed')
    require(direct_covariance[0, 1] != 0, 'offdiagonal force covariance is not exposed')
    require(all(not zero(value) for value in direct_blocks.values()),
            'each RR/RE/ER/EE comparison must have a nonzero expectation')
    return dict(ncells=n, weights=weights, polarizations=e, source_map=derived_maps[0],
                direct_force=direct_force, force_mean=force_mean,
                scalar_mean=clean(sp.trace(state*global_h)),
                differentiated_mean_scalar=differentiated_mean_scalar,
                raw_current_mean=mean(state, current_totals)/n,
                raw_current_covariance=per_cell, direct_covariance=direct_covariance,
                current_totals=current_totals, block_per_cell=block_per_cell,
                direct_blocks=direct_blocks, sectors=sector_records)


def compare_actual_consumer(reference):
    """The only actual-consumer import occurs after independent_reference returns."""
    consumer_path = HERE / 'canonical_current_noise.py'
    consumer_binding = binding(consumer_path)
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    consumer = importlib.import_module('canonical_current_noise')
    require(Path(consumer.__file__).resolve() == consumer_path.resolve(), 'wrong consumer module imported')
    e = numeric(reference['polarizations']).real
    n = reference['ncells']
    raw = numeric(reference['raw_current_covariance'])
    expected_m = numeric(reference['source_map'])
    expected_covariance = numeric(reference['direct_covariance'])
    expected_force = np.stack([numeric(value) for value in reference['direct_force']])
    total_j = np.stack([numeric(value) for value in reference['current_totals']])
    actual_m = np.asarray(consumer.current_force_matrix(e), dtype=np.complex128)
    require(actual_m.shape == expected_m.shape, 'wrong source-map shape')
    actual_covariance = np.asarray(consumer.canonical_force_covariance(raw, n, e), dtype=np.complex128)
    checks = {
        'source_vector': comparison(actual_m, expected_m),
        'source_vector_is_real': comparison(actual_m.imag, np.zeros_like(actual_m.real)),
        'differentiated_force_operator': comparison(np.einsum('ma,aij->mij', actual_m, total_j), expected_force),
        'force_mean_vs_scalar_derivative': comparison(n*actual_m@numeric(reference['raw_current_mean']),
                                                     numeric(reference['differentiated_mean_scalar'])),
        'direct_global_fock_covariance': comparison(actual_covariance, expected_covariance),
    }
    block_actual = {}
    for key, raw_block in reference['block_per_cell'].items():
        block_actual[key] = np.asarray(consumer.canonical_force_covariance(numeric(raw_block), n, e),
                                      dtype=np.complex128)
        checks['ordered_'+key] = comparison(block_actual[key], numeric(reference['direct_blocks'][key]))
    checks['ordered_split_reconstruction'] = comparison(sum(block_actual.values()), expected_covariance)
    # Each false boundary variant is evaluated and must actually fail a predicate.
    # Sign is tested on operators/means, because a covariance alone cannot see it.
    false_variants = {
        'wrong_global_force_sign': comparison(np.einsum('ma,aij->mij', -actual_m, total_j), expected_force),
        'missing_ha_ry_factor_two': comparison(actual_covariance/4, expected_covariance),
        'omitted_extensive_n': comparison(consumer.canonical_force_covariance(raw, 1, e), expected_covariance),
        'spurious_n_squared': comparison(consumer.canonical_force_covariance(raw, n*n, e), expected_covariance),
        'renormalized_partial_sector_weights': comparison(
            consumer.canonical_force_covariance(raw*n/len(reference['weights']), n, e),
            expected_covariance),
        'ignored_nontrivial_polarizations': comparison(
            consumer.canonical_force_covariance(raw, n, np.eye(3)[:len(e)]), expected_covariance),
        'dropped_ordered_imaginary_part': comparison(
            consumer.canonical_force_covariance(raw.real, n, e), expected_covariance),
        'dropped_cartesian_offdiagonal_pairs': comparison(
            consumer.canonical_force_covariance(np.diag(np.diag(raw)), n, e), expected_covariance),
    }
    passed = all(item['passes'] for item in checks.values())
    sensitive = all(not item['passes'] and item.get('finite', False)
                    and item.get('max_abs_difference', 0) > ATOL for item in false_variants.values())
    # Bind the very same bytes after execution: a moving source cannot earn a verdict.
    require(binding(consumer_path) == consumer_binding, 'consumer changed during oracle execution')
    return dict(consumer=consumer_binding, tolerances=dict(atol=ATOL, rtol=RTOL),
                actual_source_vector=complex_rows(actual_m),
                actual_force_covariance=complex_rows(actual_covariance),
                actual_ordered_blocks={key: complex_rows(value) for key, value in block_actual.items()},
                checks=checks, false_variant_computed_predicates=false_variants,
                actual_consumer_passes=passed, all_false_variants_rejected=sensitive,
                verdict=passed and sensitive)


def reference_record(reference):
    matrix_keys = ('polarizations', 'source_map', 'force_mean', 'scalar_mean',
                   'differentiated_mean_scalar', 'raw_current_mean',
                   'raw_current_covariance', 'direct_covariance')
    result = {key: exact(reference[key]) for key in matrix_keys}
    result.update(ncells=reference['ncells'], weights=[exact(w) for w in reference['weights']],
                  global_fock_modes=4, global_fock_dimension=16,
                  source_coverage='Two explicit 1/64 sectors in N64; a partial contribution, not full128 coverage.',
                  direct_force_operators=[exact(value) for value in reference['direct_force']],
                  raw_ordered_blocks={key: exact(value) for key, value in reference['block_per_cell'].items()},
                  direct_force_ordered_blocks={key: exact(value) for key, value in reference['direct_blocks'].items()},
                  sectors=reference['sectors'],
                  construction='Exact Slater-mixture Fock traces of -dH_total(q)/dq; H_total=N*sum_s(w_s*H_s(q@E)); J_Ha=-dH_Ry/da/2. No Wick covariance or proposed-map literal supplies the answer.',
                  retained_split='R=P J P, E=J-R with P=diag(1,0), not commuting with either coherent prepared density; every ordered block traced separately.')
    return result


def publish_exclusive(path, record):
    """Durably publish once; never overwrite an existing result or follow a symlink."""
    path = Path(path)
    encoded = (json.dumps(record, sort_keys=True, indent=2, allow_nan=False)+'\n').encode()
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, 'O_NOFOLLOW'):
        flags |= os.O_NOFOLLOW
    fd = os.open(path, flags, 0o644)
    with os.fdopen(fd, 'wb') as stream:
        stream.write(encoded)
        stream.flush()
        os.fsync(stream.fileno())
    directory_fd = os.open(path.resolve().parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='fresh exclusive JSON destination')
    args = parser.parse_args(argv)
    record = dict(schema='P254-A36-exact-canonical-current-force-oracle-v1',
                  created_at=datetime.now(timezone.utc).isoformat(), verdict=False,
                  earned_scope='Only the tiny exact coherent nonstationary product-Gaussian Fock scalar-to-initial-electronic-canonical-generalized-force boundary and the actual floating map on these inputs.',
                  excluded_scope=['No material source/noise recomputation or physical initialization',
                                  'No initial Pi state covariance or additional independent noise',
                                  'No total force-current joint bank or moving kernel',
                                  'No finite-q density, stationarity, FDT, irreversible loss, heat, rates, lifetime, or whole E/C/N closure'],
                  carried_not_replayed=['Existing CAR, UPF, source-restoration, and physical-source qualification proofs; this run makes no new verdict about those proofs'])
    try:
        record['bindings'], plan = authenticate_sources()
        record['original_canonical_boundary'] = {
            key: plan[key] for key in ('original_initial', 'original_boundary', 'original_preparation',
                                      'matched_aggregate', 'matched_archive', 'canonical_construction')}
        record['original_canonical_boundary']['asset_authentication_scope'] = (
            'Bindings carried from the authenticated settled preregistration; this small oracle does not reopen large original assets.')
        reference = independent_reference()
        record['independent_exact_reference'] = reference_record(reference)
        record['execution'] = compare_actual_consumer(reference)
        record['verdict'] = record['execution']['verdict']
        record['status'] = 'passed' if record['verdict'] else 'mismatch_or_insensitive_verifier'
    except Exception as error:
        record['status'] = 'failed'
        record['error'] = dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc())
    publish_exclusive(args.output, record)
    print(json.dumps(dict(output=str(args.output), verdict=record['verdict'], status=record['status']), sort_keys=True))
    return 0 if record['verdict'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
