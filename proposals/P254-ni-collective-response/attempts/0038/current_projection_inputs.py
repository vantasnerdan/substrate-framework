"""Authenticate the unchanged A37 bath and project physical J without building B.

Only the explicit current evaluator invokes physical actions. Cached H/F legs,
original S/D/R and saved E/S^-1/2 are never rebuilt or modified.
"""
from pathlib import Path
import json
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src')] + [str(HERE.parent/name) for name in
    ('0037', '0036', '0035', '0031', '0030', '0029', '0028', '0034')]

import numpy as np
from clamped_scientific_readout import encode
from expanded_bath_inputs import ExpandedBathInputs
from export_frozen_bath_comparator import car_moments, projected_factors
from canonical_current_noise import current_force_matrix
from original_noise_source import _identity, _same_binding
from substrate_framework.material_physical_operators import PhysicalFieldVertices
from substrate_framework.physical_factor_preparation import ActionFactors
from substrate_framework.streamed_material_state import authenticate

CONTRACT = dict(path=str(HERE/'implementation-contract-001.json'), bytes=15306,
    sha256='623ef2c7ebc3fec3b0166783b165458e349abaf6d4d06c4da8bb5c77433f2668')
CURRENT_NAMES = ('current_x', 'current_y', 'current_z')
CURRENT_UNITS = 'electronic_current_hartree'


def load_handoff(args):
    """Authenticate all fixed inputs and the paid basis/temporal lineage."""
    contract = json.loads(authenticate(CONTRACT).read_bytes())
    settled = {Path(row['path']).resolve(): row for row in contract['settled_source_inventory']}
    for binding in settled.values():
        authenticate(binding)

    def bound(path):
        binding = settled.get(Path(path).resolve())
        if binding is None:
            raise ValueError('fixed preregistered input path required: '+str(path))
        return dict(binding), json.loads(Path(binding['path']).read_bytes())

    bindings, data = {}, {}
    for name in ('basis_complete', 'a37_filtered_complete', 'a37_temporal_result', 'current_reference'):
        bindings[name], data[name] = bound(getattr(args, name))
    manifest = data['basis_complete']
    filtered = data['a37_filtered_complete']
    temporal = data['a37_temporal_result']
    reference = data['current_reference']
    if (manifest['schema'] != 'P254-A37-expanded-source0-basis-v1'
            or any(item['status'] != 'complete' for item in data.values())
            or temporal['schema'] != 'P254-A37-expanded-frozen-H-temporal-influence-v1'
            or reference['schema'] != 'P254-A35-original-sector-uniform-current-noise-v1'):
        raise ValueError('complete actual A37 basis/temporal and A35 current inputs required')
    for key, binding in (('basis', manifest['basis']),
                         ('a37_filtered_operators', filtered['filtered_operators'])):
        expected = settled.get(Path(binding['path']).resolve())
        if expected is None or not _same_binding(binding, expected):
            raise ValueError('unregistered bank: '+key)
        bindings[key] = dict(expected)
    if (not _same_binding(filtered['basis'], bindings['basis'])
            or not _same_binding(temporal['basis'], bindings['basis'])
            or not _same_binding(temporal['filtered_operators'], bindings['a37_filtered_operators'])
            or not _same_binding(filtered['inputs']['basis_complete'], bindings['basis_complete'])):
        raise ValueError('paid temporal result and complete basis/bank lineage differ')
    for key in ('original_state_authority', 'original_factor_identity', 'original_pw_basis_identity'):
        if manifest[key] != temporal[key]:
            raise ValueError('paid temporal identity differs: '+key)
    if filtered['normalization'] != temporal['normalization'] or temporal['normalization'] != dict(
            mesh=[4, 4, 4], ncells=64, force_insertion='unweighted F_k/sqrt(N)',
            further_source_weight_applied=False, spin_factor_applied=False,
            final_real_space_inverse_fourier_applied=False):
        raise ValueError('paid ionic insertion must remain F/sqrt(N64), exactly once')
    metadata = manifest['basis_metadata']
    if (metadata['original_rank'] != 346 or metadata['dimension'] != 1038
            or metadata['extra_rank_selected'] != 692 or metadata['extra_rank_requested'] != 692
            or manifest['ncells'] != 64 or manifest['mesh'] != [4, 4, 4]):
        raise ValueError('unchanged original346 plus extra692 bath required')
    field_path = HERE.parent/'0036'/'canonical-field-force-source0-partial-material-result-002.json'
    bindings['a36_field_reference'], data['a36_field_reference'] = bound(field_path)
    return contract, bindings, data


def restore_inputs(manifest, reference, progress):
    """Use the reviewed A37 restoration exactly once, never CachedBathInputs."""
    inputs = ExpandedBathInputs(manifest['action_packet']['path'],
        manifest['mode1_packet']['path'], manifest['filter_packet']['path'], 692, progress)
    context = inputs.context
    for field, actual in (('original_state_authority', context.state_authority),
                          ('original_factor_identity', context.factor_identity),
                          ('original_pw_basis_identity', context.pw_basis_identity)):
        if manifest[field] != actual:
            raise ValueError('original restored identity differs: '+field)
    if (manifest['original_source'] != inputs.original_source
            or manifest['restoration_source'] != inputs.restoration_source
            or not _same_binding(manifest['original_frame'], context.frame)
            or not _same_binding(manifest['original_coordinates'], context.source['coordinates'])
            or manifest['cached_original_actions'] != inputs.cached_action_records):
        raise ValueError('A37 original computation/current restoration/cached factors differ')
    match_current_reference(context, reference)
    return inputs


def match_current_reference(context, reference):
    """Compare physical state/source, not producer-specific whole context dicts."""
    if (reference['original_record'] != 0 or reference['source']['original_record'] != 0
            or reference['state_authority'] != context.state_authority
            or reference['factor_identity'] != context.factor_identity
            or reference['pw_basis_identity'] != context.pw_basis_identity
            or not _same_binding(reference['frame'], context.frame)
            or reference['quadrature_weight'] != context.quadrature_weight
            or context.quadrature_weight != 1/64 or reference['quadrature_applied'] is not False
            or reference['probe_names'] != list(CURRENT_NAMES) or reference['units'] != CURRENT_UNITS):
        raise ValueError('old current reference original state/basis/factor/weight differs')
    physical_fields = ('case', 'model_identity', 'original_model_identity', 'spin', 'source_ik',
        'mesh_index', 'k_fractional', 'operation', 'orbital_conjugation', 'reciprocal_shift',
        'gauge_policy', 'coefficient_case_identity', 'coefficient_item', 'state_convention',
        'measure_convention', 'native_wavefunction', 'selected', 'selection_recipe')
    for field in physical_fields:
        if reference['source'][field] != context.source[field]:
            raise ValueError('old current reference physical source differs: '+field)
    for field in ('source', 'initial_result', 'original_boundary', 'preparation', 'original_raw_d',
                  'source_state', 'source_archive', 'coefficient_manifest', 'coefficient_archive',
                  'saved_background_archive', 'saved_background_receipt', 'population_ledger'):
        if not _same_binding(reference['source'][field], context.source[field]):
            raise ValueError('old current reference source binding differs: '+field)


def load_projection_bank(binding, manifest, inputs):
    """Load only needed paid coordinates/forces, not H or the spectral payload."""
    context = inputs.context
    shapes = dict(extra_basis=(context.psi.shape[0], 692), inverse_sqrt_gram=(346, 346),
        projected_frame=(1038, 346), gram=(346, 346), density_coefficients=(346, 346),
        active_coordinates=context.active_coordinates.shape)
    for name in inputs.forces:
        shapes['force_'+name] = (1038, 1038)
        shapes['adjoint_force_'+name] = (1038, 1038)
    with np.load(binding['path'], allow_pickle=False) as archive:
        bank = {key: archive[key] for key in shapes}
    for key, shape in shapes.items():
        if bank[key].shape != shape or not np.isfinite(bank[key]).all():
            raise ValueError('saved finite-bath coordinate shape/data differs: '+key)
        bank[key].setflags(write=False)
    for key, actual in (('gram', context.gram), ('density_coefficients', context.density_coefficients),
                        ('active_coordinates', context.active_coordinates)):
        if not np.array_equal(bank[key], actual):
            raise ValueError('saved original S/D/R changed: '+key)
    if manifest['basis_identity'] != 'A37 finite-bath coordinate basis '+binding['sha256']:
        raise ValueError('original saved finite-bath identity differs')
    return bank


def field_map(field_reference, inputs):
    """Use the actual A36 polarizations and existing M=2e implementation."""
    normalization = field_reference['normalization']
    source = inputs.original_source
    if (field_reference['original_records'] != [0]
            or field_reference['joint_model_identity'] != source['model_identity']
            or normalization['ncells'] != inputs.ncells
            or normalization['quadrature_weight'] != inputs.context.quadrature_weight
            or not _same_binding(field_reference['original_inputs']['initial_result'], source['initial_result'])
            or not _same_binding(field_reference['original_state']['preparation'], source['preparation'])
            or not _same_binding(field_reference['original_state']['boundary'], source['original_boundary'])):
        raise ValueError('actual A36 field map original state/source normalization differs')
    polarizations = np.asarray(field_reference['field_normalization']['polarizations'], dtype=float)
    matrix = current_force_matrix(polarizations)
    if not np.array_equal(matrix, normalization['current_force_matrix']):
        raise ValueError('actual A36 saved M differs from existing canonical map')
    amplitude = np.sqrt(inputs.ncells*inputs.context.quadrature_weight)*matrix
    return polarizations.tolist(), matrix.tolist(), amplitude.tolist()


def current_specification(context, axis):
    direction = np.zeros(3)
    direction[axis] = 1.
    return dict(kind='physical-uniform-electronic-current', axis='xyz'[axis],
        field_direction=direction.tolist(), external_q=[0, 0, 0], time_inverse_ry=0,
        units=CURRENT_UNITS, current_scale=-0.5,
        field_api='PhysicalFieldVertices.prepare(real_direction).apply',
        adjoint_api='PhysicalFieldAction.adjoint_apply delegates apply',
        adjoint_reuse='real direction and checked original Hermitian UPF D; J_dagger=J',
        basis_identity=context.pw_basis_identity, source=context.source['source'],
        physical_sources=encode(list(context._coupling.factors.sources)), quadrature_applied=False)


def project_current(psi, inverse, extra, original_action, extra_action):
    """Four raw blocks of B^dagger J B, without allocating B or copied Q0."""
    matrix = np.empty((1038, 1038), dtype=np.complex128)
    matrix[:346, :346] = inverse @ (psi.conj().T @ original_action) @ inverse
    matrix[346:, :346] = (extra.conj().T @ original_action) @ inverse
    matrix[:346, 346:] = inverse @ (psi.conj().T @ extra_action)
    matrix[346:, 346:] = extra.conj().T @ extra_action
    matrix.setflags(write=False)
    return matrix


def evaluate_currents(inputs, bank, reference, progress, calls):
    """Exactly three intentional346 re-evaluations and three new692 evaluations."""
    context = inputs.context
    vertices = PhysicalFieldVertices(context._hamiltonian.frozen.h0, projector_block_size=64)
    matrices, factors, specifications = {}, {}, {}
    for axis, name in enumerate(CURRENT_NAMES):
        specification = current_specification(context, axis)
        identity = _identity(specification)
        old = reference['covariance_metadata']['operators'][name]
        if (specification != reference['source']['operators'][name]
                or old['operator_identity'] != identity or old['adjoint_identity'] != identity):
            raise ValueError('physical current specification differs from A35: '+name)
        t = perf_counter()
        prepared = vertices.prepare(specification['field_direction'])
        progress(dict(stage='current_direction_factory', name=name, status='complete',
                      elapsed_seconds=perf_counter()-t, physical_actions=0))
        actions = []
        for role, vectors in (('original346_intentional_reevaluation', context.psi),
                              ('extra692_new_evaluation', bank['extra_basis'])):
            record = dict(name=name, role=role, columns=vectors.shape[1], status='started',
                operator_identity=identity, projector_block_size=64,
                internal_column_blocks=(vectors.shape[1]+63)//64,
                independent_adjoint_action=False, qualified_adjoint_reuses=1)
            calls.append(record)
            progress(dict(stage='physical_current_action', **record))
            t = perf_counter()
            try:
                action = prepared.apply(vectors)
                action *= -0.5
                action.setflags(write=False)
                record['status'] = 'complete'
                actions.append(action)
            except BaseException as error:
                record.update(status='failed', error=repr(error))
                raise
            finally:
                record['elapsed_seconds'] = perf_counter()-t
                progress(dict(stage='physical_current_action', **record))
        original, extra_action = actions
        matrices[name] = project_current(context.psi, bank['inverse_sqrt_gram'],
            bank['extra_basis'], original, extra_action)
        factors[name] = ActionFactors(original, original, identity, identity, CURRENT_UNITS,
            context.pw_basis_identity, context.state_authority)
        specifications[name] = specification
        # Keep original346 legs only until the full-PW CAR readout; never persist them.
        del actions, extra_action, action, prepared
    return matrices, factors, specifications


def initial_comparison(inputs, bank, matrices, currents, basis_identity, reference):
    """All25 directed CAR pairs with mixed units, raw deficits and old9 comparison."""
    context = inputs.context
    full_factors = dict(inputs.forces, **currents)
    operators = {name: bank['force_'+name] for name in inputs.forces}
    adjoints = {name: bank['adjoint_force_'+name] for name in inputs.forces}
    operators.update(matrices)
    adjoints.update(matrices)  # qualified true J^dagger=J, not matrix.conj().T
    units = {name: factor.units for name, factor in full_factors.items()}
    projected = {}
    for name in operators:
        projected.update(projected_factors({name: operators[name]}, {name: adjoints[name]},
            bank['projected_frame'], basis_identity, context.state_authority, units[name]))
    full = car_moments(context.psi, context.gram, context.density_coefficients,
        full_factors, context.factor_identity, context.pw_basis_identity, context.state_authority)
    finite = car_moments(bank['projected_frame'], bank['gram'], bank['density_coefficients'],
        projected, context.factor_identity, basis_identity, context.state_authority)
    full_map = {(row['left'], row['right']): row['value'] for row in full['covariance']}
    finite_map = {(row['left'], row['right']): row['value'] for row in finite['covariance']}
    old_map = {(row['left'], row['right']): complex(row['full']['real'], row['full']['imag'])
               for row in reference['covariance']}
    expected = {(a, b) for a in full_factors for b in full_factors}
    if (set(full_map) != expected or set(finite_map) != expected
            or set(old_map) != {(a, b) for a in CURRENT_NAMES for b in CURRENT_NAMES}):
        raise ValueError('all25 initial pairs and all9 old current pairs required')
    return dict(full_pw=full, finite_bath=finite, probe_order=list(full_factors), probe_units=units,
        covariance_units='product of the two declared probe units; mixed force/current pairs',
        probe_metadata={name: {field: getattr(factor, field) for field in (
            'operator_identity', 'adjoint_identity', 'units', 'pw_basis_identity', 'source_authority')}
            for name, factor in full_factors.items()},
        finite_bath_probe_metadata={name: {field: getattr(factor, field) for field in (
            'operator_identity', 'adjoint_identity', 'units', 'pw_basis_identity', 'source_authority')}
            for name, factor in projected.items()},
        ordered_projection_deficit=[dict(left=a, right=b, full_pw_minus_finite_bath=full_map[a,b]-finite_map[a,b],
            units_left=units[a], units_right=units[b]) for a,b in full_map],
        old_current_reference_comparison=[dict(left=a, right=b, original_A35=old_map[a,b],
            reevaluated_full_pw=full_map[a,b], raw_difference=full_map[a,b]-old_map[a,b]) for a,b in old_map],
        raw_current_matrix_defects={name: dict(hermiticity_fro=float(np.linalg.norm(matrix-matrix.conj().T)),
            qualified_supplied_adjoint_minus_matrix_adjoint_fro=float(np.linalg.norm(matrix-matrix.conj().T)))
            for name, matrix in matrices.items()},
        scope='Initial unweighted source0 F/J CAR covariances and finite-bath current projection deficit only. '
            'No Hermitization, clipping, averaging, tolerance acceptance, late-time error bound, joint temporal result or loss claim.')
