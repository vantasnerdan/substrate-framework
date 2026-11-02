"""Computed finite source-local publications joined to dense coherent readouts.

Temporary native headers are explicitly finite fixtures, never material results.
Neither ordered covariance nor the coherent scalar is supplied by the parser.
"""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import os
from tempfile import TemporaryDirectory

import numpy as np
from scipy.linalg import block_diag

from run_joint_response_oracle import FiniteSystem, DenseTangent
from run_joint_partitioned_q_oracle import DiskFiniteHistory
from substrate_framework.joint_response_observables import total_force_probe, local_momentum_kick
from substrate_framework.joint_ordered_sources import JointProductPreparation
from substrate_framework.partitioned_q_response import PartitionedQResponse
from substrate_framework.partitioned_q_storage import QBudget
from substrate_framework.partitioned_joint_slots import digest
from substrate_framework.source_local_q_covariance import SourceLocalTwoSidedCovariance
from substrate_framework.streamed_material_state import file_binding
from substrate_framework.joint_spatial_readout import reconstruct
from substrate_framework.material_control_reference import bound_json
from joint_spatial_provenance import canonical, same_binding
from joint_spatial_receipts import read_entries, validate_records
from joint_spatial_committed_oracle import _native_fixture
from joint_spatial_dual_coherent import canonical_impulse, SEED_SCHEMA
from joint_spatial_dual_coherent_oracle import exercise_dual_join
from run_source_local_q_oracle import forward_matrix, terminal
from run_joint_response_oracle import pack_legs

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def exercise_source_local(*, encode, close, reject):
    with TemporaryDirectory(prefix='p254-spatial-source-local-') as temporary:
        directory = Path(temporary)

        def save(name, value):
            path = directory/name
            path.write_text(json.dumps(encode(value), allow_nan=False)+'\n')
            return file_binding(path)

        physical_source = file_binding(ROOT/'src/substrate_framework/joint_bloch_response.py')
        worker = file_binding(HERE/'joint_spatial_source_local_oracle.py')
        pair_worker = file_binding(ROOT/'src/substrate_framework/source_local_q_covariance.py')
        qualification = file_binding(HERE/'partitioned-readout-current-qualification-002.json')
        numerical = file_binding(HERE/'source-local-q-oracle-006.json')
        readout_proof = file_binding(HERE/'joint-q-readout-oracle-004.json')
        proof_sources = (bound_json(qualification)['source_bindings']
                         +bound_json(numerical)['sources']+[physical_source, pair_worker])
        source_map = {}
        for source in proof_sources:
            previous = source_map.setdefault(source['path'], source)
            if not same_binding(previous, source):
                raise ValueError('retained physical/source-local proof source disagreement')
        sources = list(source_map.values())
        baseline = save('baseline.json', dict(finite_fixture_only=True))
        aggregate = save('aggregate.json', dict(finite_fixture_only=True))
        specification = dict(scalar='finite-moving-physical-density', aggregate=aggregate,
                             implementation=[physical_source])
        model = sha256(canonical(specification).encode()).hexdigest()
        system = FiniteSystem()
        electronic = save('original-electronic.json', dict(frames=system.frames, D=system.ds))
        preparation_binding = save('original-preparation.json', dict(ionic=system.wi, field=system.wf))
        joint = save('original-boundary.json', dict(
            slots=[dict(identity='1:0,0,%d' % i) for i in range(3)],
            time_ry_inverse=0., model_identity=model))
        saved = dict(case='finite-source-local', model_identity=model,
                     specification=specification, joint_boundary=joint,
                     original_boundary=electronic, preparation=preparation_binding)
        initialization = save('initialization.json', dict(schema='P254-matched-joint-case-v1',
            status='complete', saved=saved, mode='initialize', case=saved['case']))
        q = np.array([[0., 0., k/3] for k in range(3)])
        mesh, cell = np.array([1, 1, 3]), [0, 0, 1]
        grids = ([0., .003, .009], [0., .004])
        entries, dense_left, covariance_blocks, mean_blocks = [], [], {o: [] for o in ('AB', 'BA')}, []
        seed = dict(schema=SEED_SCHEMA, initial_perturbation='ionic-momentum-only',
            case=saved['case'], model_identity=model, initialization=initialization,
            atom=1, axis='x', amplitude=1.7, cell=cell)
        seed_binding = save('original-seed.json', seed)
        discrete_checks = {}
        for index, point in enumerate(q):
            response = system.response(index == 0)
            # A complete finite three-character family. The non-Gamma members
            # deliberately retain complex defects rather than impose a repair.
            response.q_fractional = point.copy()
            dense = DenseTangent(system, response)
            probe = total_force_probe(response, [1., 0., 0.], identity='original-H1x-total-force')
            coefficient = dense.coefficient(probe)
            tl, tr = dense.transition(grids[0][-1]), dense.transition(grids[1][-1])
            kick = local_momentum_kick(response, cell, [1., 0., 0.], 1.7)
            initial = dense.physical(system.stage(0), kick)
            mean = tl@initial
            force = np.vdot(coefficient, mean)
            left_row = coefficient.conj()@tl
            right_row = coefficient.conj()@tr
            dense_left.append((left_row, right_row))
            mean_blocks.append(initial)
            for ordering in ('AB', 'BA'):
                w = dense.covariance(ordering)
                parts = []
                for start, stop in ((0, dense.ne), (dense.ne, dense.ne+6), (dense.ne+6, dense.size)):
                    selected = np.zeros_like(w)
                    selected[start:stop, start:stop] = w[start:stop, start:stop]
                    parts.append(selected)
                covariance_blocks[ordering].append(parts)
            selectors = dict(mode='pair', case=saved['case'], q=point.tolist(), time=grids[0][-1],
                second_time=grids[1][-1], subdivisions=1, atom=1, axis='x', observable='total-force',
                field_axis=None, dynamics='screened-joint', block_size=8, edge_index=0)
            runtime = dict(initialization=initialization, baseline=baseline,
                physical_qualification=qualification, source_local_oracle=numerical,
                selectors={k: v for k, v in selectors.items() if k not in ('mode', 'edge_index')},
                times=grids[0], second_times=grids[1], frame_admissions=[], observable_metadata=None)
            provenance = dict(case=saved['case'], model_identity=model, baseline=baseline,
                initialization=initialization, physical_qualification=qualification,
                source_local_numerical_proof=numerical, request=runtime, matched_aggregate=aggregate,
                scalar=specification['scalar'], observable_metadata=None, observable=probe.identity,
                all_explicit_edges=3, fixed_original_D=True, source_bindings=sources)
            identity = sha256(canonical(dict(model=model, q=point.tolist(), dynamics='screened-joint',
                sources=sources, observable_metadata=runtime)).encode()).hexdigest()
            work = directory/('q-%d' % index)
            budget = QBudget(max_actions=100000, max_seconds=900, max_bytes=512*1024**2,
                             max_memory_bytes=512*1024**2, worker_memory_bytes=256*1024**2)
            history = DiskFiniteHistory(system, response, budget)
            executor = PartitionedQResponse(response, history, kernel_adjoint=lambda n, **kw: system.k*n,
                grid_measure=1., budget=budget, directory=work/'executor', identity=identity,
                cpus=[min(os.sched_getaffinity(0))], work_directory=work/'physical-work')
            executor.provenance = provenance
            preparation = JointProductPreparation(response, system.stage(0), ionic_ordered=system.wi,
                ionic_minus_ordered=system.wi, field_ordered=system.wf if index == 0 else np.zeros((0, 0)),
                electronic_identity='finite-original-D', ionic_identity='finite-original-W',
                field_identity='finite-original-field', product_identity='original-matched-product-'+model)
            actual = SourceLocalTwoSidedCovariance(preparation, executor).pair(probe, probe, *grids, work/'pair')
            context = bound_json(actual['left_context'])
            dual_force = canonical_impulse(context, point, 3, seed)
            discrete = pack_legs(kick)
            for a, b in zip(grids[0][:-1], grids[0][1:]):
                discrete = forward_matrix(response, system, a, b-a)@discrete
            discrete_force = np.vdot(pack_legs(terminal(response, system, probe, grids[0][-1])), discrete)
            label = 'source_local_dual_discrete_forward_q%d' % index
            close(label, dual_force, discrete_force, 3e-11)
            close('source_local_dual_dense_forward_q%d' % index, dual_force, force, 3e-8)
            discrete_checks[label] = dict(passed=True, max_absolute_error=float(abs(dual_force-discrete_force)))
            transverse_seed = dict(seed, axis='y')
            transverse = local_momentum_kick(response, cell, [0., 1., 0.], 1.7)
            transverse_expected = np.vdot(coefficient, tl@dense.physical(system.stage(0), transverse))
            close('source_local_dual_distinct_seed_probe_axis_q%d' % index,
                  canonical_impulse(context, point, 3, transverse_seed), transverse_expected, 3e-8)
            if index == 1:
                impulse = 1.7*np.exp(-2j*np.pi*point@cell)/np.sqrt(3)
                raw = complex(*context['small']['p'][0])
                if abs(raw*impulse-dual_force) < 1e-9 or abs(raw.conjugate()*impulse.conjugate()-dual_force) < 1e-9:
                    raise AssertionError('dual conjugation/phase witness is not discriminating')
            for ordering in ('AB', 'BA'):
                expected = [left_row@w@right_row.conj() for w in covariance_blocks[ordering][-1]]
                close('source_local_computed_q%d_%s_origins' % (index, ordering),
                      [actual[ordering][o] for o in ('electronic', 'ionic', 'field')], expected, 3e-8)
            explicit = {k: dict(path=b['path'], sha256=b['sha256']) for k, b in (
                ('initialization', initialization), ('history', baseline),
                ('physical_qualification', qualification), ('source_local_oracle', numerical))}
            owned_request = dict(schema='P254-source-local-q-material-v1-request', selectors=selectors,
                original_time_origin=0., source_bindings=sources, explicit_inputs=explicit)
            owner = save('owner-%d.json' % index, dict(schema='P254-source-local-q-material-v1-owner', request=owned_request))
            result = dict(schema='P254-source-local-q-material-v1-result', status='complete', mode='pair',
                covariance_executed=True, required_q_for_whole_local=3, owner=owner,
                material_identity=digest(dict(schema='P254-source-local-q-material-v1', owner=owner)),
                source_bindings=sources, explicit_inputs=explicit, provenance=provenance,
                time_grid=grids[0], second_time_grid=grids[1], covariance=actual)
            request = dict(identity='finite-independent-dense-forward-'+str(index), mode='forward',
                case=saved['case'], initialization=initialization, baseline=baseline, q=point.tolist(),
                times=grids[0], second_times=grids[1], atom=1, axis='x', amplitude=1.7, cell=cell,
                subdivisions=1, block_size=8, checkpoint='None', checkpoint_sha256=None)
            request_binding = save('historical-forward-%d.json' % index, request)
            forward = dict(case=saved['case'], model_identity=model, finite_specification=specification,
                scalar=specification['scalar'], q=point.tolist(), dynamics='screened-joint',
                source_bindings=[physical_source], initialization=initialization, baseline=baseline,
                matched_aggregate=aggregate, observable=probe.identity, all_explicit_edges=3,
                coherent=dict(readout=force))
            native = _native_fixture(save, forward, request_binding, worker, index,
                                     physical=qualification, scalar=readout_proof)
            entries.append(dict(forward=dict(result=save('native-%d.json' % index, native), request=request_binding),
                                pair=dict(result=save('pair-%d.json' % index, result), owner=owner)))

        def consume(selected):
            records, requests, pairs, owners = read_entries(selected)
            return validate_records(records, requests, saved, initialization, mesh,
                                    pair_records=pairs, pair_requests=owners)

        common, forces, ab, ba = consume(entries)
        exercise_dual_join(entries, seed_binding, discrete_checks, consume,
            dict(left=forces, right=np.array([rows[1]@mean for rows, mean in
                 zip(dense_left, mean_blocks)]), AB=ab, BA=ba), save, close, reject)
        cells = np.array([[0, 0, r] for r in range(3)])
        # Dense real-space rows acting on the direct sum of independently
        # propagated physical-density/canonical blocks. No production Fourier
        # or scalar-contraction function supplies these reference matrices.
        character = np.exp(2j*np.pi*np.outer(np.arange(3), np.arange(3))/3)/np.sqrt(3)
        l = np.array([np.concatenate([character[r, k]*dense_left[k][0] for k in range(3)]) for r in range(3)])
        r = np.array([np.concatenate([character[x, k]*dense_left[k][1] for k in range(3)]) for x in range(3)])
        expected_force = l@np.concatenate(mean_blocks)
        expected = {}
        for ordering in ('AB', 'BA'):
            real = [l@block_diag(*[covariance_blocks[ordering][k][o] for k in range(3)])@r.conj().T for o in range(3)]
            expected[ordering] = np.column_stack([x[:, 1] for x in real]+[sum(real)[:, 1]])
        result = reconstruct(q, forces, ab, ba, mesh=mesh, cell_bohr=np.eye(3),
                             kick_cell=cell, reference_cell=cell,
                             probe_position_bohr=np.zeros(3), seed_position_bohr=np.zeros(3))
        close('source_local_join_dense_real_space_coherent', result['coherent_force_ry_bohr'], expected_force, 3e-8)
        for ordering in ('AB', 'BA'):
            close('source_local_join_dense_real_space_'+ordering,
                  result['connected_force_pairs_ry2_bohr2'][ordering], expected[ordering], 3e-8)
        if np.linalg.norm(expected['AB'][:, 2]) == 0:
            raise AssertionError('source-local field origin lacks numerical signal')
        reject('source_local_missing_q', lambda: consume(entries[:-1]))
        for name in ('baseline', 'probe', 'normalization', 'source', 'publication'):
            changed = deepcopy(entries)
            pair = bound_json(changed[-1]['pair']['result'])
            if name == 'baseline':
                pair['provenance']['baseline'] = aggregate
            elif name == 'probe':
                pair['provenance']['observable'] = 'original-H1y-total-force'
            elif name == 'normalization':
                pair['covariance']['convention'] = 'extra 1/N already applied'
            elif name == 'source':
                pair['source_bindings'] = [worker]
            else:
                pair['covariance']['AB']['total']['real'] += .01
            changed[-1]['pair']['result'] = save('mutated-'+name+'.json', pair)
            reject('source_local_wrong_'+name, lambda: consume(changed))
        # Numerical witnesses expose extra sqrtN / N even if someone preserves
        # schema spelling while scaling otherwise self-consistent values.
        for name, wrong, expected_value in (
                ('coherent_sqrtN', expected_force/np.sqrt(3), expected_force),
                ('covariance_1overN', expected['AB']/3, expected['AB'])):
            if np.max(np.abs(wrong-expected_value)) <= 1e-6:
                raise AssertionError('insensitive source-local '+name+' witness')
        # Mutating a source's bytes, rather than its path or source list, must
        # fail before values are admitted; the original native file is untouched.
        source_copy = save('source-copy.json', dict(finite_original_source=True))
        from joint_spatial_provenance import authenticate
        Path(source_copy['path']).write_text('changed source bytes\n')
        reject('source_local_changed_source_bytes', lambda: authenticate(source_copy))
