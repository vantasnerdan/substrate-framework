"""Finite coherent-publication and completed-noise join regressions.

References come from independent discrete-forward/dense physical-density
calculations in the source-local oracle, not from the publication wrapper.
"""
from copy import deepcopy
from pathlib import Path
from hashlib import sha256
from time import time
from types import SimpleNamespace
import os

import numpy as np

from substrate_framework.material_control_reference import bound_json
from substrate_framework.streamed_material_state import file_binding
from joint_spatial_dual_coherent import build_readout

HERE = Path(__file__).resolve().parent


def _without_noise(entry, seed, proof, snapshots):
    pair = bound_json(entry['pair']['result'])
    paths = [Path(entry['pair']['result']['path']),
             Path(pair['covariance']['context_pair_binding']['path'])]
    original = {path: path.read_bytes() for path in paths}
    try:
        for path in paths:
            path.unlink()
        return tuple(build_readout(entry['pair']['owner'], pair['covariance'][side+'_context'],
                                   seed, proof, snapshots=snapshots) for side in ('left', 'right'))
    finally:
        for path, raw in original.items():
            path.write_bytes(raw)


def _reject_incomplete(entry, seed, proof, save, reject, snapshots):
    pair = bound_json(entry['pair']['result'])
    context = bound_json(pair['covariance']['left_context'])
    context['complete'] = False
    partial = save('incomplete-one-sided-pullback.json', context)
    reject('source_local_dual_incomplete_pullback',
           lambda: build_readout(entry['pair']['owner'], partial, seed, proof, snapshots=snapshots))
    owner = bound_json(entry['pair']['owner'])
    owner['request']['selectors']['axis'] = 'z'
    foreign = save('foreign-one-sided-owner.json', owner)
    reject('source_local_dual_foreign_owner',
           lambda: build_readout(foreign, pair['covariance']['left_context'], seed, proof, snapshots=snapshots))


def _reject_changed(entries, seed_binding, consume, save, reject):
    seed = bound_json(seed_binding)
    for name in ('normalization', 'context', 'time', 'probe', 'seed'):
        changed = deepcopy(entries)
        native = bound_json(changed[-1]['coherent']['result'])
        if name == 'normalization':
            native['readout']['real'] *= 3**.5
            native['readout']['imag'] *= 3**.5
        elif name == 'context':
            pair = bound_json(changed[-1]['pair']['result'])
            native['context'] = pair['covariance']['right_context']
        elif name == 'time':
            native['time_ry_inverse'] += .001
        elif name == 'probe':
            native['observable'] = 'original-H1y-total-force'
        else:
            wrong = save('different-original-seed.json', dict(seed, amplitude=1.9))
            native['seed'] = changed[-1]['coherent']['seed'] = wrong
        changed[-1]['coherent']['result'] = save('wrong-dual-'+name+'.json', native)
        reject('source_local_dual_wrong_'+name, lambda: consume(changed))


def exercise_dual_join(entries, seed, discrete_checks, consume, reference, save, close, reject,
                       *, snapshots=(), one_sided=()):
    sources = [file_binding(HERE/name) for name in
               ('joint_spatial_dual_coherent.py', 'joint_spatial_dual_authority.py')]
    proof = save('computed-dual-oracle.json', dict(schema='P254-joint-spatial-readout-oracle-v1',
        status='passed', source_bindings=sources, checks=discrete_checks,
        scope='executed finite discrete-forward/dual numerical comparisons'))
    dual_entries = []
    for index, entry in enumerate(entries):
        native, right = _without_noise(entry, seed, proof, snapshots)
        for side, record in (('left', native), ('right', right)):
            value = complex(record['readout']['real'], record['readout']['imag'])
            close('source_local_dual_without_noise_%s_q%d' % (side, index),
                  value, reference[side][index], 3e-8)
        dual_entries.append(dict(coherent=dict(result=save('dual-readout-%d.json' % index, native),
                                               seed=seed), pair=entry['pair']))
    _, forces, ab, ba = consume(dual_entries)
    close('source_local_native_dual_join_coherent', forces, reference['left'], 3e-8)
    close('source_local_native_dual_join_AB', ab, reference['AB'], 0.)
    close('source_local_native_dual_join_BA', ba, reference['BA'], 0.)
    reject('source_local_dual_missing_q', lambda: consume(dual_entries[:-1]))
    _reject_incomplete(entries[0], seed, proof, save, reject, snapshots)
    original_pair = bound_json(entries[0]['pair']['result'])
    reject('source_local_historical_sources_require_declared_snapshots',
        lambda: build_readout(entries[0]['pair']['owner'],
            original_pair['covariance']['left_context'], seed, proof))
    _reject_changed(dual_entries, seed, consume, save, reject)
    for fixture in one_sided:
        _material_pullback(fixture, seed, proof, save, close, reject, snapshots)


def _material_pullback(fixture, seed, proof, save, close, reject, snapshots):
    """Execute the material adapter on genuine finite full-rectangular kernels."""
    from run_source_local_q_material import execute_pullback, scientific_request, validate_mode
    from source_local_material_owner import context_bindings, material_owner, material_budget, CAP_NAMES
    from source_local_material_runtime import MaterialRuntime, _selectors, executor_directory
    from run_joint_partitioned_q_oracle import DiskFiniteHistory
    from substrate_framework.partitioned_q_response import PartitionedQResponse
    from substrate_framework.source_local_q_response import SourceLocalQResponse
    from joint_spatial_provenance import canonical
    from joint_spatial_dual_coherent import validate_readout
    from run_source_local_q_oracle import exact_dual
    from run_joint_response_oracle import pack_legs

    system, response, preparation, probe, grid, directory, provenance, expected, transverse_expected = fixture
    directory.mkdir()
    p = deepcopy(provenance)
    runtime_request = p['request']
    selectors = dict(runtime_request['selectors'], mode='pullback', edge_index=0)
    del selectors['second_time']
    args = SimpleNamespace(**selectors, second_time=None, cpus=[min(os.sched_getaffinity(0))],
        resume=False, frame_admission=[], left_context=None, left_context_sha256=None,
        right_context=None, right_context_sha256=None, pullback_context=None, pullback_context_sha256=None,
        max_actions=100000, max_seconds=900., max_output_bytes=512*1024**2,
        max_memory_bytes=512*1024**2, worker_memory_bytes=256*1024**2)
    for name, runtime_name in (('initialization', 'initialization'), ('history', 'baseline'),
            ('physical_qualification', 'physical_qualification'), ('source_local_oracle', 'source_local_oracle')):
        binding = runtime_request[runtime_name]
        setattr(args, name, Path(binding['path']))
        setattr(args, name+'_sha256', binding['sha256'])
    validate_mode(args)
    caps = {name: getattr(args, name) for name in CAP_NAMES}
    runtime_request['selectors'] = _selectors(args)
    runtime_request['times'] = list(grid)
    del runtime_request['second_times']
    sources = p['source_bindings']
    request = scientific_request(args, sources, caps)
    owner, owner_binding, identity = material_owner(directory, request, resume=False, started_unix=time()-1)
    propagation = sha256(canonical(dict(model=p['model_identity'], q=response.q_fractional.tolist(),
        dynamics='screened-joint', sources=sources, observable_metadata=runtime_request)).encode()).hexdigest()

    def runtime_for(budget):
        executor = PartitionedQResponse(response, DiskFiniteHistory(system, response, budget),
            kernel_adjoint=lambda n, **kw: system.k*n, grid_measure=1., budget=budget,
            directory=executor_directory(directory/'scientific-work'), identity=propagation,
            cpus=args.cpus, work_directory=directory/'scientific-work'/'physical-work')
        executor.provenance = p
        return MaterialRuntime(SourceLocalQResponse(executor, directory/'scientific-work'/'pullback'),
            preparation, probe, np.asarray(grid), None, p, dict(mode='pullback'))

    budget = material_budget(directory, owner, caps)
    runtime = runtime_for(budget)
    result = dict(covariance_executed=False)
    execute_pullback(runtime, {}, result)
    context = bound_json(result['pullback_context'])
    label = 'material_one_sided_q'+str(response.q_fractional[-1])
    intervals = len(grid)-1
    close(label+'_only_one_pullback_no_final_electronic_replay',
          budget.snapshot()['hamiltonian_actions'],
          2*len(response.edges)*(8*intervals**2-2*intervals), 0.)
    dual = runtime.executor.restore(probe, grid, result['pullback_context'])
    close(label+'_complete_rectangular_dual', pack_legs(dual),
          pack_legs(exact_dual(response, system, probe, grid)), 3e-11)
    dual.release()
    native = build_readout(owner_binding, result['pullback_context'], seed, proof, snapshots=snapshots)
    close(label+'_original_seed_force', complex(native['readout']['real'],
          native['readout']['imag']), expected, 3e-8)
    transverse = dict(bound_json(seed), axis='y')
    transverse_binding = save(label+'-transverse-seed.json', transverse)
    independent = build_readout(owner_binding, result['pullback_context'], transverse_binding,
                               proof, snapshots=snapshots)
    close(label+'_independent_seed_selector', complex(independent['readout']['real'],
          independent['readout']['imag']), transverse_expected, 3e-8)
    if (set(result) != {'covariance_executed', 'pullback_context', 'pullback_executed', 'scientific_scope'}
            or result['covariance_executed'] or not result['pullback_executed']
            or any(path.name in ('pair', 'contraction-work') for path in directory.rglob('*'))):
        raise AssertionError('one-sided adapter performed or forged covariance work')
    original_owner = Path(owner_binding['path']).read_bytes()
    original_context = Path(result['pullback_context']['path']).read_bytes()
    args.resume = True
    args.pullback_context = Path(result['pullback_context']['path'])
    args.pullback_context_sha256 = result['pullback_context']['sha256']
    restored_owner, restored_binding, restored_identity = material_owner(directory,
        scientific_request(args, sources, caps), resume=True, started_unix=time()+50)
    if (restored_owner != owner or restored_binding != owner_binding or restored_identity != identity):
        raise AssertionError('one-sided resume refreshed the original owner/start/caps')
    resumed_budget = material_budget(directory, restored_owner, caps)
    resumed = runtime_for(resumed_budget)
    execute_pullback(resumed, context_bindings(args, directory), result)
    if (Path(owner_binding['path']).read_bytes() != original_owner
            or Path(result['pullback_context']['path']).read_bytes() != original_context):
        raise AssertionError('one-sided resume replaced original owner or context bytes')
    close(label+'_resume_preserves_original_clock', resumed_budget.started_at, owner['started_unix'], 0.)
    close(label+'_resume_preserves_actions', resumed_budget.snapshot()['hamiltonian_actions'],
          budget.snapshot()['hamiltonian_actions'], 0.)
    validate_readout(native, bound_json(seed), snapshots=snapshots)
    for name, value in (('time', args.time+.001), ('axis', 'z')):
        changed = deepcopy(args)
        setattr(changed, name, value)
        reject(label+'_changed_'+name+'_owner', lambda: material_owner(directory,
            scientific_request(changed, sources, caps), resume=True, started_unix=time()))
    changed_caps = dict(caps, max_actions=caps['max_actions']+1)
    reject(label+'_changed_caps_owner', lambda: material_owner(directory,
        scientific_request(args, sources, changed_caps), resume=True, started_unix=time()))
    reject(label+'_fresh_context_owner', lambda: context_bindings(dict_args(args, resume=False), directory))
    reject(label+'_foreign_context_namespace', lambda: context_bindings(args, directory/'unrelated'))
    reject(label+'_pair_context_selector', lambda: context_bindings(dict_args(args,
        pullback_context=None, pullback_context_sha256=None, left_context=args.pullback_context,
        left_context_sha256=args.pullback_context_sha256), directory))
    reject(label+'_meaningless_second_time', lambda: validate_mode(dict_args(args, second_time=.004)))
    reject(label+'_pair_still_requires_second_time', lambda: validate_mode(dict_args(args, mode='pair')))
    foreign = save(label+'-foreign-owner.json', owner)
    reject(label+'_foreign_context_consumer', lambda: build_readout(foreign,
        result['pullback_context'], seed, proof, snapshots=snapshots))


def dict_args(args, **changes):
    return type(args)(**dict(vars(args), **changes))
