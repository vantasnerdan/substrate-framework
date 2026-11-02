"""Finite coherent-publication and completed-noise join regressions.

References come from independent discrete-forward/dense physical-density
calculations in the source-local oracle, not from the publication wrapper.
"""
from copy import deepcopy
from pathlib import Path

from substrate_framework.material_control_reference import bound_json
from substrate_framework.streamed_material_state import file_binding
from joint_spatial_dual_coherent import build_readout

HERE = Path(__file__).resolve().parent


def _without_noise(entry, seed, proof):
    pair = bound_json(entry['pair']['result'])
    paths = [Path(entry['pair']['result']['path']),
             Path(pair['covariance']['context_pair_binding']['path'])]
    original = {path: path.read_bytes() for path in paths}
    try:
        for path in paths:
            path.unlink()
        return tuple(build_readout(entry['pair']['owner'], pair['covariance'][side+'_context'],
                                   seed, proof) for side in ('left', 'right'))
    finally:
        for path, raw in original.items():
            path.write_bytes(raw)


def _reject_incomplete(entry, seed, proof, save, reject):
    pair = bound_json(entry['pair']['result'])
    context = bound_json(pair['covariance']['left_context'])
    context['complete'] = False
    partial = save('incomplete-one-sided-pullback.json', context)
    reject('source_local_dual_incomplete_pullback',
           lambda: build_readout(entry['pair']['owner'], partial, seed, proof))
    owner = bound_json(entry['pair']['owner'])
    owner['request']['selectors']['axis'] = 'z'
    foreign = save('foreign-one-sided-owner.json', owner)
    reject('source_local_dual_foreign_owner',
           lambda: build_readout(foreign, pair['covariance']['left_context'], seed, proof))


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


def exercise_dual_join(entries, seed, discrete_checks, consume, reference, save, close, reject):
    sources = [file_binding(HERE/name) for name in
               ('joint_spatial_dual_coherent.py', 'joint_spatial_dual_authority.py')]
    proof = save('computed-dual-oracle.json', dict(schema='P254-joint-spatial-readout-oracle-v1',
        status='passed', source_bindings=sources, checks=discrete_checks,
        scope='executed finite discrete-forward/dual numerical comparisons'))
    dual_entries = []
    for index, entry in enumerate(entries):
        native, right = _without_noise(entry, seed, proof)
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
    _reject_incomplete(entries[0], seed, proof, save, reject)
    _reject_changed(dual_entries, seed, consume, save, reject)
