"""New finite-sensitive oracle only; original217.705s A44 proofs are carried.

Independent load-bearing route is the UNREDUCED existing physical per-external
four-term polynomial, not the new DD-DX/bubble/matrix-chain assembly.
All eight explicit spin/k externals, four q (nonzero and Nyquist), two bands,
unequal feature closures, all16 ordered pairs, independently complex legs.
"""
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from time import perf_counter
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'0044'))
from bare_coulomb_factors import BareCoulombFactors
from material_second_born import evaluate_material
from second_born_functional import terms_from_vertices,ordered
from all_slot_coulomb_contractions import (KEYS,ordered_q,evaluate_time_grid,
                                           validate_physical,source_bounds)
from all_slot_coulomb_checkpoint import (binding,atomic_json,AllSlotPrefix,
                                         publish_result,consume_result)


def run():
    started = perf_counter()
    rng = np.random.default_rng(4606)
    random = lambda s:(rng.normal(size=s)+1j*rng.normal(size=s))*.2
    gs = (np.array([[0,0,0],[1,0,0],[-1,0,0]]),
          np.array([[0,0,0],[1,0,0]]),
          np.array([[0,0,0],[0,0,1],[1,0,0],[-1,0,1]]),
          np.array([[0,0,0],[-1,0,0]]))
    partners = (np.array([0,2,1]),np.array([0,1]),np.array([1,0,3,2]),np.array([0,1]))
    mesh = np.array([[0,0,k] for k in range(4)])
    target = np.array([[s//4*4+(s%4+q)%4 for s in range(8)] for q in range(4)])
    density = tuple(random((8,len(g),2,2)) for g in gs)
    c = SimpleNamespace(n=4,ns=8,nq=4,m=2,density=density,slot_spin=np.repeat([1,2],4),
        slot_k=np.tile(np.arange(4),2),k_indices=mesh,q_indices=mesh,target=target,
        minus=np.array([0,3,2,1]),reciprocal=np.eye(3),tpiba=1.,feature_g=gs,
        feature_partner=partners)
    factors = BareCoulombFactors(c)
    validate_physical(factors)
    times = np.array([0.,.1,.4,.9])
    electrons = {(a,b):(random((8,2,2)),random((8,2,2))) for a in range(4) for b in range(4)}
    out = HERE/'finite-oracle-output'
    out.mkdir(exist_ok=False)
    fixture = out/'actual-independent-inputs.npz'
    np.savez(fixture,**{f'density_{q}':d for q,d in enumerate(density)},
        **{f'leg_{a}_{b}_{o}':e[o] for (a,b),e in electrons.items() for o in range(2)})
    sources = [binding(p) for p in sorted(HERE.glob('*.py'))]
    sources += [binding(HERE.parent/'0044'/name) for name in
        ('bare_coulomb_factors.py','second_born_functional.py','material_second_born.py')]
    spec = dict(schema='p254-a46-sensitive-finite-oracle-input-v1',nq=4,shape=[8,4,4,2,2],
        external_batches=[[0,1,2],[3,4,5],[6,7]],
        authenticated_sources_and_inputs=[*sources,binding(fixture)],
        times=times.tolist(),feature_counts=[len(g) for g in gs],
        original_definition='A44 four expanded terms; arbitrary independent ordered/reversed Green legs')
    checkpoint = AllSlotPrefix(out/'prefix',spec,spec['shape'],4)
    def interrupt(q):
        if q == 1:
            raise InterruptedError('intentional after first durable all-external q')
    try:
        evaluate_time_grid(factors,electrons,times,checkpoint,external_batch=3,after_commit=interrupt)
        raise AssertionError('partial interruption not exercised')
    except InterruptedError:
        pass
    partial,done = AllSlotPrefix(out/'prefix',spec,spec['shape'],4).restore()
    assert done == 1 and consume_result(out,spec) is None
    resumed = AllSlotPrefix(out/'prefix',spec,spec['shape'],4)
    actual,counts = evaluate_time_grid(factors,electrons,times,resumed,external_batch=3)
    assert counts['new_q'] == 3 and counts['reused_q'] == 1
    expected = {k:np.zeros((8,4,4,2,2),complex) for k in KEYS}
    for i in range(8):
        for a in range(4):
            for b in range(4):
                value,count = evaluate_material(factors,i,electrons[a,b],electrons[b,a],retarded=times[a]>=times[b])
                assert count['momentum_spin_quartets'] == 48
                for order,channels in value.items():
                    for channel,x in channels.items():
                        expected[order+'_'+channel][i,a,b] = x
    errors = {k:dict(absolute=float(np.linalg.norm(actual[k]-expected[k])),
        relative=float(np.linalg.norm(actual[k]-expected[k])/max(1e-30,np.linalg.norm(expected[k])))) for k in KEYS}
    assert max(x['relative'] for x in errors.values()) < 2e-11,errors
    # Exact saved result decoding exercises the public consumer's numerical branch.
    result = publish_result(out,spec,actual,counts,resumed.recoveries)
    assert all(np.array_equal(result['arrays'][k],actual[k]) for k in KEYS)
    completed,reuse = evaluate_time_grid(factors,electrons,times,
        AllSlotPrefix(out/'prefix',spec,spec['shape'],4),external_batch=3)
    assert reuse['new_q'] == 0 and reuse['reused_q'] == 4
    assert all(np.array_equal(completed[k],actual[k]) for k in KEYS)
    # Simulate a killed overwrite/corrupted newest generation: previous q3 survives.
    newest = out/'prefix'/'prefix-0.npz'
    with newest.open('ab') as f:
        f.write(b'interrupted-write-mutation')
    recovery = AllSlotPrefix(out/'prefix',spec,spec['shape'],4)
    fallback,fallback_q = recovery.restore()
    assert fallback_q == 3 and recovery.recoveries
    recovered,repair = evaluate_time_grid(factors,electrons,times,recovery,external_batch=3)
    assert repair['new_q'] == 1 and repair['reused_q'] == 3
    assert all(np.array_equal(recovered[k],actual[k]) for k in KEYS)
    # Mutations at the actual sensitive layer, no physical-state symmetry masks them.
    fwd,bwd = electrons[0,1],electrons[1,0]
    base = expected['greater_total'][:,0,1]
    wrong_reverse = sum((sum(ordered_q(factors,q,fwd[0],fwd[1],external_batch=3))
                         for q in range(4)),np.zeros_like(base))
    headless = expected['greater_direct'][:,0,1]
    no_nonzero_q = sum(ordered_q(factors,0,fwd[0],bwd[1],external_batch=3))
    original_weights = factors.weights[0].copy()
    factors.weights[0][0] = 1.
    wrong_head = sum((sum(ordered_q(factors,q,fwd[0],bwd[1],external_batch=3))
                      for q in range(4)),np.zeros_like(base))
    factors.weights[0] = original_weights
    original_partner = c.feature_partner
    c.feature_partner = (*partners[:2],np.arange(4),partners[3])
    wrong_partner = sum((sum(ordered_q(factors,q,fwd[0],bwd[1],external_batch=3))
                         for q in range(4)),np.zeros_like(base))
    c.feature_partner = original_partner
    one_spin = tuple(x.copy() for x in bwd)
    one_spin[1][4:] = 0
    wrong_spin = sum((sum(ordered_q(factors,q,fwd[0],one_spin[1],external_batch=3))
                      for q in range(4)),np.zeros_like(base))
    mutations = {name:float(np.linalg.norm(value-base)/np.linalg.norm(base)) for name,value in
        (('reverse_leg_inference',wrong_reverse),('exchange_omission',headless),
         ('nonzero_q_omission',no_nonzero_q),('bare_head_insertion',wrong_head),
         ('Nyquist_partner_omission',wrong_partner),('spin_sector_omission',wrong_spin))}
    assert min(mutations.values()) > 1e-5,mutations
    # Generic arbitrary vertices have NO physical pair license: expose false reduction.
    left = tuple((random((3,2,2)),random((3,2,2))) for _ in range(2))
    right = tuple((random((3,2,2)),random((3,2,2))) for _ in range(2))
    generic = terms_from_vertices(left,right)
    g = [random((2,2)) for _ in range(3)]
    unreduced = ordered(generic,*g)['total']
    false_reduced = ordered((generic[0],generic[1]),*g)['total']*2
    generic_counterexample = float(np.linalg.norm(unreduced-false_reduced))
    assert generic_counterexample > 1e-8
    report = dict(schema='p254-a46-new-sensitive-all-slot-oracle-v1',status='passed',
        seconds=perf_counter()-started,errors=errors,mutation_relative_defects=mutations,
        generic_vertex_false_reduction_defect=generic_counterexample,
        partial_resume=counts,completed_resume=reuse,previous_generation_recovery=repair,
        recovery_receipts=recovery.recoveries,input_and_source_specification=spec,
        saved_result=result['record']['archive'],decoded_norms=result['decoded_norms'],
        full_material_bounds=source_bounds(),
        scope='all8 finite externals/all4 q/all16 pairs/2bands/unequal3,2,4,2features/two spins; independent asymmetric complex legs; original physical polynomial; actual numeric recovery',
        carried_original_proofs_seconds=217.705,material_speed='UNMEASURED',
        findings=[dict(claim='physical all-slot sharing matches original unreduced four-term functional',status='passed'),
                  dict(claim='partial/completed/previous-generation numeric consumers recover actual values',status='passed')],
        blockers=[])
    atomic_json(HERE/'finite-oracle-report.json',report)
    print(json.dumps(report),flush=True)


if __name__ == '__main__':
    run()
