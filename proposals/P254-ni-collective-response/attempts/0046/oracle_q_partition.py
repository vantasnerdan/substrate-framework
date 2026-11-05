"""Only new exact q-list ownership, split/merge saved values and recovery boundary."""
import json
from pathlib import Path
from types import SimpleNamespace
from time import perf_counter
import numpy as np
from all_slot_coulomb_checkpoint import binding,atomic_json,AllSlotPrefix,publish_result
from all_slot_coulomb_partition import partition_specification,merge_partitions,evaluate_partition
from all_slot_coulomb_shared import evaluate_time_grid
from run_all_slot_coulomb import BareCoulombFactors
HERE=Path(__file__).resolve().parent


def run():
    start=perf_counter()
    out=HERE/'partition-oracle-001'; out.mkdir(exist_ok=False)
    fixture=HERE/'finite-oracle-output'/'actual-independent-inputs.npz'
    full=HERE/'shared-oracle-001'/'all-slots-all-pairs.npz'
    with np.load(fixture,allow_pickle=False) as a:
        density=tuple(a[f'density_{q}'] for q in range(4))
        electrons={(i,j):tuple(a[f'leg_{i}_{j}_{o}'] for o in range(2)) for i in range(4) for j in range(4)}
    gs=(np.array([[0,0,0],[1,0,0],[-1,0,0]]),np.array([[0,0,0],[1,0,0]]),
        np.array([[0,0,0],[0,0,1],[1,0,0],[-1,0,1]]),np.array([[0,0,0],[-1,0,0]]))
    partners=(np.array([0,2,1]),np.array([0,1]),np.array([1,0,3,2]),np.array([0,1]))
    mesh=np.array([[0,0,k] for k in range(4)])
    c=SimpleNamespace(n=4,ns=8,nq=4,m=2,density=density,slot_spin=np.repeat([1,2],4),
        slot_k=np.tile(np.arange(4),2),k_indices=mesh,q_indices=mesh,
        target=np.array([[s//4*4+(s%4+q)%4 for s in range(8)] for q in range(4)]),
        minus=np.array([0,3,2,1]),reciprocal=np.eye(3),tpiba=1.,feature_g=gs,feature_partner=partners)
    common=dict(nq=4,shape=[8,4,4,2,2],external_batches=[[0,1],[2,3],[4,5],[6,7]],
        authenticated_sources_and_inputs=[binding(fixture),binding(full),
            *[binding(HERE/n) for n in ('all_slot_coulomb_shared.py','all_slot_coulomb_partition.py',
                'all_slot_coulomb_checkpoint.py','oracle_q_partition.py')]])
    factors=BareCoulombFactors(c); times=[0.,.1,.4,.9]
    parts=[]; recovery=[]
    for name,owned in (('even',[0,2]),('odd',[1,3])):
        directory=out/name; directory.mkdir()
        spec=partition_specification(common,owned)
        cache=AllSlotPrefix(directory/'numeric-prefix',spec,spec['shape'],len(owned))
        def interrupt(count):
            if count==1: raise InterruptedError('new q-list prefix committed first owned globalq')
        try:
            evaluate_time_grid(factors,electrons,times,cache,external_batch=2,time_batch=4,after_commit=interrupt)
            raise AssertionError('new q-list interrupted prefix not exercised')
        except InterruptedError: pass
        arrays,done=AllSlotPrefix(directory/'numeric-prefix',spec,spec['shape'],len(owned)).restore()
        assert done==1
        side=json.loads((directory/'numeric-prefix'/'prefix-1.json').read_bytes())
        assert side['owned_q']==owned[:1] and side['exchange_inner_h_count']==4
        result,counts=evaluate_partition(factors,electrons,times,directory,spec)
        assert counts['new_q']==1 and counts['reused_q']==1
        again,reuse=evaluate_partition(factors,electrons,times,directory,spec)
        assert reuse['new_q']==0
        assert result['record']['contribution_scope']=='PARTIAL_Q_CONTRIBUTION'
        recovery.append(dict(partition=owned,resume=counts,completed_reuse=reuse,
            actual_prefix_owned_q=side['owned_q'],saved=result['record']['archive']))
        parts.append(directory)
    merged=merge_partitions(parts,out/'merged')
    with np.load(full,allow_pickle=False) as a:
        errors={k:float(np.linalg.norm(merged['arrays'][k]-a[k])/max(1e-30,np.linalg.norm(a[k]))) for k in a.files}
    assert max(errors.values())<2e-11,errors
    rejected={}
    for name,directories in (('missing',[parts[0]]),('overlap',[parts[0],parts[0]])):
        try:
            merge_partitions(directories,out/('invalid-'+name))
            raise AssertionError('invalid partition cover accepted')
        except ValueError as e: rejected[name]=str(e)
    assert merged['record']['contribution_scope']=='FULL_ALL_Q'
    report=dict(schema='p254-a46-minimal-q-partition-saved-merge-oracle-v1',status='passed',
        seconds=perf_counter()-start,relative_errors=errors,recovery=recovery,
        rejected_missing_overlap=rejected,merged_result=merged['record']['archive'],
        merged_counts=merged['record']['counts'],decoded_norms=merged['decoded_norms'],
        source_inputs=common['authenticated_sources_and_inputs'],
        scope='Only newly changed exact q-list partition/prefix/partial-label/unscaled complete merger boundary; previous full-width/physical suites carried, not replayed',
        findings=[dict(claim='disjoint noncontiguous q lists sum actual saved all-slot16-pair numeric values once, no scaling',status='passed'),
                  dict(claim='partial and completed q-list numeric recovery retains actual globalq ownership/all innerh',status='passed')],blockers=[])
    atomic_json(out/'report.json',report)
    print(json.dumps(report),flush=True)

if __name__=='__main__': run()
