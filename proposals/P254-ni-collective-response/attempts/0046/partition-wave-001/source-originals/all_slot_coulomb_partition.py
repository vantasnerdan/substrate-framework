"""Specific exact outer-q partition contract and unscaled complete saved merger."""
import argparse
import json
from pathlib import Path
import numpy as np
from all_slot_coulomb_checkpoint import (AllSlotPrefix,binding,authenticate,atomic_json,
                                         consume_result,publish_result)
from all_slot_coulomb_shared import evaluate_time_grid,bounds
from all_slot_coulomb_contractions import KEYS

PARTITION_FIELDS = ('nq','outer_q_partition','global_nq','contribution_scope')


def partition_specification(common,owned_q):
    result = dict(common,nq=len(owned_q),global_nq=common['nq'],
        outer_q_partition=list(owned_q),contribution_scope='PARTIAL_Q_CONTRIBUTION')
    if (not owned_q or list(owned_q)!=sorted(set(owned_q)) or min(owned_q)<0
            or max(owned_q)>=common['nq']):
        raise ValueError('partition must own distinct ascending original outerq')
    return result


def common_specification(part):
    result = {k:v for k,v in part.items() if k not in PARTITION_FIELDS}
    result['nq'] = part['global_nq']
    return result


def evaluate_partition(factors,electrons,times,directory,specification,*,after_commit=None):
    if specification['contribution_scope']!='PARTIAL_Q_CONTRIBUTION':
        raise ValueError('worker must label partial outer-q contribution')
    old = consume_result(directory,specification)
    if old is not None:
        return old,dict(new_q=0,reused_q=specification['nq'])
    cache = AllSlotPrefix(Path(directory)/'numeric-prefix',specification,
        specification['shape'],specification['nq'])
    arrays,counts = evaluate_time_grid(factors,electrons,times,cache,
        external_batch=2,time_batch=4,after_commit=after_commit)
    result = publish_result(directory,specification,arrays,counts,cache.recoveries)
    return result,counts


def merge_partitions(directories,output):
    """Authenticate/decode actual saved contributions; disjoint completeq only."""
    records = []
    for directory in directories:
        path = Path(directory)/'result.json'
        record = json.loads(path.read_bytes())
        spec = record['specification']
        if (record['contribution_scope']!='PARTIAL_Q_CONTRIBUTION'
                or spec['contribution_scope']!='PARTIAL_Q_CONTRIBUTION'):
            raise ValueError('merger accepts only explicit partial-q contributions')
        records.append((spec['outer_q_partition'][0],Path(directory),record,binding(path)))
    records.sort(key=lambda x:x[0])
    if not records:
        raise ValueError('no partition results')
    common = common_specification(records[0][2]['specification'])
    allq = []
    for _,_,record,_ in records:
        spec = record['specification']
        if common_specification(spec)!=common or len(spec['outer_q_partition'])!=spec['nq']:
            raise ValueError('partition source/input/basis/time/external/count contracts differ')
        allq.extend(spec['outer_q_partition'])
    if sorted(allq)!=list(range(common['nq'])) or len(set(allq))!=len(allq):
        raise ValueError('partition q ownership is not disjoint complete original0..nq-1')
    arrays = {k:np.zeros(common['shape'],complex) for k in KEYS}
    parts = []
    for _,directory,record,record_binding in records:
        saved = consume_result(directory,record['specification'])
        if record['counts']['owned_outer_q']!=record['specification']['outer_q_partition']:
            raise ValueError('saved numeric count/outerq ownership differs')
        for key in KEYS:
            arrays[key] += saved['arrays'][key]
        parts.append(dict(result=record_binding,archive=record['archive'],
            outer_q=record['specification']['outer_q_partition'],
            counts=record['counts'],decoded_norms=saved['decoded_norms']))
    spec = dict(common,contribution_scope='FULL_ALL_Q',partition_results=parts)
    spec['authenticated_sources_and_inputs'] = [*common['authenticated_sources_and_inputs'],
        *[x[3] for x in records],*[x[2]['archive'] for x in records]]
    output = Path(output); output.mkdir(parents=True,exist_ok=True)
    counts = dict(completed_q=common['nq'],owned_outer_q=list(range(common['nq'])),
        partitions=len(parts),partition_q_counts=[len(x['outer_q']) for x in parts],
        external_slots=common['shape'][0],ordered_time_pairs=common['shape'][1]*common['shape'][2],
        normalization='each unscaled outer-q contribution exactly once, no division by partitions')
    return publish_result(output,spec,arrays,counts,[])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--partition',type=Path,action='append',required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    result = merge_partitions(a.partition,a.output)
    print(json.dumps(dict(status='FULL_ALL_Q_COMPLETE',archive=result['record']['archive'],
        counts=result['record']['counts'],decoded_norms=result['decoded_norms'])),flush=True)


if __name__=='__main__':
    main()
