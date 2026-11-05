"""Decode the actual saved full-sum cubic output; never construct jets/kernels."""
import argparse
import json
from pathlib import Path
import numpy as np
from substrate_framework.streamed_material_state import authenticate,file_binding


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--result',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    source_binding=file_binding(args.source); result_binding=file_binding(args.result)
    source=json.loads(args.source.read_bytes()); result=json.loads(args.result.read_bytes())
    if (source['status']!='complete' or result['status']!='complete'
            or result['third_source']!=source_binding
            or len(source['full_sector_lineage'])!=128
            or result['external_q']!=[0,1,2]
            or result['external_rows']!=[0,1,2,93,94,95]
            or result['times_inverse_ry']!=[0.,1.,10.,30.]
            or result['internal_sums']['internal_q_count']!=64
            or result['internal_sums']['internal_first_mode_count']!=111
            or result['internal_sums']['internal_second_mode_count']!=111
            or result['internal_sums']['routed_pair_evaluations']!=3072):
        raise ValueError('actual completed full-scope original-source consumer required')
    with np.load(authenticate(source['arrays']),allow_pickle=False) as saved:
        source_norms={name:float(np.linalg.norm(saved[name])) for name in ('local','nonlocal_term')}
    quantities={}
    with np.load(authenticate(result['output']),allow_pickle=False) as saved:
        if not np.array_equal(saved['external_q'],[0,1,2]) or not np.array_equal(saved['external_rows'],[0,1,2,93,94,95]) or not np.array_equal(saved['times_inverse_ry'],[0.,1.,10.,30.]):
            raise ValueError('saved numeric axes differ')
        for name in ('greater','lesser','retarded','mean_force'):
            value=saved[name]
            shape=(4,3,6) if name=='mean_force' else (4,4,3,6,6)
            if value.shape!=shape or not np.isfinite(value).all():
                raise ValueError('nonfinite/incomplete actual saved numeric output')
            first=value.reshape(-1)[0]
            quantities[name]=dict(shape=list(value.shape),norm=float(np.linalg.norm(value)),
                max_abs=float(np.max(np.abs(value))),first_real=float(first.real),first_imag=float(first.imag))
    record=dict(schema='p254-a43-actual-saved-numerical-consumer-v1',status='complete',
        source=source_binding,result=result_binding,source_arrays=source['arrays'],kernel_arrays=result['output'],
        source_norms=source_norms,quantities=quantities,scope=dict(sectors=128,external_q=[0,1,2],external_rows=[0,1,2,93,94,95],internal_q=64,internal_modes=[111,111],times=[0.,1.,10.,30.],time_pairs=16),
        electronic_atom_indices=source['electronic_atom_indices'],
        electronic_coverage='Only declared external atoms; not an all37electronic-third source',
        consumer=file_binding(__file__),
        units=result['units'],full_physical_goal_complete=False,lifetime_or_rate_computed=False)
    with args.output.open('x') as stream:
        json.dump(record,stream,indent=2,sort_keys=True,allow_nan=False); stream.write('\n')
    print(json.dumps(record),flush=True)

if __name__=='__main__':
    main()
