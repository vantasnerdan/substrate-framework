"""Full genuine31x31 q1 EE history for owned slots; retain paid source cells."""
from pathlib import Path
import argparse
import json
import resource
import sys
from time import perf_counter
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/n) for n in ('0066','0064','0060','0052','0046','0044','0041','0045')]
from column_common import begin,authenticate,file_binding,atomic_json,BUNDLE,CHANNELS
from coulomb_material_inputs import load_coulomb_inputs
from bare_coulomb_factors import BareCoulombFactors
from general_moving_reference import SavedMovingReference
from gaussian_ee import stream_gaussian_ee
from streamed_ee import NumericSink,read_blocks


def produce(output,*,slots=(0,64),preregistration=None):
    slots=tuple(slots)
    if not slots or len(set(slots))!=len(slots) or any(type(i) is not int or not 0<=i<128 for i in slots):
        raise ValueError('unique owned original slots0..127 required')
    preregistration=Path(preregistration or HERE/'preregistration.json')
    started=perf_counter()
    out,_=begin(output,6)
    c=load_coulomb_inputs();factors=BareCoulombFactors(c)
    reference=SavedMovingReference(c,BUNDLE,c.bindings['manifest'])
    if not np.array_equal(reference.times,np.arange(31,dtype=float)):
        raise ValueError('complete genuine U31 required')
    old_record=json.loads((HERE.parent/'0064/pilot-001/result.json').read_bytes())
    column_record=json.loads(Path('/var/tmp/p254-a66-ee-first-column-001/result.json').read_bytes())
    if (c.bindings!=old_record['inputs']['material'] or
            file_binding(BUNDLE) not in old_record['inputs']['reference'] or
            reference.spec['U'] not in old_record['inputs']['reference'] or
            column_record['source_reference']!=file_binding(BUNDLE) or
            column_record['outer_q']!=1 or column_record['status']!='complete-first-column-q1'):
        raise ValueError('same-source paid endpoint/column required')
    endpoint=authenticate(next(item for item in old_record['outputs'] if item['path'].endswith('/old-all9.npz')))
    column=authenticate(column_record['output'])
    slot_index={value:index for index,value in enumerate(slots)}
    arrays={name:np.empty((len(slots),31,31,12,12),complex) for name in CHANNELS}
    seen=np.zeros((len(slots),31,31),bool)
    with np.load(column,allow_pickle=False) as saved:
        for name in CHANNELS:arrays[name][:,:,0]=saved[name+'_total'][:,slots].transpose(1,0,2,3)
    seen[:,:,0]=True
    with np.load(endpoint,allow_pickle=False) as saved:
        for name in CHANNELS:
            old=saved[name+'_total']
            for a,old_a in ((0,0),(30,1)):
                arrays[name][:,a,30]=old[list(slots),old_a,1]
    seen[:,[0,30],30]=True
    if any(not np.array_equal(seen[0],row) for row in seen):raise ValueError('different paid ownership')
    pairs=tuple((a,b) for a in range(31) for b in range(31) if not seen[0,a,b])+((15,0),)
    overlap=[];sink=NumericSink(out/'new-blocks.npyblocks');blocks=0
    atomic_json(out/'actual-inputs.json',dict(reference=reference.bindings,material=c.bindings,
        paid_column=file_binding(column),paid_endpoints=file_binding(endpoint),
        preregistration=file_binding(preregistration),slots=slots))
    input_binding=file_binding(out/'actual-inputs.json')
    committed_slots=set()
    def checkpoint_slot(i,s):
        path=out/('paid-slot-%03d.npz'%i)
        np.savez(path,times=reference.times,slot=np.array(i),outer_q=np.array(1),
            **{name:arrays[name][s] for name in CHANNELS})
        atomic_json(out/('paid-slot-%03d.json'%i),dict(schema='p254-paid-ee-slot-v1',
            status='complete-owned-history',input_binding=input_binding,
            reference_identity=file_binding(BUNDLE),slot=i,outer_q=1,
            shape=[31,31,12,12],dtype=arrays['greater'].dtype.str,
            data=file_binding(path)))
        committed_slots.add(i)
        print('EE_Q1_SLOT_%d_COMPLETE'%i,flush=True)
    def consume(block):
        nonlocal blocks
        if block.q!=1 or any(i not in slot_index for i in block.external):
            raise ValueError('unowned full-history block')
        sink(block)
        for local,i in enumerate(block.external):
            s=slot_index[i]
            for p,(a,b) in enumerate(block.pairs):
                if seen[s,a,b]:
                    if (a,b)!=(15,0):raise ValueError('repeated paid physical cell')
                    for name in CHANNELS:
                        value=block.values[name+'_total'][local,p];old=arrays[name][s,a,b]
                        delta=float(np.max(np.abs(value-old)))
                        overlap.append(dict(slot=i,pair=[a,b],channel=name,maxabs=delta,
                            passed=bool(np.allclose(value,old,rtol=3e-12,atol=3e-12))))
                    continue
                for name in CHANNELS:arrays[name][s,a,b]=block.values[name+'_total'][local,p]
                seen[s,a,b]=True
            if seen[s].all() and i not in committed_slots:
                own_overlap=[row for row in overlap if row['slot']==i]
                if len(own_overlap)==3:
                    if not all(row['passed'] for row in own_overlap):
                        raise ValueError('owned source-overlap check failed before checkpoint')
                    checkpoint_slot(i,s)
        blocks+=1
        if blocks%16==0:
            row=dict(actual_blocks=blocks,completed_owned_cells=int(seen.sum()),
                required_owned_cells=len(slots)*961,elapsed_seconds=perf_counter()-started)
            atomic_json(out/'progress.json',row);print(json.dumps(row),flush=True)
    tick=perf_counter()
    manifest=stream_gaussian_ee(factors,reference,reference.times,consume,owned_q=(1,),
        owned_external=slots,owned_pairs=pairs,pair_batch=8,inner_batch=8,internal_batch=8,
        source_bindings=reference.bindings)
    elapsed=perf_counter()-tick;stream=sink.close(manifest)
    if not seen.all() or len(overlap)!=3*len(slots) or not all(row['passed'] for row in overlap):
        raise ValueError('full owned square or actual long-time overlap comparison incomplete')
    if not all(np.isfinite(value).all() for value in arrays.values()):raise ValueError('nonfinite physical output')
    output_path=out/'full-square-q1.npz'
    np.savez(output_path,times=reference.times,slots=np.array(slots),outer_q=np.array(1),**arrays)
    actual_blocks=sum(1 for _ in read_blocks(out/'new-blocks.npyblocks'))
    if actual_blocks!=manifest['blocks']:raise ValueError('saved numeric stream incomplete')
    with np.load(output_path,allow_pickle=False) as saved:
        consumer=dict(final_retarded_rho_action_norm=float(np.linalg.norm(saved['retarded'][:,-1,0]@c.rho0[list(slots)])))
    record=dict(schema='p254-actual-owned-full31-ee-square-v1',status='complete-at-owned-slots-q1',
        output=file_binding(output_path),inputs=file_binding(out/'actual-inputs.json'),
        slots=slots,outer_q=1,inner_h=64,bands=12,internal_modes=1728,
        actual_owned_cells=int(seen.sum()),genuine_times_inverse_ry=reference.times.tolist(),
        evaluation_seconds=elapsed,elapsed_seconds=perf_counter()-started,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        manifest=manifest,stream=stream,overlap_findings=overlap,consumer=consumer,
        source_bindings=[file_binding(HERE/'produce_square.py'),file_binding(HERE.parent/'0064/gaussian_ee.py'),
                         file_binding(HERE.parent/'0064/time_gram.py')],
        checkpoint_slots=sorted(committed_slots),
        recovery='Complete paid-slot NPZ/JSON pairs retained. Relaunch only absent external slots; an uncommitted slot may be recomputed, no inner-mode/warm-history reuse promised.',
        limits=['Complete31x31 for the named external slots/q1 only, not the full outer mesh.',
                'Original bare Gaussian first map only; no dressed reconstruction, rank truncation or physical closure.',
                'Native total and separate global-direct/per-slot times are observations, not a demonstrated full-mesh speedup.'])
    atomic_json(out/'result.json',record);reference.close()
    print('A67_NATIVE_COMPLETE '+json.dumps(file_binding(out/'result.json')),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--slots',type=int,nargs='+',default=(0,64))
    parser.add_argument('--preregistration',type=Path)
    args=parser.parse_args()
    produce(args.output,slots=tuple(args.slots),preregistration=args.preregistration)
