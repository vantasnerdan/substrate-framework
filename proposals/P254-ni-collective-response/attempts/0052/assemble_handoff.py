"""Bind original A52 artifacts and close only this completed owned jobview."""
import json,os,sys,time
from pathlib import Path
from control import HERE,OUT,capture,binding,save

def main():
    started=time.monotonic();base=Path('/var/tmp/p254-a52-requested-reference-001')
    receipt=Path(str(base)+'-receipts')
    terminal=json.loads((receipt/'exit.json').read_bytes())
    if terminal['exitcode']!=0:raise RuntimeError('complete actual science required')
    foreground_code=capture('completed-foreground',['/home/administrator/.local/bin/herdr','pane','process-info','--pane','w1:p6N'])
    foreground=json.loads((OUT/'completed-foreground.stdout').read_bytes())['result']['process_info']
    if foreground_code or not all(p['argv']==['/bin/bash'] for p in foreground['foreground_processes']):raise RuntimeError('completed foreground shell required before closing owned view')
    if capture('close-completed-owned-view',['/home/administrator/.local/bin/herdr','pane','close','w1:p6N']):raise RuntimeError('owned view close failed')
    sources=[binding(p) for p in sorted(HERE.iterdir()) if p.is_file()]
    artifacts=[binding(p) for directory in (base,receipt,OUT) for p in sorted(directory.rglob('*')) if p.is_file()]
    producer=json.loads((base/'producer/result.json').read_bytes());consumer=json.loads((base/'consumer.json').read_bytes());oracle=json.loads((base/'oracle.json').read_bytes())
    runs=[json.loads(p.read_bytes()) for p in sorted(base.parent.glob('p254-a52-requested-reference-*-receipts/exit.json'))]
    handoff=dict(schema='p254-a52-original-uniform-moving-reference-handoff-v1',status='complete',blockers=[],
        outputs=dict(bundle=binding(base/'producer/bundle.json'),U=binding(base/'producer/U.npy'),hamiltonian_coefficients=binding(base/'producer/hamiltonian-chebyshev.npz'),producer=binding(base/'producer/result.json'),saved_consumer=binding(base/'consumer.json'),oracle=binding(base/'oracle.json'),input_inventory=binding(base/'input-inventory.json')),
        source_bindings=sources,artifact_inventory=artifacts,
        findings=['Genuine same saved moving-mean Hamiltonian timeordered U for all128spin-k/12bands at all31integer nodes0..30; no density-inferred U or kernel interpolation.',
            'Original full mechanical S is direct expm on demand with128MiB byte LRU plus one selected q/-q history pair. FullGamma228phase and114positions/3field selections exposed; unstable modes and literal preparation retained.',
            'Existing A41 electron/boson/transitions callback shape and existing A50 ordered_coordinates consumed actual saved reference; no U integration in saved consumer.',
            'Exact requested fractional producer times are genuinely DOP853 integrated; saved reader accepts exact retained nodes only and rejects missing nodes.',
            'New finite rotating noncommuting-H plus arbitrary paired-q/fullfield Gaussian oracle passed, then selected actual q0/1/3 all31S and positions consumed. Old proof campaigns not replayed.'],
        old_four_node_raw_U_comparison=producer['old_four_node_U_comparison'],old_four_node_raw_S_comparison=consumer['old_four_node_S_comparison'],
        old_comparison_policy=producer['comparison_policy'],oracle_summary=oracle,
        actual_S_cache=dict(expm_calls=consumer['S_expm_calls'],cache_hits=consumer['S_cache_hits'],cached_bytes=consumer['S_cached_bytes'],cache_cap=consumer['S_cache_cap']),
        actual_existing_A50_consumption=consumer['A50_existing_ordered_coordinates'],saved_U_recompute_calls=consumer['saved_U_integration_calls'],
        costs=dict(all_completed_native_receipts=runs,cumulative_completed_native_seconds=sum(r['elapsed_seconds'] for r in runs),science_stage_receipts=json.loads((base/'wave-exit.json').read_bytes()),final_metadata_current_seconds_before_handoff=time.monotonic()-started,final_metadata_receipt='bound after terminal by attempts/0052/actual-final-native-binding-002.json; not free recovery'),
        failed_native_attempts=[r for r in runs if r['exitcode']],
        allocation=dict(cpu=7,memory_max=3221225472,swap=0,threads=1,restart='no',per_attempt_cap_seconds=300,cumulative_cap_seconds=600,payload_cap_bytes=67108864),
        payload_bytes_before_handoff=sum(item['bytes'] for item in artifacts)+sum(item['bytes'] for item in sources),
        owned_jobview=dict(pane='w1:p6N',created_before_native=True,foreground_receipt=binding(OUT/'completed-foreground.exit.json'),closure_receipt=binding(OUT/'close-completed-owned-view.exit.json')),
        limits=['Bare Gaussian reference only, not electron/cubic/collision feedback.','No uniform selfenergy history, selfconsistent closure, rigorous global accuracy, physicalrates/fullE or PBE/XC counterterm claim.','New31node U integration has paid cost and changed interval boundaries at same old precision; oldnodes not forced or unitarily repaired.','Requested quantitative-verification/code-structure skills unavailable in mounted inventory; calls failed without invented contents.'])
    path=HERE/'actual-original-uniform-moving-reference-handoff-001.json';save(path,handoff)
    total=handoff['payload_bytes_before_handoff']+path.stat().st_size
    if total>=67108864:raise RuntimeError('complete payload cap exceeded')
    print(json.dumps(dict(handoff=binding(path),total_payload_bytes=total)),flush=True)
if __name__=='__main__':main()
