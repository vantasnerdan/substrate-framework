"""Aggregate actual completed findings/costs, without replaying physics work."""
import sys,json,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'0055'))
from refinement_common import binding,save

def main():
    out=HERE/'material/common-mean-001'
    result=json.loads((out/'result.json').read_bytes())
    original_result=binding(out/'result.json')
    result['reference_source_detail']=result['reference_source_kind']
    result['reference_source_kind']='bare-mechanical-with-contacts'
    result['producer_result']=original_result
    result['publication_source']=binding(__file__)
    save(out/'reference-result.json',result)
    consumer=json.loads((out/'saved-consumer.json').read_bytes())
    oracle=json.loads((out/'boundary-oracle.json').read_bytes())
    inv=HERE/'invocations'
    actual=[]
    for name in ('common-mean-001','common-mean-002'):
        folder=inv/name
        entries={p.name:binding(p) for p in folder.iterdir() if p.is_file()}
        exitrec=json.loads((folder/'launcher-exit.json').read_bytes())
        actual.append(dict(name=name,bindings=entries,launcher_exit=exitrec,native_exit=json.loads((folder/'native-exit.json').read_bytes()) if (folder/'native-exit.json').exists() else None))
    findings=[dict(category=category,name=k,value=v,evidence_tier='scoped numerical diagnostic') for category,data in (('raw_defect',result['defects']),('new_composition_boundary',oracle),('saved_consumer',consumer)) for k,v in data.items()]
    query=['systemctl','--user','show','p254-a57-common-mean-002.service','--property=ActiveState,SubState,Result,ExecMainCode,ExecMainStatus,CPUUsageNSec,MemoryPeak,MemorySwapPeak,InvocationID']
    terminal=subprocess.run(query,capture_output=True,text=True)
    payload=sum(p.stat().st_size for p in HERE.rglob('*') if p.is_file())
    record=dict(schema='p254-a57-actual-common-mean-feedback-handoff-v1',status='complete',result=binding(out/'reference-result.json'),original_producer_result=original_result,publication='Metadata-only canonical source-kind label; identical archive/source/input/material bindings. Original saved-consumer receipt remains bound to original producer result.',archive=result['archive'],saved_consumer=binding(out/'saved-consumer.json'),boundary_oracle=binding(out/'boundary-oracle.json'),execution_contract=binding(HERE/'execution-contract-002.json'),prior_execution_contract=binding(HERE/'execution-contract.json'),
        input_dataset=result['input_dataset'],coefficient_manifest=result['coefficient_manifest'],baseline_catalog=result['baseline_catalog'],source_bindings=result['source_bindings'],physical_history_bindings=result['physical_history_bindings'],
        API=dict(module=str(HERE/'mean_feedback_reference.py'),consumer='MeanFeedbackReference(c,result_path,coefficient_binding)',producer='produce_mean_feedback(dataset,baseline_catalog,output,*,rtol=2e-8,atol=2e-10)',producer_module=str(HERE/'produce_mean_feedback.py'),protocol=['c','times','u','bindings','electron(a,b)','ordered(q,a,b,selection="positions")','boson(q,a,b)','S(q,t)','transitions(q)','close()']),
        findings=findings,reference_source_kind=result['reference_source_kind'],approximation=result['approximation'],bath_scope=result['bath_scope'],limits=result['limits'],
        costs=dict(invocations=actual,preunit_failed_invocation=binding(inv/'preunit-failed-001.json'),cumulative_native_launcher_seconds=sum(x['launcher_exit']['elapsed_seconds'] for x in actual),payload_bytes_before_handoff=payload,payload_cap_bytes=67108864,native_cap_seconds=1200),
        final_unit_query=dict(argv=query,exitcode=terminal.returncode,stdout=terminal.stdout,stderr=terminal.stderr),new_oracle_replay=False,admitted_proofs_carried=True)
    save(HERE/'actual-common-mean-feedback-handoff-001.json',record)
    print(json.dumps(dict(handoff=binding(HERE/'actual-common-mean-feedback-handoff-001.json'),costs=record['costs'],terminal=record['final_unit_query'])),flush=True)
if __name__=='__main__':main()
