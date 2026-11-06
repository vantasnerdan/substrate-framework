"""Bookkeeping only: include final consumer's completed native receipt."""
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def run():
    from refinement_common import binding,save
    runs=[json.loads(p.read_bytes()) for p in sorted((HERE/'invocations').glob('*/exit.json'))]
    total=sum(r['elapsed_seconds'] for r in runs)
    ee=sum(r['elapsed_seconds'] for r in runs if r['mode']=='ee')
    cubic=sum(r['elapsed_seconds'] for r in runs if r['mode'].startswith('cubic') or r['mode']=='merge')
    grant_path=HERE/'direct-cubic-continuation-grant.json'
    cubic_cap=json.loads(grant_path.read_bytes())['amended_combined_cubic_cap_seconds'] if grant_path.exists() else 9000
    ledgers=dict(all_new_native_seconds=total,all_new_native_cap_seconds=6000,all_new_ee_native_seconds=ee,A46_original_material_seconds=6633.830,A46_material_cumulative_seconds=6633.830+ee,A46_material_cap_seconds=24000,A46_WaveI_separate_seconds=6.722,cubic_new_native_seconds=cubic,cubic_combined_original_seconds=8611.287,cubic_combined_cumulative_seconds=8611.287+cubic,cubic_combined_cap_seconds=cubic_cap,cubic_reader_accounting='all component merge/consumer invocation elapsed conservatively included')
    payload=sum(p.stat().st_size for p in HERE.rglob('*') if p.is_file())
    failures=[r for r in runs if r['exitcode']]
    ledger=dict(schema='p254-a55-final-cumulative-native-ledger-v1',all_new_invocations=runs,ledgers=ledgers,failures_included=failures,new_owned_payload_bytes=payload,all_native_argv_streams_exits=[dict(launch=binding(HERE/'invocations'/r['name']/'launch.json'),exit=binding(HERE/'invocations'/r['name']/'exit.json'),stdout=r['stdout'],stderr=r['stderr']) for r in runs])
    if total>6000 or 8611.287+cubic>cubic_cap or payload>256*1024**2:raise RuntimeError('actual cap exceeded')
    save(HERE/'actual-final-cumulative-ledger-001.json',ledger)
    handoff=HERE/'actual-genuine-midpoint-handoff-001.json'
    if handoff.exists():
        record=json.loads(handoff.read_bytes());record.update(all_new_invocations=runs,ledgers=ledgers,new_owned_payload_bytes=payload,final_ledger=binding(HERE/'actual-final-cumulative-ledger-001.json'),native_failures=failures)
        if grant_path.exists():
            record['direct_cubic_continuation_grant']=binding(grant_path)
            record['findings'].append('Root577 amended cubic ceiling9000->9500 append-only; genuine paid q0 reused and direct q1..63 continued without symmetry reduction or field/source/tolerance changes')
        save(handoff,record)
    return ledger
if __name__=='__main__':run()
