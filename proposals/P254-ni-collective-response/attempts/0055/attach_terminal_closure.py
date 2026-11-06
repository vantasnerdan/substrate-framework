"""Attach actual terminal closure receipts; no numerical replay."""
from refinement_common import *

def run():
    names=('foreground-final-001','read-final-001','close-owned-pane-001','current-final-001')
    evidence={name:[binding(HERE/'control'/(name+suffix)) for suffix in ('.start.json','.stdout','.stderr','.exit.json')] for name in names}
    for name in names:
        result=json.loads((HERE/'control'/(name+'.exit.json')).read_bytes())
        if result['exitcode']:raise RuntimeError('terminal closure control failure: '+name)
    text=(HERE/'control/read-final-001.stdout').read_text()
    if 'A55_DIRECT_FULL_CATALOG_002_TERMINAL_EXIT_0' not in text:raise RuntimeError('fresh actual final consumer marker absent from saved terminal view')
    payload=sum(p.stat().st_size for p in HERE.rglob('*') if p.is_file())
    if payload+1048576>256*1024**2:raise RuntimeError('final payload reserve exceeds cap')
    record=dict(status='complete',owned_pane='w1:p6R',final_actual_marker='A55_DIRECT_FULL_CATALOG_002_TERMINAL_EXIT_0',receipts=evidence,unrelated_original_panes_untouched=True,new_owned_payload_bytes_before_closure_attachment=payload,final_payload_bound_bytes=payload+1048576,metadata_write_reserve_bytes=1048576)
    save(HERE/'actual-terminal-view-closure-001.json',record)
    handoff=HERE/'actual-genuine-midpoint-handoff-001.json';result=json.loads(handoff.read_bytes())
    result.update(terminal_view_closure=binding(HERE/'actual-terminal-view-closure-001.json'),terminal_closure_source=binding(__file__),final_payload_bound_bytes=record['final_payload_bound_bytes'])
    save(handoff,result)
    print(json.dumps(dict(status=result['status'],handoff=binding(handoff),catalog=result.get('catalog'),ledgers=result['ledgers'],final_payload_bound_bytes=record['final_payload_bound_bytes'])),flush=True)
if __name__=='__main__':run()
