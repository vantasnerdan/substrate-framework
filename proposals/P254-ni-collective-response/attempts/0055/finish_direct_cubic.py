"""Consume completed direct cubic and existing components once in canonical ABI."""
from pathlib import Path
import json,subprocess,sys
HERE=Path(__file__).resolve().parent
PYTHON='/home/administrator/projects/substrate-framework/.venv/bin/python'
def main():
    from refinement_common import save
    exit_record=json.loads((HERE/'invocations/cubic-direct-002/exit.json').read_bytes())
    if exit_record['exitcode']:raise RuntimeError('persisted direct cubic native failure requires diagnosis')
    actual=sum(json.loads(p.read_bytes())['elapsed_seconds'] for p in (HERE/'invocations').glob('*/exit.json'))
    if actual+60>5990:raise RuntimeError('complete actual consumer allocation does not fit')
    argv=[PYTHON,str(HERE/'launch_component.py'),'merge','merge-002','6','60']
    code=subprocess.run(argv).returncode
    from finalize_receipts import run
    run()
    save(HERE/'direct-completion-exit.json',dict(exitcode=code,consumer_argv=argv,paid_before_consumer_seconds=actual))
    print('A55_DIRECT_FULL_CATALOG_002_TERMINAL_EXIT_'+str(code),flush=True);return code
if __name__=='__main__':sys.exit(main())
