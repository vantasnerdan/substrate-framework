"""After EE exits, run independent ready components with total-native reserve."""
from pathlib import Path
import json,subprocess,sys
HERE=Path(__file__).resolve().parent
PYTHON='/home/administrator/projects/substrate-framework/.venv/bin/python'

def main():
    from refinement_common import save
    exits=[json.loads(p.read_bytes()) for p in (HERE/'invocations').glob('*/exit.json')]
    ee=[x for x in exits if x['mode']=='ee']
    if len(ee)!=3 or any(x['exitcode'] for x in ee):raise RuntimeError('EE native exits must persist before remaining work')
    used=sum(x['elapsed_seconds'] for x in exits)
    allocation={'ei':1200,'joint':180,'represented':600,'cubic':380,'merge':60}
    if used+sum(allocation.values())>5990:raise RuntimeError('remaining component upper allocations do not fit totalNative6000')
    save(HERE/'remaining-wave-decision.json',dict(actual_paid_new_native_seconds=used,remaining_seconds=6000-used,allocated_caps_seconds=allocation,action='launch independent EI/core6 joint/core7 cubic/core8; represented follows joint oncore7',no_core_overlap=True))
    def command(mode,name,cpu,cap):return [PYTHON,str(HERE/'launch_component.py'),mode,name,str(cpu),str(cap)]
    ei=subprocess.Popen(command('ei','ei-001',6,allocation['ei']))
    joint=subprocess.Popen(command('joint','joint-001',7,allocation['joint']))
    cubic=subprocess.Popen(command('cubic','cubic-001',8,allocation['cubic']))
    codes={'joint':joint.wait()}
    represented=subprocess.Popen(command('represented','represented-001',7,allocation['represented']))
    codes.update(ei=ei.wait(),cubic=cubic.wait(),represented=represented.wait())
    print(json.dumps(codes),flush=True)
    codes['merge']=subprocess.run(command('merge','merge-001',6,allocation['merge'])).returncode
    save(HERE/'remaining-wave-exit.json',dict(exitcodes=codes))
    from finalize_receipts import run
    run()
    print('A55_REMAINING_WAVE_001_TERMINAL',flush=True);return int(any(codes.values()))
if __name__=='__main__':sys.exit(main())
