"""Launch actual single-core units with append-only allocation receipts."""
import os,subprocess,sys,time,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
PYTHON='/home/administrator/projects/substrate-framework/.venv/bin/python'

def main(mode,name,cpu,cap,*extra):
    if os.environ.get('HERDR_ENV')!='1':raise RuntimeError('live Herdr context required')
    cpu=int(cpu);cap=float(cap)
    if cpu not in (6,7,8):raise ValueError('outside owned CPU allocation')
    from refinement_common import binding,save
    out=HERE/'invocations'/name;out.mkdir(parents=True,exist_ok=False)
    manager=['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus']
    runtime=json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']['python_environment_assignment']
    argv=manager+['systemd-run','--user','--wait','--pipe','--unit=p254-a55-'+name+'.service','--property=RuntimeMaxSec='+str(cap),'--property=MemoryMax=3221225472','--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity='+str(cpu),'env',runtime]+[k+'='+v for k,v in os.environ.items() if k.startswith('HERDR_')]+[k+'=1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')]+['taskset','-c',str(cpu),PYTHON,'-u',str(HERE/'native_dispatch.py'),mode,name,str(cpu),str(cap),*extra]
    save(out/'launch.json',dict(argv=argv,mode=mode,cpu=cpu,cap_seconds=cap,start_unix=time.time(),source=binding(__file__),preregistration=binding(HERE/'preregistration.json')))
    began=time.monotonic()
    with (out/'stdout').open('xb') as stdout,(out/'stderr').open('xb') as stderr:
        p=subprocess.run(argv,stdout=stdout,stderr=stderr)
    elapsed=time.monotonic()-began
    receipt=dict(name=name,mode=mode,cpu=cpu,exitcode=p.returncode,elapsed_seconds=elapsed,stdout=binding(out/'stdout'),stderr=binding(out/'stderr'))
    save(out/'exit.json',receipt)
    print(json.dumps(receipt),flush=True);print('A55_TERMINAL_'+name+'_EXIT_'+str(p.returncode),flush=True)
    return p.returncode
if __name__=='__main__':sys.exit(main(*sys.argv[1:]))
