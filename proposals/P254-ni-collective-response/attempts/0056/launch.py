"""Actual single-thread systemd launcher with immutable per-invocation receipts."""
from dressed_reference import *
import subprocess,time
PYTHON='/home/administrator/projects/substrate-framework/.venv/bin/python'
def main(mode,name,cpu,cap,dataset,catalog,output,*extra):
    if os.environ.get('HERDR_ENV')!='1':raise RuntimeError('actual Herdr context required')
    cpu=int(cpu);cap=float(cap)
    if cpu not in (6,7,8):raise ValueError('not owned CPU')
    receipt=HERE/'invocations'/name;receipt.mkdir(parents=True,exist_ok=False)
    runtime=json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']['python_environment_assignment']
    argv=['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus','systemd-run','--user','--wait','--pipe','--unit=p254-a56-'+name+'.service','--property=RuntimeMaxSec='+str(cap),'--property=MemoryMax=3221225472','--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity='+str(cpu),'env',runtime]+[k+'='+v for k,v in os.environ.items() if k.startswith('HERDR_')]+[k+'=1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')]+['taskset','-c',str(cpu),PYTHON,'-u',str(HERE/'native.py'),mode,name,str(cpu),str(cap),dataset,catalog,output,*extra]
    save(receipt/'launch.json',dict(argv=argv,mode=mode,cpu=cpu,cap_seconds=cap,start_unix=time.time(),source=binding(__file__),preregistration=binding(HERE/'preregistration.json')))
    began=time.monotonic()
    with (receipt/'stdout').open('xb') as out,(receipt/'stderr').open('xb') as err:
        p=subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=err)
        for line in p.stdout:
            out.write(line);out.flush();print(line.decode(errors='replace'),end='',flush=True)
        p.wait()
    value=dict(name=name,mode=mode,cpu=cpu,exitcode=p.returncode,elapsed_seconds=time.monotonic()-began,stdout=binding(receipt/'stdout'),stderr=binding(receipt/'stderr'))
    save(receipt/'exit.json',value);print(json.dumps(value),flush=True);print('A56_TERMINAL_'+name+'_EXIT_'+str(p.returncode),flush=True);return p.returncode
if __name__=='__main__':sys.exit(main(*sys.argv[1:]))
