"""One actual CPU7/3GiB/noSwap bounded invocation, cumulative all-attempt ledger."""
import json,os,sys,time,subprocess
from pathlib import Path
from control import HERE,binding,save
PYTHON='/home/administrator/projects/substrate-framework/.venv/bin/python'
BASE=Path('/var/tmp/p254-a52-requested-reference')
def main(generation,mode='science'):
    if os.environ.get('HERDR_ENV')!='1':raise RuntimeError('actual HERDR_ENV=1 required')
    receipts=Path(str(BASE)+'-'+generation+'-receipts');receipts.mkdir(exist_ok=False)
    previous=sum(json.loads(p.read_bytes())['elapsed_seconds'] for p in BASE.parent.glob(BASE.name+'-*-receipts/exit.json'))
    remaining=600-previous;cap=min(300,int(remaining))
    if cap<20:raise RuntimeError('cumulative Native allocation exhausted')
    output=Path(str(BASE)+'-'+generation);unit='p254-a52-reference-'+generation+'.service'
    contract=json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())
    manager=['env','XDG_RUNTIME_DIR=/run/user/1000','DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus']
    script='native_wave.py' if mode=='science' else 'assemble_handoff.py'
    arguments=[str(output)] if mode=='science' else []
    inherited_herdr=[key+'='+value for key,value in sorted(os.environ.items()) if key.startswith('HERDR_')]
    argv=manager+['systemd-run','--user','--wait','--pipe','--unit='+unit,'--property=RuntimeMaxSec='+str(cap),'--property=MemoryMax=3221225472','--property=MemorySwapMax=0','--property=Restart=no','--property=CPUAffinity=7','env',contract['common_runtime']['python_environment_assignment'],*inherited_herdr]+[key+'=1' for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')]+['taskset','-c','7',PYTHON,'-u',str(HERE/script),*arguments]
    began=time.monotonic();save(receipts/'argv.json',dict(argv=argv,start_unix=time.time(),HERDR_ENV=os.environ['HERDR_ENV'],prior_native_seconds=previous,remaining_seconds=remaining,source=binding(HERE/'launch_native.py'),preregistration=binding(HERE/'implementation-contract-001.json')))
    with (receipts/'stdout').open('xb') as out,(receipts/'stderr').open('xb') as err:
        proc=subprocess.Popen(argv,stdout=out,stderr=err);time.sleep(1)
        query=manager+['systemctl','--user','show',unit,'--property=Id,InvocationID,MainPID,MemoryMax,MemorySwapMax,CPUAffinity,Restart,RuntimeMaxUSec,ActiveState,SubState']
        save(receipts/'manager-query-argv.json',dict(argv=query,start_unix=time.time()))
        result=subprocess.run(query,capture_output=True,check=False)
        (receipts/'manager-startup.stdout').write_bytes(result.stdout);(receipts/'manager-startup.stderr').write_bytes(result.stderr)
        save(receipts/'manager-startup.json',dict(argv=query,exitcode=result.returncode,stdout=binding(receipts/'manager-startup.stdout'),stderr=binding(receipts/'manager-startup.stderr')))
        code=proc.wait()
    elapsed=time.monotonic()-began
    save(receipts/'exit.json',dict(argv=argv,exitcode=code,elapsed_seconds=elapsed,cumulative_native_seconds=previous+elapsed,stdout=binding(receipts/'stdout'),stderr=binding(receipts/'stderr')))
    if mode=='finalize':
        handoff=HERE/'actual-original-uniform-moving-reference-handoff-001.json'
        runs=[json.loads(p.read_bytes()) for p in sorted(BASE.parent.glob(BASE.name+'-*-receipts/exit.json'))]
        evidence=[binding(p) for p in sorted(receipts.iterdir()) if p.is_file()]
        save(HERE/'actual-final-native-binding-002.json',dict(schema='p254-a52-final-terminal-binding-v1',status='complete' if code==0 else 'failed',handoff=binding(handoff) if handoff.exists() else None,all_native_runs=runs,cumulative_native_seconds=sum(r['elapsed_seconds'] for r in runs),final_native_evidence=evidence,source_bindings=[binding(p) for p in sorted(HERE.glob('*.py'))],limits='uniform reference only, not uniform selfenergy or closure'))
    print((receipts/'stdout').read_text(),end='',flush=True);print((receipts/'stderr').read_text(),end='',file=sys.stderr,flush=True)
    print('A52_TERMINAL_'+generation+'_EXIT_'+str(code),flush=True);return code
if __name__=='__main__':sys.exit(main(sys.argv[1],sys.argv[2] if len(sys.argv)>2 else 'science'))
