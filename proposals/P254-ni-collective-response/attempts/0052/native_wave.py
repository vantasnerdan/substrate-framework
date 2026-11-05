"""A52 one-thread CPU7 durable prerequisite/stage/exit evidence."""
import os,json,sys,time,resource,subprocess
from pathlib import Path
from general_moving_reference import HERE,atomic_json,file_binding,authenticate

def main(output):
    output=Path(output);output.mkdir(exist_ok=False);start=time.monotonic();stages=[]
    startup=dict(pid=os.getpid(),affinity=sorted(os.sched_getaffinity(0)),cgroup=Path('/proc/self/cgroup').read_text(),
        status=Path('/proc/self/status').read_text(),environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')},
        sources=[file_binding(p) for p in sorted(HERE.iterdir()) if p.is_file()],start_unix=time.time())
    atomic_json(output/'startup.json',startup)
    if startup['affinity']!=[7] or any(v!='1' for v in startup['environment'].values()):raise RuntimeError('actual CPU7/oneThread required')
    handoff_path=HERE.parent/'0045/actual-complete-original-handoff-001.json'
    handoff=json.loads(handoff_path.read_bytes());original=handoff['original_input_bindings']
    bindings=[file_binding(handoff_path),original['baseline'],original['material_manifest'],*original['settled_sources'],*handoff['source_bindings'],handoff['outputs']['bundle'],handoff['outputs']['raw_transitions']]
    for binding in bindings:authenticate(binding)
    inventory=dict(schema='p254-a52-authenticated-prerequisite-inventory-v1',bindings=bindings,
        source_inventory=startup['sources'],preregistration=file_binding(HERE/'implementation-contract-001.json'),
        original_initial_policy='unchanged raw original coherentCAR and prepared Gaussian; all banks read-only',
        carried_proofs=handoff['unaffected_carried_proofs'],authenticated_before_new_oracle=True)
    atomic_json(output/'input-inventory.json',inventory)
    def stage(label,args,cap):
        began=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
        argv=[sys.executable,'-u',str(HERE/args[0]),*map(str,args[1:])]
        atomic_json(output/(label+'.argv.json'),dict(argv=argv,cap_seconds=cap,start_unix=time.time(),source=file_binding(HERE/args[0]),prerequisite=file_binding(output/'input-inventory.json')))
        timedout=False
        with (output/(label+'.stdout')).open('xb') as out,(output/(label+'.stderr')).open('xb') as err:
            try:code=subprocess.run(argv,stdout=out,stderr=err,timeout=cap,check=False).returncode
            except subprocess.TimeoutExpired:timedout=True;code=-9
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        receipt=dict(argv=argv,exitcode=code,timedout=timedout,elapsed_seconds=time.monotonic()-began,user_seconds=after.ru_utime-before.ru_utime,system_seconds=after.ru_stime-before.ru_stime,stdout=file_binding(output/(label+'.stdout')),stderr=file_binding(output/(label+'.stderr')))
        atomic_json(output/(label+'.exit.json'),receipt);stages.append(receipt)
        print((output/(label+'.stdout')).read_text(),end='',flush=True);print((output/(label+'.stderr')).read_text(),end='',file=sys.stderr,flush=True)
        print(json.dumps(receipt),flush=True)
        if code:raise RuntimeError('actual stage failure: '+label)
    try:
        stage('oracle',['oracle_requested_reference.py',output/'oracle.json'],30)
        stage('producer',['run_requested_reference.py',output/'producer'],140)
        stage('consumer',['consume_requested_reference.py',output/'producer/bundle.json',output/'consumer.json'],110)
    finally:
        total=sum(p.stat().st_size for p in output.rglob('*') if p.is_file())
        record=dict(stages=stages,elapsed_seconds=time.monotonic()-start,output_bytes_before_receipt=total,output_cap_bytes=64*1024**2,within_cap=total<64*1024**2)
        atomic_json(output/'wave-exit.json',record)
        if total>=64*1024**2:raise RuntimeError('hard new payload cap exceeded')
if __name__=='__main__':main(sys.argv[1])
