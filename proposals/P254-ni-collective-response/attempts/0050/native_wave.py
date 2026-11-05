"""Bounded science stages with durable argv/raw streams/all-exit accounting."""
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
from joint_gamma_adapter import HERE
from moving_reference_inputs import atomic_json,file_binding


def main(output,mode):
    output=Path(output);output.mkdir(exist_ok=False)
    start=time.monotonic()
    startup=dict(pid=os.getpid(),affinity=sorted(os.sched_getaffinity(0)),
        cgroup=Path('/proc/self/cgroup').read_text(),status=Path('/proc/self/status').read_text(),
        environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')},
        implementation=[file_binding(p) for p in sorted(HERE.glob('*.py'))],wall_start=time.time())
    atomic_json(output/'startup.json',startup)
    if startup['affinity']!=[9] or any(v!='1' for v in startup['environment'].values()):
        raise RuntimeError('CPU9 and all single BLASthread required')
    stages=[]
    def stage(label,args,cap):
        before=resource.getrusage(resource.RUSAGE_CHILDREN);began=time.monotonic()
        argv=[sys.executable,'-u',str(HERE/args[0]),*map(str,args[1:])]
        atomic_json(output/(label+'.argv.json'),dict(argv=argv,cap_seconds=cap,start_unix=time.time()))
        timedout=False
        with (output/(label+'.stdout')).open('xb') as stdout,(output/(label+'.stderr')).open('xb') as stderr:
            try:
                result=subprocess.run(argv,stdout=stdout,stderr=stderr,timeout=cap,check=False)
                code=result.returncode
            except subprocess.TimeoutExpired:
                timedout=True;code=-9
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        receipt=dict(argv=argv,exitcode=code,timedout=timedout,elapsed_seconds=time.monotonic()-began,
            user_seconds=after.ru_utime-before.ru_utime,system_seconds=after.ru_stime-before.ru_stime,
            stdout=file_binding(output/(label+'.stdout')),stderr=file_binding(output/(label+'.stderr')))
        atomic_json(output/(label+'.exit.json'),receipt);stages.append(receipt)
        print((output/(label+'.stdout')).read_text(),end='',flush=True)
        print((output/(label+'.stderr')).read_text(),end='',file=sys.stderr,flush=True)
        print(json.dumps(receipt),flush=True)
        if code:raise RuntimeError('actual stage failure: '+label)
    try:
        if mode=='all':
            stage('oracle',['oracle_joint_gamma.py',output/'oracle.json'],60)
            stage('common',['run_joint_gamma.py',output/'material','common'],450)
            stage('represented',['run_joint_gamma.py',output/'material','represented'],890)
            stage('consumer',['consume_joint_gamma.py',output/'material/result.json',output/'consumer.json'],60)
        elif mode=='recover':
            stage('oracle',['oracle_joint_gamma.py',output/'oracle.json'],60)
            stage('represented',['run_joint_gamma.py',output/'material','represented',
                '/var/tmp/p254-a50-joint-gamma-002/material/common.npz'],890)
            stage('consumer',['consume_joint_gamma.py',output/'material/result.json',output/'consumer.json'],60)
        elif mode=='handoff':
            stage('handoff',['assemble_handoff.py','003'],60)
        elif mode=='oracle':
            stage('oracle',['oracle_joint_gamma.py',output/'oracle.json'],60)
        else:
            raise ValueError('mode all/recover/oracle/handoff required')
    finally:
        total=sum(p.stat().st_size for p in output.rglob('*') if p.is_file())
        atomic_json(output/'wave-exit.json',dict(stages=stages,cumulative_elapsed_seconds=time.monotonic()-start,
            bytes_before_receipt=total,within_output_cap=total<=128*1024**2))


if __name__=='__main__':main(sys.argv[1],sys.argv[2])
