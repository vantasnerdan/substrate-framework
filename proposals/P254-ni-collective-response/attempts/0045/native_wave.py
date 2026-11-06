"""One bounded Native wave, durable original argv/stdout/stderr/all-exit receipts."""
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
from moving_reference_inputs import HERE, ROOT, CONTRACT, atomic_json, file_binding, own_sources


def main(output):
    output=Path(output);output.mkdir(exist_ok=False)
    start=time.monotonic()
    startup=dict(pid=os.getpid(),affinity=sorted(os.sched_getaffinity(0)),
        cgroup=Path('/proc/self/cgroup').read_text(),status=Path('/proc/self/status').read_text(),
        environment={key:os.environ.get(key) for key in ('PYTHONPATH','OPENBLAS_NUM_THREADS',
        'OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')},
        contract=CONTRACT,implementation=own_sources(),wall_start=time.time())
    atomic_json(output/'actual-startup.json',startup)
    if startup['affinity'] != [7]:
        raise RuntimeError('CPU7 singleton affinity required')
    manager=subprocess.run(['env','XDG_RUNTIME_DIR=/run/user/1000',
        'DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus','systemctl','--user','show',
        'p254-a45-moving-wave-001.service','--property=Id,InvocationID,MainPID,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,CPUAffinity,WatchdogUSec,ActiveState,SubState'],
        capture_output=True,check=False)
    atomic_json(output/'manager-startup.json',dict(argv=manager.args,exitcode=manager.returncode,
        stdout=manager.stdout.decode(),stderr=manager.stderr.decode()))
    stages=[]
    def stage(label,argv,cap):
        before=resource.getrusage(resource.RUSAGE_CHILDREN)
        launched=time.monotonic();wall=time.time()
        cap=min(cap,1785-(launched-start))
        if cap<=0:raise RuntimeError('whole1800 budget exhausted before '+label)
        atomic_json(output/(label+'-argv.json'),dict(argv=argv,start_unix=wall,hard_cap_seconds=cap,
            cumulative_native_elapsed_before=launched-start))
        timedout=False
        with (output/(label+'-stdout.txt')).open('wb') as stdout, (output/(label+'-stderr.txt')).open('wb') as stderr:
            try:
                result=subprocess.run(argv,stdout=stdout,stderr=stderr,timeout=cap,check=False)
                exitcode=result.returncode
            except subprocess.TimeoutExpired:
                timedout=True;exitcode=-9
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        receipt=dict(label=label,argv=argv,start_unix=wall,exitcode=exitcode,timedout=timedout,
            elapsed_seconds=time.monotonic()-launched,cumulative_native_seconds=time.monotonic()-start,
            child_user_seconds=after.ru_utime-before.ru_utime,child_system_seconds=after.ru_stime-before.ru_stime,
            stdout=file_binding(output/(label+'-stdout.txt')),stderr=file_binding(output/(label+'-stderr.txt')))
        atomic_json(output/(label+'-exit.json'),receipt);stages.append(receipt)
        print(json.dumps(receipt),flush=True)
        if exitcode:
            atomic_json(output/'wave-result.json',dict(status='failed',stages=stages,
                failure=label,cumulative_native_seconds=time.monotonic()-start))
            raise RuntimeError('Native stage failed: '+label)
        return receipt
    stage('new-oracle',[sys.executable,'-u',str(HERE/'oracle_moving_reference.py'),str(output/'oracle')],60)
    atomic_json(output/'first-producer-allocation.json',dict(decision='continue authorized first actual producer',
        expected_first_useful_checkpoint='actual U(1,0), expected <=120s including original authentication and affine h',
        expected_producer_seconds=300,hard_cap_seconds=600,whole_budget_seconds=1800,
        trigger='3x first-checkpoint expectation (360s) requires stop/change allocation, not unchanged renewal',
        remaining_budget_seconds=1800-(time.monotonic()-start),wall_unix=time.time()))
    stage('actual-producer',[sys.executable,'-u',str(HERE/'run_moving_reference.py'),'produce',
        '--output',str(output/'producer')],600)
    stage('actual-moving-EI',[sys.executable,'-u',str(HERE.parent/'0041/run_material_electron_ion_kernel.py'),
        '--output',str(output/'moving-EI'),'--saved-transitions',str(output/'producer/bundle.json'),
        '--sigma-slots','0,64','--pi-q','0,1,2'],1200)
    stage('actual-kernel-recovery',[sys.executable,'-u',str(HERE/'run_moving_reference.py'),'recover',
        '--kernel-output',str(output/'moving-EI'),'--bundle',str(output/'producer/bundle.json'),
        '--output',str(output/'kernel-recovery.json')],120)
    files=[file_binding(path) for path in sorted(output.rglob('*')) if path.is_file()]
    total=sum(item['bytes'] for item in files)
    if total>402653184:
        raise RuntimeError('combined output exceeded authorized384MiB')
    atomic_json(output/'wave-result.json',dict(schema='p254-a45-native-wave-v1',status='complete',
        stages=stages,cumulative_native_seconds=time.monotonic()-start,
        total_output_bytes_before_final_receipt=total,files=files,
        scope='genuine original moving-mean U plus original bare S and conditional moving-reference E-I kernels only',
        findings=[],blockers=[]))
    print('[dft-done] A45 original moving-reference wave complete',flush=True)

if __name__=='__main__':main(sys.argv[1])
