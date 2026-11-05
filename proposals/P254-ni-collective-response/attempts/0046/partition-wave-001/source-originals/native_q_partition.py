"""Root565 exact outer-q partition worker; no copied material banks.

One q owns all128externals/all16pairs/all64innerh. Recovery granularity is
ONLY committed outerq; currentq can be lost (~4-6min expected), never a
claimed fine-grained withinq recovery. Prefix storage remains two generations.
"""
import argparse
import contextlib
import datetime
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time
import traceback
HERE=Path(__file__).resolve().parent


class SafePause(Exception): pass


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--name',choices=('A','B','C'),required=True)
    p.add_argument('--q-start',type=int,required=True); p.add_argument('--q-stop',type=int,required=True)
    p.add_argument('--cpu',type=int,required=True); p.add_argument('--output',type=Path,required=True)
    p.add_argument('--resume',action='store_true')
    a=p.parse_args()
    out=a.output.resolve(); out.mkdir(parents=True,exist_ok=True)
    unit='p254-a46-partition-'+a.name.lower()+'-001.service'
    start=time.monotonic()
    def save(name,value):
        from all_slot_coulomb_checkpoint import atomic_json
        atomic_json(out/name,value)
    show=['systemctl','--user','show',unit,'--property=MainPID,ControlGroup,MemoryMax,MemorySwapMax,Restart,RuntimeMaxUSec,CPUAffinity,ActiveState,SubState,InvocationID']
    m=subprocess.run(show,capture_output=True,text=True)
    props=dict(line.split('=',1) for line in m.stdout.splitlines() if '=' in line)
    cg=Path('/sys/fs/cgroup')/props['ControlGroup'].lstrip('/')
    kernel={}
    for name in ('memory.max','memory.swap.max','cpuset.cpus.effective'):
        try: kernel[name]=dict(available=True,value=(cg/name).read_text().strip())
        except OSError as e: kernel[name]=dict(available=False,reason=str(e))
    save('actual-startup.json',dict(unit=unit,argv=sys.argv,manager_argv=show,
        manager_exit=m.returncode,manager_stdout=m.stdout,manager_stderr=m.stderr,
        properties=props,kernel=kernel,affinity=sorted(os.sched_getaffinity(0)),
        utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    if sorted(os.sched_getaffinity(0))!=[a.cpu]:
        raise RuntimeError('actual partition CPU affinity differs')
    code,status=0,'complete'
    def deadline(sig,frame):
        raise SafePause('bounded owned deadline; durable completeq prefixes retained, current uncommittedq may replay')
    signal.signal(signal.SIGALRM,deadline)
    signal.setitimer(signal.ITIMER_REAL,780)
    with (out/'science.rawstdout.txt').open('w',buffering=1) as stdout,(out/'science.rawstderr.txt').open('w',buffering=1) as stderr:
        with contextlib.redirect_stdout(stdout),contextlib.redirect_stderr(stderr):
            try:
                import numpy as np
                import run_all_slot_coulomb as run
                from all_slot_coulomb_checkpoint import AllSlotPrefix,binding,consume_result
                from all_slot_coulomb_partition import partition_specification,evaluate_partition
                from all_slot_coulomb_shared import bounds
                from material_electron_ion_inputs import SavedTransitions
                owned=list(range(a.q_start,a.q_stop+1))
                c=run.load_coulomb_inputs(density_cap_bytes=900*1024**2)
                times=np.array([0.,1.,10.,30.])
                bundle='/var/tmp/p254-a45-moving-reference-001/producer/bundle.json'
                reference=SavedTransitions(c,times,bundle,c.bindings['manifest'])
                common=run.specification(c,reference,times,2,[binding(__file__)])
                common['bounds']=dict(bounds(external_batch=2,time_batch=4),
                    conservative_worker_ram_bytes=3*1024**3,global_output_cap_bytes=768*1024**2,
                    time_batch=4,external_batch=2,interruption_loss='current outerq only, no withinq history')
                common['partition_wave']='Root565 exact all64q shared-superoperator wave'
                spec=partition_specification(common,owned)
                save('source-input-specification.json',spec)
                result_dir=out/'result'
                if result_dir.exists() and not a.resume:
                    raise FileExistsError('partition output exists; explicit exact-source --resume required')
                result_dir.mkdir(exist_ok=True)
                electrons={(ti,si):reference.electron(ti,si) for ti in range(4) for si in range(4)}
                last=[time.monotonic()]
                progress=[]
                def committed(count):
                    now=time.monotonic(); seconds=now-last[0]; last[0]=now
                    cache=AllSlotPrefix(result_dir/'numeric-prefix',spec,spec['shape'],len(owned))
                    saved,restored=cache.restore()
                    if restored!=count: raise ValueError('actual committedq consumer count differs')
                    norms={k:float(np.linalg.norm(x)) for k,x in saved.items()}
                    side=json.loads((result_dir/'numeric-prefix'/f'prefix-{count%2}.json').read_bytes())
                    remaining_q=owned[count:]
                    feature_counts=[len(g) for g in c.feature_g]
                    forecast=seconds*sum(feature_counts[q] for q in remaining_q)/feature_counts[owned[count-1]]
                    event=dict(partition=a.name,completed_count=count,completed_global_q=owned[:count],
                        just_completed_global_q=owned[count-1],currentq_seconds=seconds,
                        elapsed_native_seconds=now-start,remaining_partition_seconds_forecast=forecast,
                        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                        actual_cgroup_memory_bytes=int((cg/'memory.current').read_text()),
                        numeric_prefix=side['archive'],decoded_norms=norms,
                        approximation='raw original A45 movingU, no repair; not globalaccuracycertification',
                        contribution_scope='PARTIAL_Q_CONTRIBUTION')
                    progress.append(event); save('actual-q-progress.json',progress)
                    print(json.dumps(event),flush=True)
                    print(json.dumps(event),file=sys.__stdout__,flush=True)
                    if count==1:
                        save('first-complete-q.json',event)
                        print('A46_PARTITION_'+a.name+'_FIRST_Q_COMPLETE',flush=True)
                        print('A46_PARTITION_'+a.name+'_FIRST_Q_COMPLETE',file=sys.__stdout__,flush=True)
                        signal.setitimer(signal.ITIMER_REAL,max(1.,8900-(time.monotonic()-start)))
                    if event['actual_cgroup_memory_bytes']>int(2.8*1024**3):
                        raise SafePause('own numeric allocations near3GiB; pause after committedq, no duplicateworker')
                result,counts=evaluate_partition(run.BareCoulombFactors(c),electrons,times,
                    result_dir,spec,after_commit=committed)
                save('complete-partial-result-consumption.json',dict(status='PARTIAL_Q_CONTRIBUTION',
                    archive=result['record']['archive'],counts=counts,
                    decoded_norms=result['decoded_norms'],record=binding(result_dir/'result.json')))
                reference.close()
            except SafePause as e:
                code,status=2,'checkpoint-safe-paused'
                print(str(e),file=sys.stderr,flush=True)
            except BaseException:
                code,status=1,'failed'
                traceback.print_exc()
            finally:
                signal.setitimer(signal.ITIMER_REAL,0)
    terminal=dict(unit=unit,partition=a.name,status=status,exit=code,seconds=time.monotonic()-start,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    save('actual-native-exit.json',terminal)
    print(json.dumps(terminal),flush=True)
    return code


if __name__=='__main__': sys.exit(main())
