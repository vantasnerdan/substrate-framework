"""Bounded post-terminal cache publication; no all64 replay."""
from pathlib import Path
import os,sys,subprocess,json,time,resource,traceback
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'0055'))
from refinement_common import save,binding
NAME='paid-cache-001';UNIT='p254-a57-'+NAME+'.service';RECEIPTS=HERE/'invocations'/NAME

def main():
    RECEIPTS.mkdir(parents=True,exist_ok=True)
    start=time.monotonic()
    if len(sys.argv)==1:
        prior=json.loads((HERE/'actual-common-mean-feedback-handoff-001.json').read_bytes())
        save(HERE/'paid-cache-execution-contract.json',dict(source=binding(HERE/'paid_mean_feedback_cache.py'),runner=binding(__file__),input_result=prior['result'],prior_cumulative_native_seconds=prior['costs']['cumulative_native_launcher_seconds'],plan='Export paid full force plus old Gamma term and Phi000 from unchanged actual result/source after completed Native0. Only internalq0 contracted; no nonGamma replay, no new symmetry proof. Exercise saved cache on original actual dataset, require exact full-force byte recovery.',caps=dict(cpu=10,blas=1,memory_bytes=3221225472,swap_bytes=0,native_seconds=120,new_total_payload_bytes=67108864)))
        runtime=json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
        argv=runtime['manager_env']+['systemd-run','--user','--wait','--pipe','--unit='+UNIT,'--property=CPUAffinity=10','--property=MemoryMax=3G','--property=MemorySwapMax=0','--property=Restart=no','--property=RuntimeMaxSec=120','env',runtime['python_environment_assignment']]+[k+'='+v for k,v in os.environ.items() if k.startswith('HERDR_')]+[k+'=1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')]+[runtime['python'],'-u',str(Path(__file__).resolve()),'native']
        save(RECEIPTS/'invocation.json',dict(argv=argv,started_unix=time.time()))
        with (RECEIPTS/'stdout-stderr.log').open('w') as log:
            p=subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            for line in p.stdout:log.write(line);log.flush();print(line,end='',flush=True)
            code=p.wait()
        save(RECEIPTS/'launcher-exit.json',dict(exitcode=code,elapsed_seconds=time.monotonic()-start,output=binding(RECEIPTS/'stdout-stderr.log')))
        print('A57_PAID_CACHE_ALL_EXIT_'+str(code),flush=True);return code
    props=subprocess.run(['systemctl','--user','show',UNIT,'--property=MainPID,ControlGroup,MemoryMax,MemorySwapMax,CPUAffinity,Restart,RuntimeMaxUSec,InvocationID'],capture_output=True,text=True)
    fields=dict(line.split('=',1) for line in props.stdout.splitlines() if '=' in line);cg=Path('/sys/fs/cgroup')/fields['ControlGroup'].lstrip('/')
    limits={k:(cg/k).read_text().strip() for k in ('memory.max','memory.swap.max')}
    banner=dict(argv=sys.argv,pid=os.getpid(),affinity=sorted(os.sched_getaffinity(0)),properties=fields,kernel=limits,HERDR={k:v for k,v in os.environ.items() if k.startswith('HERDR_')},threads={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')})
    save(RECEIPTS/'native-start.json',banner);print(json.dumps(dict(stage='A57_PAID_CACHE_NATIVE_START',**banner)),flush=True)
    code=0
    try:
        if banner['affinity']!=[10] or limits!={'memory.max':'3221225472','memory.swap.max':'0'}:raise RuntimeError('owned resource cap mismatch')
        from paid_mean_feedback_cache import export_paid_mean_cache
        result=export_paid_mean_cache(HERE/'material/common-mean-001/reference-result.json',HERE/'material/paid-cache-001')
        save(RECEIPTS/'result.json',result);print(json.dumps(result),flush=True)
    except BaseException:traceback.print_exc();code=1
    save(RECEIPTS/'native-exit.json',dict(exitcode=code,elapsed_seconds=time.monotonic()-start,maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,cgroup_peak_bytes=(cg/'memory.peak').read_text().strip()))
    print('A57_PAID_CACHE_NATIVE_EXIT_'+str(code),flush=True);return code
if __name__=='__main__':raise SystemExit(main())
