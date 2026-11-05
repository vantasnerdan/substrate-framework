"""One owned bounded native invocation with original argv/stdout/stderr/exit."""
import os,sys,json,time,resource,subprocess,traceback
from pathlib import Path
from mean_feedback_functional import HERE,binding,save

UNIT='p254-a57-common-mean-002.service'

def main():
    start=time.monotonic();receipts=HERE/'invocations/common-mean-002';receipts.mkdir(parents=True,exist_ok=True)
    query=['systemctl','--user','show',UNIT,'--property=MainPID,ControlGroup,MemoryMax,MemorySwapMax,CPUAffinity,Restart,RuntimeMaxUSec,InvocationID,ActiveState,SubState']
    props=subprocess.run(query,capture_output=True,text=True)
    fields=dict(line.split('=',1) for line in props.stdout.splitlines() if '=' in line)
    cg=Path('/sys/fs/cgroup')/fields['ControlGroup'].lstrip('/')
    limits={k:(cg/k).read_text().strip() if (cg/k).exists() else 'unavailable-not-delegated' for k in ('memory.max','memory.swap.max','cpuset.cpus.effective')}
    banner=dict(argv=sys.argv,pid=os.getpid(),affinity=sorted(os.sched_getaffinity(0)),properties=fields,kernel=limits,query_argv=query,query_exit=props.returncode,query_stdout=props.stdout,query_stderr=props.stderr,HERDR={k:v for k,v in os.environ.items() if k.startswith('HERDR_')},threads={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')})
    save(receipts/'native-start.json',banner);print(json.dumps(dict(stage='A57_NATIVE_START',**banner)),flush=True)
    code=0
    try:
        if banner['affinity']!=[10] or limits['memory.max']!='3221225472' or limits['memory.swap.max']!='0':raise RuntimeError('actual owned cap differs')
        from produce_mean_feedback import produce_mean_feedback
        value=produce_mean_feedback('/var/tmp/p254-a53-fixed-kernel-material-003/actual-material-step.npz',HERE.parent/'0053/actual-input-catalog-003.json',HERE/'material/common-mean-001')
        save(receipts/'result.json',value)
    except BaseException:
        traceback.print_exc();code=1
    exitrec=dict(exitcode=code,elapsed_seconds=time.monotonic()-start,maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,cgroup_peak_bytes=(cg/'memory.peak').read_text().strip())
    save(receipts/'native-exit.json',exitrec)
    print('A57_NATIVE_COMMON_MEAN_002_EXIT_'+str(code),flush=True)
    return code

if __name__=='__main__':sys.exit(main())
