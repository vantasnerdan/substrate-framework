"""Execute one bounded actual component and persist native startup/exit first."""
import os,sys,time,resource,traceback,subprocess,json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def main(mode,name,cpu,cap,*args):
    start=time.monotonic();cpu=int(cpu);cap=float(cap)
    receipt=HERE/'invocations'/name;receipt.mkdir(parents=True,exist_ok=True)
    from refinement_common import save,binding
    unit='p254-a55-'+name+'.service'
    query=['systemctl','--user','show',unit,'--property=MainPID,ControlGroup,MemoryMax,MemorySwapMax,CPUAffinity,Restart,RuntimeMaxUSec,InvocationID,ActiveState,SubState']
    result=subprocess.run(query,capture_output=True,text=True)
    props=dict(line.split('=',1) for line in result.stdout.splitlines() if '=' in line)
    cg=Path('/sys/fs/cgroup')/props['ControlGroup'].lstrip('/')
    limits={}
    for key in ('memory.max','memory.swap.max','cpuset.cpus.effective'):
        try:limits[key]=(cg/key).read_text().strip()
        except OSError as error:limits[key]=dict(unavailable=str(error))
    save(receipt/'native-start.json',dict(argv=sys.argv,pid=os.getpid(),affinity=sorted(os.sched_getaffinity(0)),properties=props,kernel=limits,query_argv=query,query_exit=result.returncode,query_stdout=result.stdout,query_stderr=result.stderr,HERDR={k:v for k,v in os.environ.items() if k.startswith('HERDR_')},threads={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')}))
    code=0
    try:
        if sorted(os.sched_getaffinity(0))!=[cpu] or limits['memory.max']!='3221225472' or limits['memory.swap.max']!='0':raise RuntimeError('actual owned resource cap differs')
        if mode=='reference':
            from produce_midpoint import run
            value=run(HERE/'reference')
        elif mode=='ee':
            from refine_ee import run
            value=run(name,int(args[0]),int(args[1]),cap)
        elif mode=='cubic':
            from refine_cubic import run
            value=run(HERE/'material/cubic',wall_cap=cap-8)
        elif mode=='cubic-resume':
            from resume_cubic import run
            value=run(HERE/'material/cubic-complete',wall_cap=cap-8)
        elif mode in ('ei','joint','represented'):
            import refine_ei_joint as module
            value=getattr(module,'run_'+mode)(HERE/'material'/mode)
        elif mode=='merge':
            from assemble_midpoint import run
            value=run()
        else:raise ValueError(mode)
        save(receipt/'result.json',value)
        print(json.dumps(value),flush=True)
    except BaseException:
        traceback.print_exc();code=1
    save(receipt/'native-exit.json',dict(mode=mode,name=name,exitcode=code,elapsed_seconds=time.monotonic()-start,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,cgroup_peak_bytes=(cg/'memory.peak').read_text().strip()))
    print('A55_NATIVE_'+name+'_EXIT_'+str(code),flush=True);return code
if __name__=='__main__':sys.exit(main(*sys.argv[1:]))
