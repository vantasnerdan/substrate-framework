"""One bounded one-core invocation, raw native banner and all-exit receipts."""
from dressed_reference import *
import time,resource,traceback,subprocess

def main(mode,name,cpu,cap,dataset,catalog,output,*extra):
    start=time.monotonic();cpu=int(cpu);cap=float(cap);receipt=HERE/'invocations'/name;unit='p254-a56-'+name+'.service'
    query=['systemctl','--user','show',unit,'--property=MainPID,ControlGroup,MemoryMax,MemorySwapMax,CPUAffinity,Restart,RuntimeMaxUSec,InvocationID,ActiveState,SubState']
    result=subprocess.run(query,capture_output=True,text=True);props=dict(line.split('=',1) for line in result.stdout.splitlines() if '=' in line);cg=Path('/sys/fs/cgroup')/props['ControlGroup'].lstrip('/')
    limits={}
    for key in ('memory.max','memory.swap.max','cpuset.cpus.effective'):
        try:limits[key]=(cg/key).read_text().strip()
        except OSError as error:limits[key]=dict(unavailable=str(error))
    banner=dict(argv=sys.argv,pid=os.getpid(),affinity=sorted(os.sched_getaffinity(0)),properties=props,kernel=limits,query_exit=result.returncode,query_stdout=result.stdout,query_stderr=result.stderr,HERDR={k:v for k,v in os.environ.items() if k.startswith('HERDR_')},threads={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')})
    banner['implementation_bindings']=[binding(p) for p in sorted(HERE.glob('*.py'))]
    save(receipt/'native-start.json',banner);print('A56_NATIVE_START '+json.dumps(banner),flush=True);code=0
    try:
        if banner['affinity']!=[cpu] or limits['memory.max']!='3221225472' or limits['memory.swap.max']!='0':raise RuntimeError('actual allocation differs')
        output=Path(output)
        if mode=='assemble':
            from evaluate_dressed import assemble
            value=assemble(dataset,catalog,output)
        else:
            import dressed_functionals as f
            if mode=='oracle':value=f.oracle(dataset,catalog,output/'boundary')
            elif mode=='ee':value=f.ee(dataset,catalog,output/('-'.join(name.split('-')[:2])),int(extra[0]),int(extra[1]),cap-8)
            elif mode=='mixed':value=f.mixed(dataset,catalog,output/'mixed')
            elif mode=='cubic':value=f.cubic(dataset,catalog,output/'cubic')
            else:raise ValueError(mode)
        save(receipt/'result.json',value)
    except BaseException:traceback.print_exc();code=1
    save(receipt/'native-exit.json',dict(mode=mode,name=name,exitcode=code,elapsed_seconds=time.monotonic()-start,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,cgroup_peak_bytes=(cg/'memory.peak').read_text().strip()))
    print('A56_NATIVE_'+name+'_EXIT_'+str(code),flush=True);return code
if __name__=='__main__':sys.exit(main(*sys.argv[1:]))
