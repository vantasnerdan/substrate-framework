"""Native launch receipt and live terminal output; no detached hidden work."""
import os,subprocess,time,json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'0055'))
from refinement_common import binding,save

def main():
    receipts=HERE/'invocations/common-mean-002';receipts.mkdir(parents=True,exist_ok=True)
    source_paths=[HERE/n for n in ('mean_feedback_functional.py','mean_feedback_reference.py','produce_mean_feedback.py','mean_feedback_oracle.py','native_feedback.py','launch_feedback.py')]
    save(HERE/'execution-contract-002.json',dict(schema='p254-a57-common-feedback-preregistration-v1',source_inventory=[binding(p) for p in source_paths],
        input_dataset=binding('/var/tmp/p254-a53-fixed-kernel-material-003/actual-material-step.npz'),baseline_catalog=binding(HERE.parent/'0053/actual-input-catalog-003.json'),
        authoritative_cubic=binding('/var/tmp/p254-a47-full37-source-001/source.json'),authoritative_cubic_arrays=binding('/var/tmp/p254-a47-full37-source-001/source-assembly/assembled.npz'),
        times_inverse_ry=[0.,.5,1.],physical_composition='A29 scalar derivative with original mean/rho0/W0, contacts/passive complement/active K; A47 Gaussian E3 full64 covariance force, Phi000 coherent force and phase curvature; no field cubic, no XC counterterm, no represented-bubble contact subtraction.',
        new_oracle='Once after implementation settles: actual raw-density/Gamma ordered boundary, complex scalar/mean/Hamiltonian JVP, same-source cubic tensor/gradient/Hessian/covariance derivative, sign/K/conjugation sensitivities. Carry admitted A29/A47 proofs; no broad replay.',
        feature='Actual node-force contractions once per q/node then simultaneous genuine x/U/S DOP853 integration with order1 forcing only; save full original phase mean/U/GammaS and exercise persisted electron/fullGamma provider without reintegration.',
        resource_budget=dict(cpu=10,blas_threads=1,memory_max_bytes=3221225472,swap_max_bytes=0,cumulative_native_seconds=1200,new_payload_max_bytes=67108864),
        verdict='Scoped actual bounded common feedback provider and raw defects/sensitivity only; temporal forcing error unresolved, original nonGamma bath and physical calibration fixed; no conserving closed-system/selfconsistent-convergence/T30/rate/lifetime/heat claim.'))
    runtime=json.loads((HERE.parent/'0045/implementation-contract-001.json').read_bytes())['common_runtime']
    argv=runtime['manager_env']+['systemd-run','--user','--unit=p254-a57-common-mean-002','--wait','--pipe','--property=CPUAffinity=10','--property=MemoryMax=3G','--property=MemorySwapMax=0','--property=Restart=no','--property=RuntimeMaxSec=1196','/usr/bin/env',runtime['python_environment_assignment']]+[k+'='+v for k,v in os.environ.items() if k.startswith('HERDR_')]+['OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1','BLIS_NUM_THREADS=1','NUMEXPR_NUM_THREADS=1',runtime['python'],'-u',str(HERE/'native_feedback.py')]
    start=time.monotonic();save(receipts/'invocation.json',dict(argv=argv,started_unix=time.time(),launcher_pid=os.getpid()))
    with (receipts/'stdout-stderr.log').open('w') as log:
        process=subprocess.Popen(argv,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        for line in process.stdout:
            log.write(line);log.flush();print(line,end='',flush=True)
        code=process.wait()
    save(receipts/'launcher-exit.json',dict(exitcode=code,elapsed_seconds=time.monotonic()-start,output=binding(receipts/'stdout-stderr.log')))
    print('A57_ALL_EXIT_002_'+str(code),flush=True)
    return code
if __name__=='__main__':raise SystemExit(main())
