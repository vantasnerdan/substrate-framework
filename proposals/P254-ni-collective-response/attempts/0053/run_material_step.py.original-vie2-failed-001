"""Actual fixed-kernel intrinsic Ni/H Dyson prefix; no new kernel evaluations."""
from pathlib import Path
import argparse
import json
import os
import resource
from time import perf_counter
import numpy as np
from material_step_inputs import load_step, file_binding
from joint_gamma_adapter import ordered_coordinates
from supplied_vie2 import solve_supplied, load_solution
from consume_material_step import consume, json_value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalog',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    out = args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    started = perf_counter()
    print(json.dumps(dict(stage='actual-material-native-start',pid=os.getpid(),
        affinity=sorted(os.sched_getaffinity(0)),cgroup=Path('/proc/self/cgroup').read_text(),
        threads={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')})),flush=True)
    c, reference, catalog, sigma, pi, represented = load_step(args.catalog)
    try:
        times = reference.times;nt=len(times);order=min(5,nt-1)
        e_ret=np.zeros((nt,nt,c.ns,c.m,c.m),complex);e_less=np.zeros_like(e_ret)
        for a in range(nt):
            for b in range(nt):
                greater,lesser=reference.electron(a,b)
                e_less[a,b]=lesser
                if a>=b:e_ret[a,b]=greater-lesser
        original_e_less0=e_less[0,0].copy()
        binary=catalog['solver_binary']['path']
        for slot in range(c.ns):
            path=solve_supplied(times,e_ret[:,:,slot],e_less[:,:,slot],
                sigma['retarded'][:,:,slot],sigma['lesser'][:,:,slot],
                statistics=-1,order=order,binary=binary,
                workdir=out/'fermions'/('slot-%03d'%slot))
            solved=load_solution(path)
            e_ret[:,:,slot],e_less[:,:,slot]=solved['ret'],solved['less']
        d=c.d+c.f
        b_ret=np.zeros((nt,nt,d,d),complex);b_less=np.zeros_like(b_ret)
        for a in range(nt):
            for b in range(nt):
                greater,lesser=ordered_coordinates(reference,a,b)
                b_less[a,b]=lesser
                if a>=b:b_ret[a,b]=greater-lesser
        outputs={}
        for label,kernels in [('intrinsic',pi),('represented-control',represented)]:
            path=solve_supplied(times,b_ret,b_less,kernels['retarded'],kernels['lesser'],
                statistics=1,order=order,binary=binary,workdir=out/label)
            outputs[label]=load_solution(path)
        dataset=out/'actual-material-step.npz'
        np.savez(dataset,times=times,n=c.n,electron_ret=e_ret,electron_less=e_less,
            source_electron_less0=original_e_less0,
            gamma_ret=outputs['intrinsic']['ret'],gamma_less=outputs['intrinsic']['less'],
            control_ret=outputs['represented-control']['ret'],control_less=outputs['represented-control']['less'],
            source_gamma_ret=b_ret,source_gamma_less=b_less)
        # Consumer runs after all numeric outputs are persisted.
        readout=consume(dataset)
        (out/'readout.json').write_text(json.dumps(readout,default=json_value,indent=2)+'\n')
        result=dict(schema='p254-actual-fixed-kernel-material-step-v1',status='complete',
            input_catalog=file_binding(args.catalog),dataset=file_binding(dataset),
            readout=file_binding(out/'readout.json'),times_inverse_ry=times.tolist(),
            actual_solved_fermion_blocks=c.ns,actual_Gamma_coordinates=d,
            actual_existing_solver_invocations=c.ns+2,new_physical_kernel_evaluations=0,
            fixed_kernel=True,selfconsistent=False,mean='original prescribed moving mean',
            cutover='Bare mechanical source keeps contacts and original covariance; common electronic Pi replaces represented active response/noise, never stacks it. Separate represented-only control.',
            cubic='Full111 ionic Pi embedded in original114 Gamma; no field-cubic invention.',
            source_reference_bindings=reference.bindings,
            elapsed_seconds=perf_counter()-started,
            maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            limits=['Original0/1 uniformprefix only, not full0..30history.',
                    'No mean/kernel selfconsistent feedback or fullspatial bath update.',
                    'No physicalcalibration/XCdoublecountcure/lifetime/rate/heat/fullEclaim.'])
        (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(dict(stage='actual-material-step-complete',result=file_binding(out/'result.json'),readout=readout),default=json_value),flush=True)
    finally:
        reference.close()


if __name__=='__main__':main()
