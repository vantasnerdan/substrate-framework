"""Actual original all12831node U, coefficients, and raw four-node comparison."""
import sys,json
from pathlib import Path
from time import perf_counter
import numpy as np
from general_moving_reference import (HERE,inputs,AffineSavedHamiltonian,produce_u,
    SETTINGS,atomic_json,file_binding,authenticate)

def main(output):
    output=Path(output);output.mkdir(exist_ok=False);began=perf_counter()
    handoff_path=HERE.parent/'0045/actual-complete-original-handoff-001.json'
    handoff=json.loads(handoff_path.read_bytes())
    c,baseline,coefficient_binding,history=inputs()
    h=AffineSavedHamiltonian(baseline)
    times=np.arange(31,dtype=float)
    u,work=produce_u(h,times,np.broadcast_to(np.eye(c.m),(c.ns,c.m,c.m)))
    upath=output/'U.npy'
    with upath.open('xb') as f:np.save(f,u,allow_pickle=False)
    coeff=output/'hamiltonian-chebyshev.npz'
    with coeff.open('xb') as f:np.savez_compressed(f,coefficients=np.asarray(h.coefficients),segment_times=h.times)
    old_binding=handoff['outputs']['raw_transitions'];oldpath=authenticate(old_binding)
    with np.load(oldpath,allow_pickle=False) as old:
        old_u=old['U'];comparisons=[]
        for index,t in enumerate((0,1,10,30)):
            delta=u[t]-old_u[index]
            comparisons.append(dict(time=t,maxabs=float(np.max(np.abs(delta))),frobenius=float(np.linalg.norm(delta)),relative_frobenius=float(np.linalg.norm(delta)/np.linalg.norm(old_u[index]))))
    sources=[file_binding(p) for p in sorted(HERE.glob('*.py'))]
    bundle=dict(schema='p254-a52-requested-time-reference-v1',times_inverse_ry=times.tolist(),
        coefficient_manifest=coefficient_binding,U=file_binding(upath),hamiltonian_coefficients=file_binding(coeff),
        physical_history_bindings=[history,file_binding(handoff_path),handoff['outputs']['bundle'],old_binding],source_bindings=sources,
        initial_conditions='literal original coherent CAR rho0 and original prepared ionic/field W0 unchanged; literal U0=S0=I',
        settings=SETTINGS,approximation='genuine timeordered original saved moving mean U; original constant bare full Cartesian Gaussian generators, unstable directions retained; no feedback',
        provider_api=dict(module=str(HERE/'general_moving_reference.py'),class_name='SavedMovingReference',constructor='(c,bundle,coefficient_binding,times=None,cache_bytes=128*1024**2)',callbacks=['electron(a,b)','boson(q,a,b)','transitions(q)','S(q,time)','connected(q,t,s,selection)','ordered(q,a,b,selection)','position_indices(q,selection)'],positions='ionic u:0:d; Gamma field a:2d:2d+f; full Gamma phase228; positions114',fractional='produce_u integrates exact requested fractions with existing DOP853; reader only selects exactly present nodes, rejects absent times; never interpolates U',cache='128MiB byte-LRU S plus one q/-q history pair (<50MiB at31nodes); no disk S bank',A41_A50='same electron/boson/transitions/times/c/u callbacks; A50 ordered_coordinates works unchanged'),
        saved_coefficients='same existing AffineSavedHamiltonian exact-degree7 Lobatto construction on saved A42 mean segments; coefficients retained losslessly')
    atomic_json(output/'bundle.json',bundle)
    result=dict(schema='p254-a52-uniform-original-U-producer-v1',status='complete',bundle=file_binding(output/'bundle.json'),
        U_shape=list(u.shape),times=times.tolist(),work=work,hamiltonian_construction_seconds=h.construction_seconds,
        elapsed_seconds=perf_counter()-began,old_four_node_U_comparison=comparisons,
        comparison_policy='raw differences retained: identical DOP853 rtol2e-8/atol2e-10, now restarted each requested integer instead of old boundaries0/1/10/30; altered adaptive steps and accumulated rounding, no forced matching/unitary repair; this is paid new integration, not free recovery',
        U_unitarity_frobenius_max=float(max(np.linalg.norm(x@x.conj().swapaxes(-1,-2)-np.eye(c.m)) for x in u)),
        input_bindings=dict(coefficient_manifest=coefficient_binding,baseline=history,old_transitions=old_binding,old_handoff=file_binding(handoff_path)),
        source_bindings=sources,scope='uniform reference legs only, no uniform selfenergy history or selfconsistent closure/global accuracy claim')
    atomic_json(output/'result.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':main(sys.argv[1])
