"""Actual Gamma/all111/full64q on genuine A45 S and original preparation.

Only original0/1/10/30 pairs; not uniform history or coupled evolution.
Checkpoint contains only the running exact ordered-q sum, never Phi/chunk banks.
"""
import json
from pathlib import Path
import resource
import sys
import time
import numpy as np
from cubic_contour import CubicContour,ordered_matrix,ordered_action
from capture_command import binding,save

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/name) for name in ('0041','0043')]
from material_electron_ion_inputs import load_material,SavedTransitions
from cubic_cartesian_reference import CartesianReference
from oracle_cubic_contour import main as oracle

class ActualReference(SavedTransitions):
    ordered_grid=CartesianReference.ordered_grid


def checkpoint(output,k,arrays,timings):
    path=output/'Gamma-Pi.npz'
    np.savez(path,**arrays,times_inverse_ry=[0.,1.,10.,30.],external_q=[0],external_rows=np.arange(111),completed_internal_q=np.arange(k+1))
    metadata=dict(status='complete' if k==63 else 'partial',completed_internal_q=list(range(k+1)),
        internal_coordinate_count=111,external_rows=list(range(111)),external_q=[0],routed_time_pairs=16*(k+1),
        archive=binding(path),timings=timings)
    tmp=output/'scope.tmp.json';tmp.write_text(json.dumps(metadata,indent=2)+'\n');tmp.replace(output/'scope.json')
    if sum(p.stat().st_size for p in output.iterdir() if p.is_file())>32<<20:raise RuntimeError('hard new payload cap exceeded')


def consume(output):
    scope=json.loads((output/'scope.json').read_bytes())
    if binding(output/'Gamma-Pi.npz')!=scope['archive']:raise ValueError('saved archive differs')
    with np.load(output/'Gamma-Pi.npz',allow_pickle=False) as saved:
        pg=saved['greater'];pl=saved['lesser'];pr=saved['retarded']
        times=saved['times_inverse_ry'];completed=saved['completed_internal_q']
    if pg.shape!=(4,4,111,111) or pl.shape!=pg.shape or pr.shape!=pg.shape:raise ValueError('actual Gamma Pi axes differ')
    theta=(times[:,None]>=times[None,:])[:,:,None,None]
    error=float(np.max(np.abs(pr-theta*(pg-pl))))
    rng=np.random.default_rng(254510);x=rng.normal(size=114)+1j*rng.normal(size=114)
    y=rng.normal(size=114)+1j*rng.normal(size=114)
    embedded=np.zeros((114,114),complex);embedded[:111,:111]=pr[3,1]
    action=embedded@x;dual=embedded.T@y
    noise=.5j*(pg+pl)
    record=dict(status='passed',scope=scope['status'],completed_internal_q=completed.tolist(),archive=scope['archive'],
        times_inverse_ry=times.tolist(),causal_saved_array_error=error,
        original114_embedding_bilinear_error=float(abs(y@action-dual@x)),
        actual_retarded_action_norm=float(np.linalg.norm(action)),actual_noise_norm=float(np.linalg.norm(noise)),
        actual_greater_norm=float(np.linalg.norm(pg)),actual_lesser_norm=float(np.linalg.norm(pl)),
        actual_complex_sample=[float(pg[3,1,0,93].real),float(pg[3,1,0,93].imag)],
        ionic_block_only=True,field_cubic_terms_added=False,coupled_evolution=False)
    if error>2e-12 or record['original114_embedding_bilinear_error']>2e-10 or not all(np.isfinite(a).all() for a in (pg,pl,pr)):raise AssertionError(record)
    save(output/'saved-consumer.json',record);print(json.dumps(record),flush=True)
    return record


def main(output,wall_cap):
    started=time.monotonic();output=Path(output);output.mkdir(exist_ok=False)
    oracle(output)
    c,manifest,coefficient_binding=load_material()
    source_path='/var/tmp/p254-a47-full37-source-001/source.json'
    provider=CubicContour.from_source(source_path)
    bundle=HERE.parent/'0045/actual-original-bundle-001.json'
    reference=ActualReference(c,[0.,1.,10.,30.],bundle,coefficient_binding)
    arrays={name:np.zeros((4,4,111,111),complex) for name in ('greater','lesser','retarded')}
    timings=[];decisions=[];smoke={}
    for k in range(64):
        unit=time.monotonic();l,left,right=provider.pair(0,k)
        wg1,wl1=reference.ordered_grid(k);wg2,wl2=reference.ordered_grid(l)
        vertex_leg_seconds=time.monotonic()-unit
        for a in range(4):
            for b in range(4):
                for name,first,second in (('greater',wg1,wg2),('lesser',wl1,wl2)):
                    arrays[name][a,b]+=ordered_matrix(left,right,-1j*first[a,b],-1j*second[a,b])
        arrays['retarded'][:]=(reference.times[:,None]>=reference.times[None,:])[:,:,None,None]*(arrays['greater']-arrays['lesser'])
        elapsed=time.monotonic()-unit
        timings.append(dict(internal_q=k,seconds=elapsed,vertex_and_actual_legs_seconds=vertex_leg_seconds))
        checkpoint(output,k,arrays,timings)
        print(json.dumps(dict(completed_internal_q=k,work_unit_seconds=elapsed,total_seconds=time.monotonic()-started)),flush=True)
        if k==0:
            x=np.linspace(.1,.9,111)*(1+.37j);y=x[::-1]*(.3-.8j)
            matrix=ordered_matrix(left,right,-1j*wg1[3,1],-1j*wg2[3,1])
            action=ordered_action(left,right,-1j*wg1[3,1],-1j*wg2[3,1],x)
            dual=ordered_action(left,right,-1j*wg1[3,1],-1j*wg2[3,1],y,True)
            smoke=dict(matrixfree_error=float(np.max(np.abs(action-matrix@x))),transpose_error=float(np.max(np.abs(dual-matrix.T@y))),
                matrix_norm=float(np.linalg.norm(matrix)),original_Gamma_field_in_S=reference.transitions(0)[0].shape[-1]==228)
            if max(smoke['matrixfree_error'],smoke['transpose_error'])>2e-10:raise AssertionError(smoke)
            forecast=(time.monotonic()-started)+63*elapsed+15
            decision=dict(first_actual_work_unit=timings[0],forecast_seconds=forecast,wall_cap_seconds=wall_cap,
                fits=forecast<wall_cap-5,action='continue' if forecast<wall_cap-5 else 'preserve actual partial; grant insufficient')
            decisions.append(decision);save(output/'first-unit-decision.json',decision)
            if not decision['fits']:break
        if time.monotonic()-started+max(elapsed*1.5,10)>wall_cap-10:
            decisions.append(dict(action='preserve partial at actual remaining cap',completed_internal_q=k));break
    reference.close();consumer=consume(output)
    record=dict(schema='p254-a51-full-coordinate-cubic-contour-v1',status='complete' if len(timings)==64 else 'partial',
        source=binding(source_path),coefficient_manifest=coefficient_binding,actual_reference_bindings=reference.bindings,
        input_inventory=binding(HERE/'input-inventory-001.json'),oracle=binding(output/'oracle.json'),
        archive=binding(output/'Gamma-Pi.npz'),actual_saved_consumer=binding(output/'saved-consumer.json'),
        completed_internal_q=list(range(len(timings))),external_q=[0],external_coordinates=111,
        original_time_pair_count=16,source_electronic_sectors=128,source_electronic_atoms=37,
        source_proofs='A47 complete full37 common cubic source + A43 passed coefficient/Fock/SymPy; unchanged and carried, not replayed',
        decisions=decisions,actual_material_API_smoke=smoke,elapsed_seconds=time.monotonic()-started,
        maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        payload_bytes=sum(p.stat().st_size for p in output.iterdir() if p.is_file()),
        source_bindings=[binding(p) for p in sorted(HERE.glob('*.py'))],
        conventions=dict(D='-i ordered connected Cartesian covariance',Pi='-i connected cubic force covariance',theta_at_zero=1,
            canonical_N='1/sqrt64 in each vertex; no extra volume/spin factor',units='Pi Ry^2/bohr^2, time inverse Ry'),
        limits=['Only Gamma external material archive exercised; API supports arbitrary original q',
            'Only original0/1/10/30 samples; not uniform-time history','Leading cubic ionic111 block only; no field cubic invention',
            'No coupled evolution, passive dynamics, physical accuracy, damping or lifetime claim'])
    save(output/'result.json',record);print(json.dumps(record),flush=True)

if __name__=='__main__':main(sys.argv[1],float(sys.argv[2]))
