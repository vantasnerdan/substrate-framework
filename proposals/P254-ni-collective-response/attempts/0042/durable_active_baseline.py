"""Durably save/read one explicitly NEW original-settings baseline generation.

Only public solve_ivp result.t/y/sol are used. Each accepted DOP853 state
interpolant is degree7; h is derived by the EXISTING baseline_hamiltonian.
No U, old dense history recovery, finite-dose campaign or free resume is claimed.
"""
from pathlib import Path
import argparse
import json
import sys
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src')] + [str(HERE.parent/n) for n in ('0041','0030','0029','0028')]
import numpy as np
from canonical_force_filter import _TRANSFORM, _NODES
from saved_active_baseline_inputs import _execution_bindings
from reduced_electronic_dynamics import ActiveBandDynamics, BaselineState
from reduced_projection_cache import atomic_json
from material_electron_ion_inputs import load_material
from substrate_framework.streamed_material_state import authenticate, file_binding
from time import perf_counter

CONTRACT=dict(path=str(HERE/'implementation-contract-001.json'),bytes=7438,
    sha256='6cb22f829d8fd3450d0dba2606a57f50695b5382d213d4d06757182312695755')


class SavedMovingBaseline:
    def __init__(self, receipt_path, coefficients):
        self.binding = file_binding(receipt_path)
        self.receipt = json.loads(authenticate(self.binding).read_bytes())
        self.metadata = self.receipt
        if self.receipt['schema'] != 'P254-A42-durable-moving-active-baseline-v1' or self.receipt['status'] != 'complete':
            raise ValueError('complete durable A42 reconstructed baseline required')
        self.model = ActiveBandDynamics(coefficients)
        self.gamma = coefficients.q(0)
        self.arrays = {key: np.load(authenticate(binding), mmap_mode='r', allow_pickle=False)
                       for key,binding in self.receipt['arrays'].items()}
        self.times = self.arrays['boundary_times']
        width = coefficients.ns*coefficients.m**2+coefficients.mechanical_sizes[0]
        count = len(self.times)-1
        if (self.times[0] != 0 or np.any(np.diff(self.times)<=0)
            or self.arrays['chebyshev'].shape != (count,8,width)
            or self.arrays['boundary_states'].shape != (count+1,width)):
            raise ValueError('durable baseline shape/intervals differ')
        specification = self.receipt['coefficient_specification']
        if (list(coefficients.mesh)!=specification['mesh'] or coefficients.m!=specification['bands']
            or coefficients.ns!=len(specification['indexing']['slot_spin'])):
            raise ValueError('baseline numeric coefficient layout differs')
        initial = self.model.unpack_baseline(self.arrays['boundary_states'][0])
        if not np.array_equal(initial.rho,coefficients.rho0) or np.any(initial.x):
            raise ValueError('original unexcited initial baseline changed')

    def state(self, time):
        time=float(time)
        if not np.isfinite(time) or not self.times[0]<=time<=self.times[-1]:
            raise ValueError('time outside saved reconstructed moving baseline')
        right=int(np.searchsorted(self.times,time,side='left'))
        if right<len(self.times) and self.times[right]==time:
            value=self.arrays['boundary_states'][right]
        else:
            left,end=self.times[right-1:right+1]
            x=2*(time-left)/(end-left)-1
            value=np.polynomial.chebyshev.chebval(x,self.arrays['chebyshev'][right-1])
        return self.model.unpack_baseline(value)

    def __call__(self,time):
        state=self.state(time)
        return state.rho,self.model.baseline_hamiltonian(state,self.gamma)


def produce(original_result, output, *, max_bytes=256*1024**2):
    """One authorized reconstruction, original initial state/H/settings, no study.

    Authenticate original execution/source history; borrow complete union banks.
    Existing ActiveBandDynamics owns the sole original-settings baseline solve.
    No derivative Gamma view or copied old bank is constructed.
    A pre-write exact payload guard stops rather than truncating intervals/state.
    No existing output directory is accepted and no resume/retry is automatic.
    """
    directory=Path(output).resolve();directory.mkdir(parents=True,exist_ok=False)
    started=perf_counter()
    original_binding=file_binding(original_result)
    original=json.loads(authenticate(original_binding).read_bytes())
    config=original['configuration']
    if (original['status']!='complete' or config['bands']!=12 or config['features']!=27
        or config['times']!=[0.,1.,10.,30.] or config['baseline_only']):
        raise ValueError('original complete loaded12/27 finite campaign result required')
    specification=config['coefficient_specification']
    snapshot=json.loads(authenticate(original['original_coefficient_manifest_snapshot']).read_bytes())
    authority,execution=_execution_bindings(original,specification,snapshot)
    contract=json.loads(authenticate(CONTRACT).read_bytes())
    manifest_binding=contract['settled_material']['manifest']
    c,manifest,manifest_binding=load_material(manifest_binding['path'])
    keys=(['slot/%03d'%s for s in range(128)]+['mechanical']+
          ['edge/00/%03d'%s for s in range(128)]+['kernel/00'])
    for key in keys:
        if manifest['items'][key]['archive']!=snapshot['items'][key]['archive']:
            raise ValueError('settled union differs from original baseline input: '+key)
    initial_binding=original['completed'][0]['baseline_state']
    with np.load(authenticate(initial_binding),allow_pickle=False) as archive:
        initial=BaselineState(archive['rho'],archive['x'])
    if not np.array_equal(initial.rho,c.rho0) or np.any(initial.x):
        raise ValueError('literal original unexcited baseline0 changed')
    model=ActiveBandDynamics(c)
    baseline=model.integrate_baseline((0.,30.),initial,rtol=config['rtol'],atol=config['atol'],method='DOP853')
    result=baseline.result
    source=dict(original_result=original_binding,baseline0=initial_binding,
        original_input_manifest=original['original_input_manifest'],
        original_coefficient_manifest_snapshot=original['original_coefficient_manifest_snapshot'],
        borrowed_coefficient_manifest=manifest_binding,original_implementation_authorities=authority,
        execution_implementation=execution,original_boundary=specification['original_boundary'])
    work=dict(new_baseline_integrations=1,finite_dose_replays=0,preparation_replays=0,
        oldbank_copies=0,baseline_nfev=int(result.nfev),solve_seconds=perf_counter()-started)
    times=np.asarray(result.t);width=result.y.shape[0];count=len(times)-1
    payload=((count*8+count+1)*width*16+len(times)*8)
    if payload>max_bytes:
        raise ValueError(f'complete baseline payload {payload} exceeds authorized {max_bytes}; no partial baseline published')
    bindings={}
    for key,value in (('boundary_times',times),('boundary_states',result.y.T)):
        path=directory/(key+'.npy')
        with path.open('xb') as stream:np.save(stream,value,allow_pickle=False)
        bindings[key]=file_binding(path)
    path=directory/'chebyshev.npy'
    coefficients=np.lib.format.open_memmap(path,mode='w+',dtype=np.complex128,shape=(count,8,width))
    for segment,(left,right) in enumerate(zip(times[:-1],times[1:])):
        nodes=np.array([right if j==0 else left if j==7 else
                        left+(right-left)*(node+1)/2 for j,node in enumerate(_NODES)])
        coefficients[segment]=_TRANSFORM@result.sol(nodes).T
    coefficients.flush();del coefficients
    bindings['chebyshev']=file_binding(path)
    # Original receipt and exact integer metadata remain raw Python JSON values.
    receipt=dict(schema='P254-A42-durable-moving-active-baseline-v1',status='complete',
        arrays=bindings,source=source,work=work,
        coefficient_specification=specification,
        settings=dict(method='DOP853',rtol=config['rtol'],atol=config['atol']),
        generation='NEW one-generation necessary original-settings reconstruction; NOT original execution/cache/free resume',
        state='literal original unexcited rho/x, unchanged active Hamiltonian/settings',
        interpolation='public degree7 DOP853 state interpolant re-expressed by8 Lobatto samples; exact accepted boundaries',
        h_policy='existing ActiveBandDynamics.baseline_hamiltonian of restored rho/x; no duplicated H/G bank',
        transition_propagators_saved=False,payload_bytes=int(payload),
        implementation=[file_binding(path) for path in (Path(__file__),
            HERE.parent/'0041/material_electron_ion_inputs.py',
            HERE.parent/'0030/canonical_force_filter.py',
            HERE.parent/'0030/saved_active_baseline_inputs.py')])
    atomic_json(directory/'baseline-complete.json',receipt)
    return file_binding(directory/'baseline-complete.json')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original-result',required=True)
    parser.add_argument('--output',required=True)
    parser.add_argument('--max-bytes',type=int,default=256*1024**2)
    args=parser.parse_args()
    print(json.dumps(produce(args.original_result,args.output,max_bytes=args.max_bytes),allow_nan=False))
