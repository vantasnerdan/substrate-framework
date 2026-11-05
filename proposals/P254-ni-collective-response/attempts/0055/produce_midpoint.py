"""Integrate saved A52 Hamiltonian coefficients; preserve original paid U0,U1."""
from refinement_common import *
from moving_reference_propagation import AffineSavedHamiltonian,integrate_u

def run(output):
    output=Path(output);output.mkdir(exist_ok=False)
    old_spec=json.loads(authenticate(catalog()['reference_bundle']).read_bytes())
    coeff_binding=old_spec['hamiltonian_coefficients']
    h=AffineSavedHamiltonian.__new__(AffineSavedHamiltonian)
    with np.load(authenticate(coeff_binding),allow_pickle=False) as z:
        h.times=z['segment_times'];h.coefficients=z['coefficients']
    old_binding=next(x for x in old_spec['physical_history_bindings'] if x['path'].endswith('raw-transitions.npz')) if any(x['path'].endswith('raw-transitions.npz') for x in old_spec['physical_history_bindings']) else None
    if old_binding is None:
        handoff=json.loads((HERE.parent/'0045/actual-complete-original-handoff-001.json').read_bytes())
        old_binding=handoff['outputs']['raw_transitions']
    with np.load(authenticate(old_binding),allow_pickle=False) as z:
        old=z['U'][[0,1]].copy()
    midpoint,work_half=integrate_u(h,0.,.5,old[0])
    raw_one,work_one=integrate_u(h,.5,1.,midpoint)
    u=np.array([old[0],midpoint,old[1]])
    path=output/'U.npy'
    with path.open('xb') as f:np.save(f,u,allow_pickle=False)
    bundle=dict(old_spec,times_inverse_ry=TIMES.tolist(),U=binding(path),
        physical_history_bindings=[*old_spec['physical_history_bindings'],catalog()['reference_bundle'],coeff_binding,old_binding],
        source_bindings=[*old_spec['source_bindings'],binding(__file__),binding(HERE/'refinement_common.py')],
        approximation='literal old U0/U1; genuinely integrated midpoint from saved A52 Hamiltonian coefficients, original lazy Gaussian S; no rho factorization/repair')
    save(output/'bundle.json',bundle)
    receipt=dict(status='complete',times_inverse_ry=TIMES.tolist(),U=binding(path),bundle=binding(output/'bundle.json'),
        coefficients=coeff_binding,old_U_bank=old_binding,work=[work_half,work_one],
        literal_U0_equal=bool(np.array_equal(u[0],old[0])),literal_U1_equal=bool(np.array_equal(u[2],old[1])),
        raw_reintegrated_U1_maxabs=float(np.max(np.abs(raw_one-old[1]))),
        raw_reintegrated_U1_frobenius=float(np.linalg.norm(raw_one-old[1])),
        unitarity_defect=float(np.linalg.norm(midpoint@midpoint.conj().swapaxes(-1,-2)-np.eye(midpoint.shape[-1]))),
        repair_performed=False,scope='genuine U(.5), raw comparison retained; not coupled evolution/global accuracy')
    save(output/'result.json',receipt);return receipt
