"""Exact raw-force gate, independent norms, print enclosure and COM false green."""
import copy
import numpy as np
import sympy as sp
from scipy.linalg import norm
from substrate_framework.material_primary_boundary import squared_norms_decimal, raw_gate, force_rounding_upper, com_project


def run(ledger,c):
    out=c['native']; p=c['published']['data']; s=c['scf']['data']; factor=c['units']['force_factor']
    native=out['forces']*factor; published=np.array(p['forces_ev_per_angstrom']); gate=c['settings']['optimizer']['fmax_ev_per_angstrom']
    eps=64*np.finfo(float).eps
    allowance=eps*max(1,float(np.max(abs(native))))
    predicate=lambda f:np.asarray(f).shape==(37,3) and np.isfinite(f).all() and np.max(abs(np.asarray(f)-native))<=allowance
    changed=published.copy(); changed[0,1]+=0.001
    sign=published.copy(); sign[0,1]*=-1
    ledger.mutation_sensitive('independent native total force conversion publication',predicate,published,[changed,sign,published[::-1]])
    ledger.check('SCF and relaxation force arrays agree exactly',s['forces_ev_per_angstrom']==p['forces_ev_per_angstrom'])
    ledger.check('SCF snapshot force serialization matches eight decimals',np.max(abs(c['snapshot']['forces']-native))<=0.5e-8+allowance)
    squared=squared_norms_decimal(published)
    exact_norm=np.array([float(sp.sqrt(x)) for x in squared])
    scipy_norm=np.array([norm(row) for row in published])
    ledger.check('all 37 exact decimal norms agree with independent SciPy BLAS norms',np.max(abs(exact_norm-scipy_norm))<=allowance)
    ledger.mutation_sensitive('raw stationarity frozen 0.005 gate',lambda f:raw_gate(f,gate),published,
        [published+np.eye(37,3)*0.006, np.full((37,3),0.006)])
    fmax=float(max(exact_norm))
    ledger.mutation_sensitive('published fmax matches independent maximum',lambda x:abs(x-fmax)<=allowance,p['raw_fmax_ev_per_angstrom'],[p['raw_fmax_ev_per_angstrom']+1e-4])
    upper=force_rounding_upper(native,0.5e-8*factor)
    ledger.check('raw force gate remains true throughout stdout rounding enclosure',max(upper)<=gate)
    ledger.check('source-derived print enclosure exposes threshold-crossing input',max(force_rounding_upper(native*2,0.5e-8*factor))>gate)
    masses=np.array(c['final']['masses_amu']); projected=com_project(published,masses)
    from ase import Atoms
    from ase.constraints import FixCom
    atoms=Atoms(c['final']['symbols'],masses=masses)
    ase_project=published.copy(); FixCom().adjust_forces(atoms,ase_project)
    ledger.check('independent COM projection agrees with source implementation',np.max(abs(projected-ase_project))<=allowance)
    # A force in the removed direction is a substantive optimizer false-green.
    false_green=masses[:,None]*np.array([[0.001,0,0]])
    ledger.check('COM gauge cannot substitute for raw stationarity',raw_gate(com_project(false_green,masses),gate) and not raw_gate(false_green,gate))
    state=c['optimizer']['state']
    ledger.mutation_sensitive('published optimizer convergence is necessary',lambda v:v is True,p['optimizer_converged'],[False])
    ledger.check('optimizer accepted step count equals durable checkpoint',p['steps_completed']==state['steps_completed'] and 0<=state['steps_completed']<=c['settings']['optimizer']['max_steps'])
    ledger.check('raw acceptance with only translation gauge',p['translation_gauge']=='centre_of_mass' and raw_gate(projected,gate) and raw_gate(published,gate))
    energy=out['energy_ry']*c['units']['rydberg_ev']; energy_allowance=64*np.finfo(float).eps*abs(energy)
    ledger.mutation_sensitive('native free energy publication conversion',lambda e:abs(e-energy)<=energy_allowance,p['energy_ev'],[p['energy_ev']+1.])
    ledger.check('SCF snapshot and publication energy identity',s['energy_ev']==p['energy_ev']==c['snapshot']['energy']==c['snapshot']['free_energy'])
    ledger.mutation_sensitive('native convergence completion and version',lambda x:x['version']=='7.5' and x['iterations']==38 and x['done_count']==1,out,[{**out,'done_count':0},{**out,'version':'6.0'}])
    ledger.check('SCF published solver status exact',s['scf_converged'] and s['job_done'] and s['returncode']==0 and s['qe_version']==out['version'])
    return {'per_atom':[{'index_1based':i+1,'species':c['final']['symbols'][i],'force_ev_per_angstrom':row.tolist(),
               'squared_norm_exact_decimal':str(squared[i]),'norm_ev_per_angstrom':float(exact_norm[i]),
               'stdout_rounding_upper_ev_per_angstrom':float(upper[i])} for i,row in enumerate(published)],
            'raw_fmax_ev_per_angstrom':fmax,'gate_ev_per_angstrom':gate,'margin_ev_per_angstrom':gate-fmax,
            'max_rounding_upper_ev_per_angstrom':float(max(upper)),
            'force_component_half_quantum_ry_per_bohr':0.5e-8,'arithmetic_allowance_ev_per_angstrom':allowance,
            'raw_force_sum_ev_per_angstrom':np.sum(published,axis=0).tolist(),
            'com_projected_fmax_ev_per_angstrom':float(max(norm(row) for row in projected)),
            'com_projection_change_max_ev_per_angstrom':float(np.max(abs(projected-published))),
            'energy_ry':out['energy_ry'],'energy_ev':energy,'energy_print_half_quantum_ry':0.5e-8,
            'native_scf_iterations':out['iterations'],'optimizer_steps':state['steps_completed'],
            'scope':'Gate on the published primary-resolution raw forces; source rounding only, not DFT discretization uncertainty.'}
