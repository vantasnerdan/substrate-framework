"""New-sensitive-block proof of full-column moment/D association only.

Existing physical SymPy/Fock/routing proofs are carried, not replayed. The direct
PW-square scalar here is deliberately tiny and independent of production's
projected-moment contraction; no actual material sector is recomputed.
"""
from itertools import permutations
from pathlib import Path
from types import SimpleNamespace
import json
import numpy as np
from cubic_electronic_donor import full_column_projector_third,TRIPLES
from substrate_framework.streamed_material_state import file_binding
from cubic_material_source import atomic_json

HERE=Path(__file__).resolve().parent
rng=np.random.default_rng(4302026)
def complex_array(shape):
    return rng.normal(size=shape)+1j*rng.normal(size=shape)
psi=complex_array((7,5)); state=complex_array((5,5)); p=rng.normal(size=(7,3))
beta=np.asfortranarray(complex_array((7,6))); d=complex_array((3,3))
species=SimpleNamespace(atom_indices=(1,2),active_d=d,beta=beta)
coupling=SimpleNamespace(ionic=SimpleNamespace(inputs=SimpleNamespace(reference=SimpleNamespace(nat=2)),p=p),factors=SimpleNamespace(species=(species,)))
actual=full_column_projector_third(psi,state,coupling)
rho=psi@state@psi.conj().T
expected=np.zeros((2,3,3,3),complex)
for atom in range(2):
    b=beta[:,atom*3:(atom+1)*3]
    for axes in TRIPLES:
        value=0j
        for mask in range(8):
            bra=[axes[leg] for leg in range(3) if mask&(1<<leg)]
            ket=[axes[leg] for leg in range(3) if not mask&(1<<leg)]
            mb=np.prod(p[:,bra],axis=1) if bra else np.ones(len(p))
            mk=np.prod(p[:,ket],axis=1) if ket else np.ones(len(p))
            value+=(-1j)**len(bra)*(1j)**len(ket)*np.einsum('gh,hi,ij,gj,h,g->',rho,b,d,b.conj(),mb,mk,optimize=False)
        for axes_order in set(permutations(axes)):
            expected[(atom,*axes_order)]=value
error=float(np.linalg.norm(actual-expected)/max(1.,np.linalg.norm(expected)))
mutation=float(np.linalg.norm(full_column_projector_third(psi,state.conj(),coupling)-expected)/max(1.,np.linalg.norm(expected)))
if not error<2e-12 or not mutation>1e-3:
    raise AssertionError(dict(relative_error=error,conjugated_state_mutation=mutation))
original_manifest=HERE/'actual-prepare-001.partial-manifest.original.json'
previous=json.loads(original_manifest.read_bytes())
record=dict(schema='p254-a43-contraction-association-admission-v1',status='passed',
    scope='Only associative full-column projector moment/D contraction; physical donor/all-eight legs/state/weights unchanged. Existing physical oracles carried.',
    relative_error=error,tolerance=2e-12,conjugated_state_mutation=mutation,
    original_partial_cache=file_binding('/var/tmp/p254-a43-full-state-third-002/sector-thirds/manifest.json'),
    original_partial_cache_snapshot=file_binding(original_manifest),
    original_source_bindings=previous['specification']['source_bindings'],
    original_source_snapshots=[file_binding(HERE/(name+'.original-actual-001')) for name in ('cubic_material_source.py','cubic_electronic_donor.py')],
    current_source_bindings=[file_binding(HERE/name) for name in ('cubic_material_source.py','cubic_electronic_donor.py')],oracle=file_binding(__file__))
atomic_json(HERE/'actual-contraction-association-admission-001.json',record)
print(json.dumps(record),flush=True)
