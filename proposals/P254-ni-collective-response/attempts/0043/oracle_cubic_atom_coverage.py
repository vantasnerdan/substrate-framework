"""Expose exact onsite atom selection, cached slicing and unavailable readouts."""
import json
from types import SimpleNamespace
import numpy as np
from cubic_material_source import HERE,atomic_json
from cubic_electronic_donor import full_column_projector_third,local_third
from cubic_periodic_material import PeriodicThird
from substrate_framework.streamed_material_state import file_binding,authenticate

rng=np.random.default_rng(4304026)
def random(shape):
    return rng.normal(size=shape)+1j*rng.normal(size=shape)
psi=random((7,5)); state=random((5,5)); p=rng.normal(size=(7,3)); beta=np.asfortranarray(random((7,6))); d=random((3,3))
species=SimpleNamespace(atom_indices=(1,2),active_d=d,beta=beta)
base=SimpleNamespace(ionic=SimpleNamespace(inputs=SimpleNamespace(reference=SimpleNamespace(nat=2)),p=p),factors=SimpleNamespace(species=(species,)))
full=full_column_projector_third(psi,state,base)
selected_species=SimpleNamespace(atom_indices=(2,),active_d=d,beta=np.asfortranarray(beta[:,3:]))
selected=SimpleNamespace(ionic=base.ionic,factors=SimpleNamespace(species=(selected_species,)),cubic_atom_indices=(1,))
partial=full_column_projector_third(psi,state,selected)
errors={'nonlocal_selected_original_atom':float(np.linalg.norm(partial-full[[1]])/max(1.,np.linalg.norm(full[[1]])))}
miller=np.indices((2,2,2)).reshape(3,-1); g=miller.T.astype(float)
reference=SimpleNamespace(nat=2,atom_names=('Ni','H'),positions_bohr=rng.normal(size=(3,2)))
ionic=SimpleNamespace(inputs=SimpleNamespace(reference=reference,density=SimpleNamespace(miller=miller)),g=g,vloc={'Ni':random((8,)),'H':random((8,))})
fft=SimpleNamespace(nnr=8,physical=(2,2,2),omega=7.)
coupling=SimpleNamespace(ionic=ionic,factors=SimpleNamespace(fft=fft))
density=random((8,)); full_local=local_third(density,coupling)
coupling.cubic_atom_indices=(1,)
partial_local=local_third(density,coupling)
errors['local_selected_original_atom']=float(np.linalg.norm(partial_local-full_local[[1]]))
full_provider=object.__new__(PeriodicThird); full_provider.n=64; full_provider.electronic_atoms=tuple(range(37))
partial_provider=object.__new__(PeriodicThird); partial_provider.n=64; partial_provider.electronic_atoms=(0,31)
tensor=random((37,3,3,3)); rows=[0,1,2,93,94,95]
expected=full_provider._onsite(tensor,rows)
actual=partial_provider._onsite(tensor[[0,31]],rows)
errors['full111_mode_onsite_vertex_readout']=float(np.linalg.norm(actual-expected))
failed_closed=False
try:
    partial_provider._onsite(tensor[[0,31]],[3])
except ValueError:
    failed_closed=True
if max(errors.values())>=2e-12 or not failed_closed or actual.shape!=(6,111,111):
    raise AssertionError(dict(errors=errors,failed_closed=failed_closed,shape=actual.shape))
snapshot=HERE/'actual-prepare-005.partial-manifest.original.json'
previous=json.loads(snapshot.read_bytes())
paid=[]
for key,item in previous['items'].items():
    with np.load(authenticate(item['archive']),allow_pickle=False) as saved:
        if saved['local'].shape!=(37,3,3,3) or saved['nonlocal_term'].shape!=(37,3,3,3):
            raise ValueError('original paid full-atom archive differs')
        paid.append(dict(key=key,archive=item['archive'],original_shape=list(saved['local'].shape),selected_shape=list(saved['local'][[0,31]].shape)))
record=dict(schema='p254-a43-contraction-association-admission-v1',status='passed',scope='Only requested onsite external atoms0,31; all128 original full-column complex states/weights and all37atom ionic/full111internal modes unchanged. No all37electronic-complete claim. Paid full37archives reused by authenticated slicing, never rewritten.',
    electronic_atom_indices=[0,31],errors=errors,tolerance=2e-12,uncomputed_external_atom_failed_closed=failed_closed,
    saved_paid_sector_reuse=paid,carried_physical_moment_proof=file_binding(HERE/'actual-contraction-association-admission-002.json'),
    rejected_star_transport=file_binding(HERE/'actual-paid-source-star-qualification-001.json'),
    original_partial_cache=file_binding('/var/tmp/p254-a43-full-state-third-002/sector-thirds/manifest.json'),
    original_partial_cache_snapshot=file_binding(snapshot),original_source_bindings=previous['specification']['source_bindings'],
    original_source_snapshots=[file_binding(HERE/(name+'.original-actual-005')) for name in ('cubic_material_source.py','cubic_electronic_donor.py','cubic_position_moments.py')],
    current_source_bindings=[file_binding(HERE/name) for name in ('cubic_material_source.py','cubic_electronic_donor.py','cubic_position_moments.py')],oracle=file_binding(__file__))
atomic_json(HERE/'actual-contraction-association-admission-003.json',record)
print(json.dumps(record),flush=True)
