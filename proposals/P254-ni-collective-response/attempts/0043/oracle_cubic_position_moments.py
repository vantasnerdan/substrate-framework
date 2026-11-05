"""Expose only the changed exact radial/polynomial/FFT block on actual material."""
import json
from pathlib import Path
import numpy as np
from cubic_material_source import restore_original,restore_basis,atomic_json,HERE
from material_electron_ion_inputs import load_material
from reduced_projection_storage import ImmutableItems
from substrate_framework.material_physical_operators import build_physical_factors
from substrate_framework.streamed_material_state import file_binding
from cubic_position_moments import position_factors,compressed_overlaps,position_fft
from cubic_electronic_donor import moment_overlaps,MOMENTS

c,manifest,binding=load_material()
control,inputs,rows,states,original,lower,context=restore_original(manifest)
reader=ImmutableItems(manifest['items'])
slot=rows[0][1]
basis=restore_basis(reader.get('preparation/000',('miller','k_cart_bohr_inv','reciprocal_bohr_inv','row_map','source_miller')),reader.items['preparation/000']['metadata'],slot,control.grid)
old=build_physical_factors(inputs,basis)
new,compact=position_factors(inputs,basis)
rng=np.random.default_rng(4303026)
frame=(rng.normal(size=(len(basis.miller.T),5))+1j*rng.normal(size=(len(basis.miller.T),5)))/np.sqrt(len(basis.miller.T))
p=(basis.reciprocal_bohr_inv@basis.miller+basis.k_cart_bohr_inv[:,None]).T
errors={};mutation=0.
for previous,current in zip(old.species,new.species):
    shape=(len(previous.atom_indices),len(previous.active_d),5)
    errors[previous.name+'_physical_beta']=float(np.linalg.norm(previous.beta-current.beta)/max(1.,np.linalg.norm(previous.beta)))
    expected=moment_overlaps(previous.beta,frame,p,shape)
    actual=compressed_overlaps(current,frame,p,compact[current.name],MOMENTS)
    errors[previous.name+'_all20_moments']=max(float(np.linalg.norm(actual[axes]-expected[axes])/max(1.,np.linalg.norm(expected[axes]))) for axes in MOMENTS)
    radial,upf,positions=compact[current.name]
    wrong=compressed_overlaps(current,frame,p,(radial,upf,-positions),MOMENTS)
    mutation=max(mutation,max(float(np.linalg.norm(wrong[axes]-expected[axes])/max(1.,np.linalg.norm(expected[axes]))) for axes in MOMENTS))
expected=old.fft.from_g(frame)
actual=position_fft(new.fft,frame)
errors['original_full_grid_fft']=float(np.linalg.norm(actual-expected)/max(1.,np.linalg.norm(expected)))
if max(errors.values())>=2e-12 or mutation<=1e-3:
    raise AssertionError(dict(errors=errors,atomic_phase_mutation=mutation))
snapshot=HERE/'actual-prepare-003.partial-manifest.original.json'
previous=json.loads(snapshot.read_bytes())
record=dict(schema='p254-a43-contraction-association-admission-v1',status='passed',
    scope='Only exact angular-polynomial sharing, exact unique squared-modulus radial reuse, in-place single-thread FFT; original actual physical placed-beta/all20 moment/full-grid normalization exposed. No original full-state higher jet replay.',
    errors=errors,tolerance=2e-12,atomic_phase_mutation=mutation,
    carried_association_proof=file_binding(HERE/'actual-contraction-association-admission-001.json'),
    original_partial_cache=file_binding('/var/tmp/p254-a43-full-state-third-002/sector-thirds/manifest.json'),
    original_partial_cache_snapshot=file_binding(snapshot),original_source_bindings=previous['specification']['source_bindings'],
    original_source_snapshots=[file_binding(HERE/(name+'.original-actual-003')) for name in ('cubic_material_source.py','cubic_electronic_donor.py')],
    current_source_bindings=[file_binding(HERE/name) for name in ('cubic_material_source.py','cubic_electronic_donor.py','cubic_position_moments.py')],
    oracle=file_binding(__file__),coefficient_manifest=binding,actual_basis_identity=basis.identity)
atomic_json(HERE/'actual-contraction-association-admission-002.json',record)
print(json.dumps(record),flush=True)
