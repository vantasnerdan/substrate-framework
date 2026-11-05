"""Qualify changed ionic association on the paid actual q0/k0 full-mode kernel."""
import json
from pathlib import Path
import numpy as np
from cubic_material_source import HERE,atomic_json
from material_electron_ion_inputs import load_material
from cubic_cartesian_reference import CartesianReference
from cubic_periodic_material import PeriodicThird
from cubic_material_kernel import _chunk
from run_material_cubic_kernel import SOURCES
from reduced_projection_storage import ImmutableItems
from substrate_framework.streamed_material_state import file_binding,authenticate

snapshot=HERE/'actual-kernel-001.partial-manifest.original.json'
previous=json.loads(snapshot.read_bytes())
spec=previous['specification']
c,manifest,binding=load_material()
source_path=authenticate(spec['third_source'])
provider=PeriodicThird(source_path,manifest,binding)
reference=CartesianReference(c)
actual,metadata,timings=_chunk(provider,reference,0,0,spec['external_rows'])
expected=ImmutableItems(previous['items']).get('q/00/00')
errors={}
for name in ('mean_force','greater','lesser','retarded'):
    errors[name]=dict(relative=float(np.linalg.norm(actual[name]-expected[name])/max(1.,np.linalg.norm(expected[name]))),
        max_abs=float(np.max(np.abs(actual[name]-expected[name]))))
if any(error['relative']>=2e-12 or error['max_abs']>=2e-12 for error in errors.values()):
    raise AssertionError(errors)
record=dict(schema='p254-a43-ionic-pair-contraction-admission-v1',status='passed',
    scope='Changed ionic association compared to already paid actual source-bound q0/k0 kernel, both full111mode legs and all16original timepairs. Original paid chunk retained, not overwritten. General-q/opposite physical block comparison carried unchanged; only cache capacity16to512 differs from that qualified calculation.',
    errors=errors,tolerance=2e-12,actual_new_chunk_timings=timings,actual_new_chunk_metadata=metadata,
    paid_original_chunk=previous['items']['q/00/00']['archive'],source_binding=spec['third_source'],coefficient_manifest=binding,
    compact_pair_cache=dict(entries=len(provider.ionic._reciprocal_cache),array_bytes=sum(value.nbytes for pair in provider.ionic._reciprocal_cache.values() for value in pair),maximum_entries=512,all64_raw_array_bytes=2*2*37*27*16*64,maximum512_raw_array_bytes=2*2*37*27*16*512,persistent_phi_bank=False),
    general_q_physical_block_proof=file_binding(HERE/'actual-ionic-pair-contraction-admission-001.json'),
    original_partial_cache=file_binding('/var/tmp/p254-a43-material-cubic-kernel-002/q-kernels/manifest.json'),original_partial_cache_snapshot=file_binding(snapshot),
    original_source_bindings=previous['specification']['source_bindings'],current_source_bindings=[file_binding(HERE/name) for name in SOURCES],
    original_source_snapshots=[file_binding(HERE/(name+'.original-actual-kernel-001')) for name in ('cubic_ionic_donor.py','cubic_material_kernel.py','run_material_cubic_kernel.py')]+[file_binding(HERE/'cubic_ionic_donor.py.qualified-ionic-pairs-proof-001')],oracle=file_binding(__file__))
atomic_json(HERE/'actual-ionic-pair-contraction-admission-002.json',record)
print(json.dumps(record),flush=True)
