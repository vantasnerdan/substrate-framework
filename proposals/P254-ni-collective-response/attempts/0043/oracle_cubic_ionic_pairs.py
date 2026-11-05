"""Compare only new ionic pair association/pruning to archived qualified donor."""
import importlib.machinery
import importlib.util
import json
from pathlib import Path
import numpy as np
from cubic_ionic_donor import IonicThird
from cubic_material_source import HERE,atomic_json
from run_material_cubic_kernel import SOURCES
from substrate_framework.streamed_material_state import file_binding

path=HERE/'cubic_ionic_donor.py.original-actual-kernel-001'
loader=importlib.machinery.SourceFileLoader('original_qualified_ionic_donor',str(path))
spec=importlib.util.spec_from_loader(loader.name,loader)
module=importlib.util.module_from_spec(spec); loader.exec_module(module)
cell=np.diag([4.1,4.7,5.3]); positions=np.array([[.21,.38,.67],[1.17,1.56,2.41],[2.32,.83,3.76]])
charges=np.array([2.,4.,7.]); mesh=(4,4,4); rows=[0,1,2,6,7,8]
old=module.IonicThird(cell,positions,charges,mesh)
new=IonicThird(cell,positions,charges,mesh)
triples=[np.zeros((3,3),int),np.array([[0,0,0],[1,2,3],[3,2,1]]),np.array([[1,0,2],[2,3,1],[1,1,1]])]
errors={}
for index,q in enumerate(triples):
    for piece in ('real','reciprocal'):
        expected=old.block(q,rows,piece=piece)
        actual=new.block(q,rows,piece=piece)
        errors[str(index)+'_'+piece]=float(np.linalg.norm(actual-expected)/max(1.,np.linalg.norm(expected)))
    expected=old.block(-q,rows,piece='total')
    actual=new.block(-q,rows,piece='total')
    errors[str(index)+'_independent_opposite_total']=float(np.linalg.norm(actual-expected)/max(1.,np.linalg.norm(expected)))
forward,reverse=new._reciprocal_pairs(np.array([1,2,3]),(0,2))
same=new._reciprocal_pairs(np.array([1,2,3]),(0,2))
cache_hit_reuses_actual_arrays=same[0] is forward and same[1] is reverse
orientation_mutation=float(np.linalg.norm(forward[0,1]-reverse[1,0]))
if max(errors.values())>=2e-12 or not cache_hit_reuses_actual_arrays or orientation_mutation<=1e-3 or len(new._reciprocal_cache)>16:
    raise AssertionError(dict(errors=errors,cache_hit=cache_hit_reuses_actual_arrays,orientation_mutation=orientation_mutation))
snapshot=HERE/'actual-kernel-001.partial-manifest.original.json'
previous=json.loads(snapshot.read_bytes())
record=dict(schema='p254-a43-ionic-pair-contraction-admission-v1',status='passed',scope='Exact original physical reciprocal spheres and all8subset signs; independent pair orientations; external-row early pruning preserves both complete internal mode legs. Bounded16exact-shift LRU, not an allq/Phi bank. Existing physical ionic SymPy/finite-BvK proof carried.',
    errors=errors,tolerance=2e-12,cache_hit_reuses_actual_arrays=cache_hit_reuses_actual_arrays,pair_orientation_mutation=orientation_mutation,
    original_partial_cache=file_binding('/var/tmp/p254-a43-material-cubic-kernel-002/q-kernels/manifest.json'),original_partial_cache_snapshot=file_binding(snapshot),
    original_source_bindings=previous['specification']['source_bindings'],current_source_bindings=[file_binding(HERE/name) for name in SOURCES],
    original_source_snapshots=[file_binding(HERE/(name+'.original-actual-kernel-001')) for name in ('cubic_ionic_donor.py','run_material_cubic_kernel.py')],oracle=file_binding(__file__))
atomic_json(HERE/'actual-ionic-pair-contraction-admission-001.json',record)
print(json.dumps(record),flush=True)
