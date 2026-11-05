"""A45 original state/history authority; no new preparation or mean integration."""
import json
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src')] + [str(HERE.parent/n) for n in ('0045','0042','0041','0030','0029')]
import numpy as np
from substrate_framework.streamed_material_state import authenticate, file_binding
from material_electron_ion_inputs import load_material, SavedTransitions
from durable_active_baseline import SavedMovingBaseline

CONTRACT = dict(path=str(HERE/'implementation-contract-001.json'), bytes=11007,
    sha256='849c17684f6c04cf2ebf23a4a67ccce24877c82cfe79b832ae82aba261924b85')
TIMES = [0.,1.,10.,30.]

def inputs():
    contract = json.loads(authenticate(CONTRACT).read_bytes())
    for source in contract['settled_sources']:
        authenticate(source)
    c, manifest, binding = load_material(authenticate(contract['settled_material_manifest']))
    history = contract['settled_original_moving_baseline']
    baseline = SavedMovingBaseline(authenticate(history), c)
    if baseline.receipt['settings'] != dict(method='DOP853',rtol=2e-8,atol=2e-10):
        raise ValueError('original integration tolerances changed')
    for source in baseline.receipt['implementation']:
        authenticate(source)
    source = baseline.receipt['source']
    for key in ('original_result','baseline0','original_input_manifest',
                'original_coefficient_manifest_snapshot','borrowed_coefficient_manifest'):
        authenticate(source[key])
    if source['borrowed_coefficient_manifest'] != binding:
        raise ValueError('dense mean and original Gaussian use different coefficient authority')
    return c, baseline, binding, history

def own_sources():
    return [file_binding(HERE/name) for name in ('implementation-contract-001.json',
        'moving_reference_inputs.py','moving_reference_propagation.py',
        'oracle_moving_reference.py','run_moving_reference.py','native_wave.py')]

def atomic_json(path, value):
    path=Path(path); temporary=path.with_suffix(path.suffix+'.tmp')
    with temporary.open('w') as stream:
        json.dump(value,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')
        stream.flush()
        import os
        os.fsync(stream.fileno())
    temporary.replace(path)
