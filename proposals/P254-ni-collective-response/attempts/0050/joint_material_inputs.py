"""Authenticate only original fixed/Gamma banks actually consumed by A50."""
import json
from pathlib import Path
import numpy as np
from joint_gamma_adapter import HERE
from material_electron_ion_inputs import (MANIFEST_SHA, MANIFEST_BYTES, DEFAULT_MANIFEST,
    ReadOnlyQBanks, ReducedCoefficients, ClampedArray, ImmutableItems, file_binding, authenticate,
    SavedTransitions)
from durable_active_baseline import SavedMovingBaseline

BUNDLE = HERE.parent/'0045/actual-original-bundle-001.json'
HISTORY = '/var/tmp/p254-a42-original-moving-baseline-001/baseline-complete.json'


def load_inputs():
    binding = file_binding(DEFAULT_MANIFEST)
    if (binding['sha256'], binding['bytes']) != (MANIFEST_SHA, MANIFEST_BYTES):
        raise ValueError('not original full12/27 manifest')
    manifest = json.loads(Path(DEFAULT_MANIFEST).read_bytes())
    spec = manifest['specification']
    if manifest['status'] != 'complete' or manifest['failures']:
        raise ValueError('incomplete original material')
    reader = ImmutableItems(manifest['items'])
    fixed = [reader.get('slot/%03d'%s, ('h0','rho0','field','field_reverse')) for s in range(128)]
    arrays = reader.get('mechanical')
    reader.get('passive/12')
    authenticate(spec['original_boundary'])
    mechanical = ClampedArray(**arrays)
    c = ReducedCoefficients(mesh=tuple(spec['mesh']),
        **{key:np.asarray(spec['indexing'][key]) for key in
           ('k_indices','q_indices','slot_spin','slot_k','target','minus')},
        **{key:np.asarray([item[key] for item in fixed]) for key in ('h0','rho0','field','field_reverse')},
        mechanical=mechanical, q_provider=ReadOnlyQBanks(reader,spec),
        field_polarizations=arrays['polarizations'])
    if (c.ns,c.n,c.d,c.f,c.m)!=(128,64,111,3,12):
        raise ValueError('full original Gamma scope differs')
    baseline = SavedMovingBaseline(HISTORY,c)
    reference = SavedTransitions(c,[0.,1.,10.,30.],BUNDLE,binding)
    keys = ['mechanical','passive/12','kernel/00'] + ['slot/%03d'%s for s in range(128)] + ['edge/00/%03d'%s for s in range(128)]
    bindings = dict(manifest=binding,bundle=file_binding(BUNDLE),history=file_binding(HISTORY),
                    consumed_items={key:manifest['items'][key] for key in keys},
                    original_boundary=spec['original_boundary'],transitions=reference.bindings,
                    original_preparation_policy='literal W0/rho0 and full original mechanical blocks unchanged')
    return c,baseline,reference,bindings
