"""Authenticated fixed-kernel assembly; no interpolated component histories."""
from pathlib import Path
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), *(str(HERE.parent/n) for n in
    ('0052','0050','0049','0045','0042','0041','0030','0029'))]
from material_electron_ion_inputs import load_material
from general_moving_reference import SavedMovingReference
from substrate_framework.streamed_material_state import authenticate, file_binding


def load_step(catalog_path):
    catalog = json.loads(Path(catalog_path).read_bytes())
    if catalog['schema'] != 'p254-actual-fixed-kernel-component-catalog-v1':
        raise ValueError('unsupported actual component catalog')
    times = np.asarray(catalog['selected_times_inverse_ry'], float)
    source_times = np.asarray(catalog['component_times_inverse_ry'], float)
    indices = []
    for t in times:
        found = np.flatnonzero(source_times == t)
        if len(found) != 1:
            raise ValueError('actual component history absent; no interpolation')
        indices.append(int(found[0]))
    if len(times) < 2 or times[0] != 0 or not np.array_equal(
            times, np.arange(len(times))*(times[1]-times[0])):
        raise ValueError('actual uniform complete solver history required')
    files = {key: authenticate(value) for key,value in catalog['archives'].items()}
    authenticate(catalog['reference_bundle'])
    authenticate(catalog['solver_binary'])
    c, manifest, binding = load_material()
    reference = SavedMovingReference(c, catalog['reference_bundle']['path'],
                                     binding, times=times)
    pair = np.ix_(indices, indices)
    sigma, pi, represented = {}, {}, {}
    with np.load(files['ee'], allow_pickle=False) as ee, \
         np.load(files['ei'], allow_pickle=False) as ei, \
         np.load(files['joint'], allow_pickle=False) as joint, \
         np.load(files['represented'], allow_pickle=False) as old, \
         np.load(files['cubic'], allow_pickle=False) as cubic:
        if not np.array_equal(cubic['times_inverse_ry'], source_times):
            raise ValueError('cubic actual sample identity differs')
        if int(cubic['external_q'].item()) != 0 or not np.array_equal(
                cubic['external_rows'], np.arange(c.d)):
            raise ValueError('full original Gamma cubic rows required')
        if not np.array_equal(cubic['completed_internal_q'], np.arange(c.nq)):
            raise ValueError('complete original internal cubic q required')
        for name in ('greater','lesser','retarded'):
            direct_exchange = ee[name+'_total'].transpose(1,2,0,3,4)[pair]
            ionic = ei[name][pair]
            replacement = (joint['joint_sigma_'+name][pair]
                           - joint['ionic_sigma_'+name][pair])
            sigma[name] = direct_exchange + ionic + replacement
            pi[name] = joint['joint_pi_'+name][pair].copy()
            pi[name][:,:,:c.d,:c.d] += cubic[name][pair]
            represented[name] = old['pi_'+name][pair]
    expected = (len(times),len(times),c.ns,c.m,c.m)
    if any(value.shape != expected for value in sigma.values()):
        raise ValueError('complete canonical fermion layout differs')
    if any(value.shape != (len(times),len(times),c.d+c.f,c.d+c.f)
           for value in [*pi.values(),*represented.values()]):
        raise ValueError('full joint Gamma layout differs')
    return c, reference, catalog, sigma, pi, represented
