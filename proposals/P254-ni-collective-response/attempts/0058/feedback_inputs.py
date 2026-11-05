"""Saved-state/common-functional boundary for an actual material feedback step."""
from pathlib import Path
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), *(str(HERE.parent/n) for n in
    ('0057', '0054', '0053', '0052', '0050', '0047', '0045', '0042', '0041', '0030', '0029'))]
from mean_feedback_reference import MeanFeedbackReference
from material_electron_ion_inputs import load_material
from substrate_framework.streamed_material_state import authenticate, file_binding


def load_feedback(components, mean_result, previous, baseline_catalog):
    """Keep common Pi and represented-only control on distinct bare references."""
    catalog = json.loads(Path(components).read_bytes())
    baseline = json.loads(Path(baseline_catalog).read_bytes())
    mean_spec = json.loads(Path(mean_result).read_bytes())
    if catalog['schema'] != 'p254-dressed-functional-components-v1':
        raise ValueError('complete dressed functional catalog required')
    if catalog['reference_source_kind'] != 'bare-mechanical-with-contacts':
        raise ValueError('this step requires a bare mechanical reference')
    if catalog['pi_convention'] != 'full-common-joint-plus-cubic; represented-only-control-separate':
        raise ValueError('full common Pi required; no represented-channel subtraction')
    authenticate(catalog['input_dataset'])
    if catalog['input_dataset'] != file_binding(previous):
        raise ValueError('kernels and previous dressed state differ')
    authenticate(mean_spec['input_dataset'])
    if mean_spec['input_dataset'] != catalog['input_dataset']:
        raise ValueError('mean and kernels come from different dressed iterates')
    if mean_spec['reference_source_kind'] != catalog['reference_source_kind']:
        raise ValueError('mean and collision reference conventions differ')
    archive = authenticate(catalog['archive'])
    authenticate(baseline['solver_binary'])
    times = np.asarray(catalog['times_inverse_ry'], float)
    if (len(times) < 2 or times[0] != 0 or not np.array_equal(
            times, np.arange(len(times))*(times[1]-times[0]))):
        raise ValueError('complete genuine uniform history required')
    c, manifest, binding = load_material()
    with np.load(previous, allow_pickle=False) as old:
        if not np.array_equal(old['times'], times):
            raise ValueError('previous state time history differs')
        retained = {key: old[key].copy() for key in (
            'electron_ret', 'electron_less', 'gamma_ret', 'gamma_less',
            'source_electron_less0', 'control_ret', 'control_less',
            'source_gamma_ret', 'source_gamma_less')}
    sigma, pi = {}, {}
    with np.load(archive, allow_pickle=False) as data:
        for name in ('greater', 'lesser', 'retarded'):
            sigma[name] = data['sigma_'+name].copy()
            pi[name] = data['pi_'+name].copy()
    expected_e = (len(times), len(times), c.ns, c.m, c.m)
    expected_b = (len(times), len(times), c.d+c.f, c.d+c.f)
    if any(x.shape != expected_e or not np.isfinite(x).all() for x in sigma.values()):
        raise ValueError('complete finite canonical fermion kernels required')
    if any(x.shape != expected_b or not np.isfinite(x).all() for x in pi.values()):
        raise ValueError('complete finite joint Gamma kernels required')
    reference = MeanFeedbackReference(c, mean_result, binding)
    if not np.array_equal(reference.times, times):
        reference.close()
        raise ValueError('mean/reference and collision time histories differ')
    return c, reference, catalog, baseline, retained, sigma, pi
