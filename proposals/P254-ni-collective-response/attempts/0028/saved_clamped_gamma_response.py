"""Same-seed Gamma counterfactual from the existing saved frozen-electron scalar.

Matter coordinates are physical per-cell amplitudes; field coordinates remain
GLOBAL. This finite quadratic evaluation is not a coupled trajectory, covariance,
local-space reconstruction, finite-dose energy, heat or stationary screening.
"""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT / 'src'), str(HERE)]

import numpy as np
from scipy.linalg import expm
from substrate_framework.material_clamped_array import ClampedArray
from substrate_framework.material_control_reference import bound_arrays, bound_json
from substrate_framework.streamed_material_state import file_binding
from mixed_response_storage import MappedMixedState


def array(value):
    return np.asarray(value['real']) + 1j * np.asarray(value['imag'])


def encoded(value):
    value = np.asarray(value)
    if not np.isfinite(value).all():
        raise ValueError('nonfinite canonical control; no repair')
    return dict(real=value.real.tolist(), imag=value.imag.tolist())


def canonical(vector, d, n):
    return dict(u=encoded(vector[:d] / np.sqrt(n)),
                p=encoded(vector[d:2*d] / np.sqrt(n)),
                coordinate=encoded(vector[2*d:2*d+3]),
                momentum=encoded(vector[2*d+3:]))


def build_readout(registration_binding):
    registration = bound_json(registration_binding)
    inputs = registration['source_inputs']
    saved = bound_json(inputs['metadata'])
    identity = saved['identity']
    if (identity['matched_joint_identity'] != registration['model_identity']
            or identity['control'] != 'periodic'
            or identity['mass']['variant'] != 'original'
            or identity['field_order'] != 'a_xyz,Pi_xyz'
            or identity['extra_delta_curvature'] != 'exact zero'):
        raise ValueError('original matched periodic scalar and masses required')
    aggregate = bound_json(identity['aggregate'])
    if aggregate['archive'] != identity['archive']:
        raise ValueError('original scalar archive differs')
    arrays = bound_arrays(identity['archive'])
    prep = identity['preparation']
    masses = array(prep['mass_ry'])
    d = len(masses)
    lattice = arrays['common_force_lattice_ry_bohr2']
    n = int(np.prod(lattice.shape[:3]))
    if n != len(array(identity['cell_indices'])):
        raise ValueError('original cell normalization differs')
    model = ClampedArray(lattice, masses, arrays['valence_charges'],
        np.zeros((d, d)), arrays['contact'],
        np.asarray(arrays['mixed_ionic_field']).reshape(d, 3),
        array(prep['force_per_cell_ry_bohr']), array(prep['field_gradient_per_cell']),
        np.eye(3), prep['field_mass'], array(prep['ionic_ordered']),
        array(prep['field_ordered']), array(prep['trap_omega_ry']))
    initial_mean = bound_json(inputs['mean_initial'])
    endpoint_mean = bound_json(inputs['mean_endpoint'])
    if (initial_mean['model_identity'] != registration['model_identity']
            or endpoint_mean['model_identity'] != registration['model_identity']
            or initial_mean['time_ry_inverse'] != 0):
        raise ValueError('same-model original mean and endpoint time required')
    time = endpoint_mean['time_ry_inverse']
    with MappedMixedState(inputs['mixed_initial']) as initial:
        if initial.identity != registration['model_identity'] or initial.time != 0 or np.any(initial.q):
            raise ValueError('original Gamma mixed seed required')
        expected = np.zeros(d, complex)
        expected[0] = 1 / n
        for component in ('plus', 'minus'):
            small = initial.small[component]
            if (not np.array_equal(small.p, expected) or np.any(small.u)
                    or np.any(small.coordinate) or np.any(small.momentum)):
                raise ValueError('original Ni1/x J=1 cell0 seed required')
        if any(np.any(getattr(initial.small['mixed'], key))
               for key in ('u', 'p', 'coordinate', 'momentum')):
            raise ValueError('original zero mixed initialization required')
        vectors = {key: np.concatenate((np.sqrt(n)*small.u,
            np.sqrt(n)*small.p, small.coordinate, small.momentum))
            for key, small in initial.small.items()}
        seed = initial.record['provenance']['original_local_pulse']
    source_paths = {Path(__file__).resolve()} | {
        Path(module.__file__).resolve() for module in sys.modules.values()
        if getattr(module, '__file__', None)
        and Path(module.__file__).resolve().is_relative_to(ROOT)}
    source_bindings = [file_binding(path) for path in sorted(source_paths)]
    generator = model.generators[0]
    transition = expm(time * generator)
    augmented = np.zeros((len(generator)+1, len(generator)+1), complex)
    augmented[:-1, :-1] = generator
    augmented[:-1, -1] = model.j[0] @ model.linear
    drift = expm(time * augmented)[:-1, -1]
    mean = np.concatenate((np.sqrt(n)*np.asarray(initial_mean['u_bohr']),
        np.sqrt(n)*np.asarray(initial_mean['p_canonical']),
        initial_mean['field_coordinate'], initial_mean['field_momentum']))
    if source_bindings != [file_binding(path) for path in sorted(source_paths)]:
        raise ValueError('consumer source changed during finite control evaluation')
    return dict(schema='P254-saved-clamped-Gamma-response-v1', status='complete',
        model_identity=registration['model_identity'], time_ry_inverse=time,
        q=[0, 0, 0], ncells=n, original_seed=seed, source_inputs=inputs,
        aggregate=identity['aggregate'], archive=identity['archive'],
        initial={key: canonical(vector, d, n) for key, vector in vectors.items()},
        endpoint={key: canonical(transition @ vector, d, n) for key, vector in vectors.items()},
        own_unexcited_baseline=canonical(transition @ mean + drift, d, n),
        source_bindings=source_bindings,
        software=dict(numpy=np.__version__, scipy=sys.modules['scipy'].__version__),
        interpretation='Existing frozen-electron quadratic scalar at the coupled time and original canonical seed; physical per-cell matter, one global Cartesian field. Own baseline separate; not a stationary screening correction, full spatial response, coupled trajectory, finite-dose energy, covariance or heat.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration', type=Path, required=True)
    parser.add_argument('--registration-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    binding = file_binding(args.registration)
    if binding['sha256'] != args.registration_sha256:
        raise ValueError('registration digest differs')
    report = build_readout(binding)
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status=report['status'], output=file_binding(args.output))))
