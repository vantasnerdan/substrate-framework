"""Five genuine midpoint pairs; four literal paid pairs, never reevaluated.

Independent scheduler entry points are run_ei, run_joint and run_represented.
The electron consumer is EI_raw + joint_Sigma - ionic_Gamma_Sigma, exactly
once. The raw archives deliberately remain separate: represented subtraction
is Pi-only and is never applied to electron Sigma. Scientific admission and
cutover limitations are inherited unchanged from A41/A48/A50/A52/A53.
"""
from pathlib import Path
from dataclasses import replace

from refinement_common import (HERE, TIMES, NEW_PAIRS, OLD_MAP, binding, save,
                               catalog, authenticate, reference)
import numpy as np
from joint_material_inputs import load_inputs
from joint_gamma_adapter import common_pair, vertices, ordered_coordinates
from electron_ion_functional import EdgePair, evaluate_components
from material_electron_ion_inputs import edge_pair
from represented_electron_channel import RepresentedElectronChannel, CUTOVER
from substrate_framework.partitioned_q_storage import release_mapped_pages

CHANNELS = ('greater', 'lesser', 'retarded')


def _output(output):
    output = Path(output).resolve()
    if HERE / 'material' not in output.parents and output != HERE / 'material':
        raise ValueError('owned numeric payload must be under 0055/material')
    output.mkdir(parents=True, exist_ok=True)
    return output


def _inputs():
    c, baseline, old_reference, bindings = load_inputs()
    old_reference.close()
    ref = reference(c, bindings['manifest'])
    if not np.array_equal(ref.times, TIMES):
        raise ValueError('exact three-node saved reference required')
    if (c.ns, c.n, c.m, c.d, c.f, len(c.q_indices)) != (128, 64, 12, 111, 3, 64):
        raise ValueError('full original CAR/q/coordinate scope required')
    return c, baseline, ref, bindings


def _old(kind):
    entry = catalog()['archives'][kind]
    return np.load(authenticate(entry), allow_pickle=False), entry


def _merge(old, new):
    """Time-first layout, retaining the original four entries bit for bit."""
    result = np.empty((3, 3, *old.shape[2:]), dtype=old.dtype)
    for a, oa in OLD_MAP.items():
        for b, ob in OLD_MAP.items():
            result[a, b] = old[oa, ob]
    for k, (a, b) in enumerate(NEW_PAIRS):
        result[a, b] = new[k]
    return result


def _finite(value):
    value = np.asarray(value)
    if not np.isfinite(value).all():
        raise FloatingPointError('nonfinite actual refinement; no clipping or repair')
    return value


def _scalar(value):
    value = complex(value)
    return [float(value.real), float(value.imag)]


def _consumer(path):
    """Read persisted values and actually apply every kernel to a fixed vector."""
    result = {}
    with np.load(path, allow_pickle=False) as archive:
        for key in archive.files:
            value = _finite(archive[key])
            vector = np.linspace(0.5, 1.5, value.shape[-1])
            action = _finite(value @ vector)
            result[key] = dict(shape=list(value.shape),
                norm=float(np.linalg.norm(value)), action_norm=float(np.linalg.norm(action)),
                action_sum=_scalar(action.sum()))
    return result


def _finish(output, kind, archive, c, ref, bindings, old_binding, sensitivity, **scope):
    consumer = _consumer(archive)
    record = dict(schema='p254-a55-five-new-pair-' + kind + '-v1', status='complete',
        archive=binding(archive), times_inverse_ry=TIMES.tolist(),
        evaluated_pairs=[list(pair) for pair in NEW_PAIRS],
        literal_old_mapping={str(k): v for k, v in OLD_MAP.items()},
        literal_old_pair_count=4, original_archive=old_binding,
        coefficient_manifest=bindings['manifest'], saved_reference=ref.bindings,
        source=binding(__file__), explicit_slots=list(range(c.ns)), bands=c.m,
        sensitivity=sensitivity, consumer=consumer,
        inherited_scientific_scope='A41/A48/A50/A52/A53 unchanged; fixed kernels only; no selfconsistent closure, physical rates, new proof, or XC double-count prescription',
        **scope)
    save(output / 'consumer.json', record)
    return record


def run_ei(output):
    """All64 internal q, all128 target Sigma; only external Gamma ionic Pi.

    A48 sigma-0.npz contains unprefixed greater/lesser/retarded arrays with
    [time,time,slot,band,band] axes, not a slot-first e-e archive.
    """
    output = _output(output)
    c, baseline, ref, bindings = _inputs()
    sigma = {name: np.zeros((5, c.ns, c.m, c.m), complex) for name in CHANNELS}
    pi = {name: np.zeros((5, c.d, c.d), complex) for name in CHANNELS}
    electrons = {pair: ref.electron(*pair) for pair in NEW_PAIRS}
    sensitivity = None
    count = 0
    try:
        for q in range(len(c.q_indices)):
            data = c.q(q)
            try:
                for k, (a, b) in enumerate(NEW_PAIRS):
                    boson = ref.boson(q, a, b)
                    for source in range(c.ns):
                        target = int(c.target[q, source])
                        pair = edge_pair(data, source, target, electrons[a, b],
                                         electrons[b, a], boson)
                        values = evaluate_components(pair, n=c.n,
                            t=float(TIMES[a]), s=float(TIMES[b]),
                            components=('sigma', 'pi') if q == 0 else ('sigma',))
                        for name in CHANNELS:
                            sigma[name][k, target] += values['sigma_' + name]
                            if q == 0:
                                pi[name][k] += values['pi_' + name]
                        count += 1
                        if q == 0 and k == 0 and source == 0:
                            # One affected orientation comparison only, not a campaign.
                            changed = evaluate_components(
                                replace(pair, right_s=pair.right_s.swapaxes(-1, -2)),
                                n=c.n, t=float(TIMES[a]), s=float(TIMES[b]),
                                components=('sigma',))
                            delta = _finite(changed['sigma_greater'] - values['sigma_greater'])
                            sensitivity = dict(pair=[a, b], q=0, source=0,
                                comparison='stored reverse band orientation versus transposed reverse at right_s only',
                                reference_norm=float(np.linalg.norm(values['sigma_greater'])),
                                affected_difference_norm=float(np.linalg.norm(delta)),
                                affected_action=_scalar(delta.sum()))
            finally:
                for name in ('ionic', 'ionic_reverse', 'density', 'potential', 'kernel'):
                    release_mapped_pages(getattr(data, name))
        old, old_binding = _old('ei')
        # Gamma ionic Pi comes literally from the catalog's common archive;
        # no old Pi-q banks are copied or generated.
        old_joint, old_joint_binding = _old('joint')
        try:
            arrays = {name: _merge(old[name], _finite(sigma[name])) for name in CHANNELS}
            arrays.update({'pi_' + name: _merge(old_joint['ionic_pi_' + name], _finite(pi[name]))
                           for name in CHANNELS})
        finally:
            old.close()
            old_joint.close()
        path = output / 'ei.npz'
        np.savez(path, **arrays)
        return _finish(output, 'ei', path, c, ref, bindings, old_binding, sensitivity,
            internal_q=list(range(64)), pi_external_q=[0], ionic_coordinates=c.d,
            edge_pair_evaluations=count, literal_pi_archive=old_joint_binding,
            normalization='unscaled independently stored forward/reverse vertices; exactly 1/N=1/64; no spin factor',
            consumer_combination='Sigma_EI_final = ei[channel] + joint[joint_sigma_channel] - joint[ionic_sigma_channel], exactly once; pi_* here is Gamma ionic diagnostic, not an additional insertion')
    finally:
        ref.close()


def _joint_sensitivity(c, ref):
    """One original edge with only the canonical field vertex mask changed."""
    a, b = NEW_PAIRS[0]
    forward, reverse = vertices(c)
    dg, dl = ordered_coordinates(ref, a, b)
    gg, gl = ref.electron(a, b)
    bg, bl = ref.electron(b, a)
    source = 0
    target = int(c.target[0, source])
    left = forward[source]
    right = reverse[source].transpose(0, 2, 1)
    pair = EdgePair(left, right, left, right, gg[source], gl[source],
                    gg[target], gl[target], bg[source], bl[source], dg, dl)
    full = evaluate_components(pair, n=1, t=float(TIMES[a]), s=float(TIMES[b]))
    masked_left, masked_right = left.copy(), right.copy()
    masked_left[c.d:] = 0
    masked_right[c.d:] = 0
    masked = evaluate_components(replace(pair, left_t=masked_left, left_s=masked_left,
        right_t=masked_right, right_s=masked_right), n=1,
        t=float(TIMES[a]), s=float(TIMES[b]))
    return dict(pair=[a, b], source=source,
        comparison='one-edge direct canonical-field vertex mask; original ionic vertices retained',
        affected_difference_norms={key: float(np.linalg.norm(_finite(full[key] - masked[key])))
                                   for key in full},
        affected_sigma_action=_scalar((full['sigma_greater'] - masked['sigma_greater']).sum()))


def run_joint(output):
    """A50 all114 joint and ionic Gamma components, new five pairs only."""
    output = _output(output)
    c, baseline, ref, bindings = _inputs()
    new = {}
    try:
        for k, (a, b) in enumerate(NEW_PAIRS):
            for prefix, values in (('joint_', common_pair(c, ref, a, b)),
                                   ('ionic_', common_pair(c, ref, a, b, ionic_only=True))):
                for name, value in values.items():
                    key = prefix + name
                    if key not in new:
                        new[key] = np.empty((5, *value.shape), complex)
                    new[key][k] = _finite(value)
        sensitivity = _joint_sensitivity(c, ref)
        old, old_binding = _old('joint')
        try:
            arrays = {key: _merge(old[key], value) for key, value in new.items()}
        finally:
            old.close()
        path = output / 'common.npz'
        np.savez(path, **arrays)
        return _finish(output, 'joint', path, c, ref, bindings, old_binding, sensitivity,
            external_q=[0], joint_coordinates=c.d+c.f, ionic_coordinates=c.d,
            common_pair_evaluations=10,
            normalization='A50 normalized ionic/sqrtN and canonical field vertices, n=1; both stored gradient orientations, no second 1/N',
            consumer_combination='Use joint_pi_* once for Gamma. Final electron EI Sigma = raw all64 EI + joint_sigma_* - ionic_sigma_* exactly once. Do not add both joint and ionic Gamma Sigma; represented subtraction is Pi-only.')
    finally:
        ref.close()


def run_represented(output):
    """Genuine density/K-dressed pullbacks; no old response recomputation.

    Terminal0,.5,1 duals at zero are required for the five ordered noises.
    Terminal1 is integrated through .5 and zero but its old response entries
    (1,0),(1,1) in physical times are never evaluated or replaced.
    """
    output = _output(output)
    c, baseline, ref, bindings = _inputs()
    channel = RepresentedElectronChannel(c, baseline)
    dimension = c.d+c.f
    zeros = []
    responses = np.empty((5, dimension, dimension), complex)
    # Strictly future-source response is exactly causal zero, not missing data.
    for k, (a, b) in enumerate(NEW_PAIRS):
        if a < b:
            responses[k] = 0
    work = []
    sensitivity = None
    try:
        for a, terminal in enumerate(TIMES):
            at_zero = np.empty((*channel.shape, dimension), complex)
            batches = []
            needed = [(k, b) for k, (left, b) in enumerate(NEW_PAIRS)
                      if left == a and b <= a]
            for first in range(0, dimension, 19):
                indices = np.arange(first, min(first+19, dimension))
                pulled, receipt = channel.pullback(float(terminal), TIMES[:a+1], indices)
                at_zero[:, :, :, indices] = _finite(pulled[0.])
                for k, b in needed:
                    value = channel.response_from_dual(pulled[float(TIMES[b])], float(TIMES[b]))
                    responses[k, indices, :] = _finite(value)
                    if a == 1 and b == 0 and first == 0:
                        # One small bilinear-orientation comparison at this boundary.
                        coordinates = np.zeros((dimension, 1), complex)
                        coordinates[0, 0] = 1
                        forced = channel.forcing(0., coordinates)
                        dual = pulled[0.][:, :, :, :1]
                        bilinear = np.sum(dual * forced)
                        wrong_adjoint = np.sum(dual.conj() * forced)
                        sensitivity = dict(terminal=0.5, source_time=0., output_coordinate=0,
                            input_coordinate=0, comparison='response versus explicit bilinear forcing contraction; conjugated-dual orientation sensitivity',
                            response=_scalar(value[0, 0]), bilinear_action=_scalar(bilinear),
                            identity_defect=float(abs(value[0, 0] - bilinear)),
                            affected_orientation_difference=float(abs(bilinear - wrong_adjoint)))
                batches.append(receipt)
                del pulled
            zeros.append(at_zero)
            work.append(dict(terminal=float(terminal), batches=batches,
                             evaluated_response_sources=[float(TIMES[b]) for _, b in needed]))
        ordered = np.asarray([_finite(channel.ordered_noise(zeros[a], zeros[b]))
                              for a, b in NEW_PAIRS])
        positions = {pair: k for k, pair in enumerate(NEW_PAIRS)}
        reverse = np.asarray([ordered[positions[b, a]].T for a, b in NEW_PAIRS])
        new = dict(pi_retarded=_finite(responses), pi_greater=-1j*ordered,
                   pi_lesser=-1j*reverse, ordered_ab=ordered, ordered_ba=reverse)
        old, old_binding = _old('represented')
        try:
            arrays = {key: _merge(old[key], value) for key, value in new.items()}
        finally:
            old.close()
        path = output / 'represented.npz'
        np.savez(path, **arrays)
        return _finish(output, 'represented', path, c, ref, bindings, old_binding, sensitivity,
            external_q=[0], joint_coordinates=dimension, work=work,
            new_response_pairs=[[a, b] for a, b in NEW_PAIRS if a >= b],
            causal_new_zero_pairs=[[a, b] for a, b in NEW_PAIRS if a < b],
            pullback_terminals=TIMES.tolist(), persisted_duals=False,
            response_scope='A50 RepresentedElectronChannel with genuine moving-baseline density->K->delta_h tangent/transpose; all128 slots and all114 input/output coordinates',
            consumer_combination='Pi insertion is new_joint_Pi minus this represented Pi in identical coordinates; never subtract represented Sigma',
            cutover=CUTOVER)
    finally:
        ref.close()


def run(output):
    """Optional serial orchestration; native scheduler may use each entry alone."""
    output = _output(output)
    return dict(ei=run_ei(output/'ei'), joint=run_joint(output/'joint'),
                represented=run_represented(output/'represented'))
