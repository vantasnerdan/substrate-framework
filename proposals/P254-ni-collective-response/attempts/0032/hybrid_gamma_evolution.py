"""Checkpoint-only NEW moving Gamma mean and nonconjugating Gaussian duals.

No old adaptive history is recovered. Only requested mean boundaries are retained;
the reverse mean is never reset to (or repaired with) a forward checkpoint.
"""
from dataclasses import dataclass
from pathlib import Path
import os
import sys
from time import monotonic, perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0029')]

import numpy as np
from scipy.integrate import DOP853
from reduced_ode_initial_step import bounded_first_step
from reduced_coherence_readout import normalized_ordered
from substrate_framework.streamed_material_state import file_binding


TIMES = (0., 1., 10., 30.)
RECOVERY = ('Actual interaction coordinates and completed forward/backward boundaries '
            'survive failure with input-bound receipts. No adaptive/dense step history '
            'is saved; no automatic restart or free-resume claim. This runner accepts fresh output.')


def save_npz_exclusive(path, **arrays):
    """Never overwrite a scientific artifact; only receipts admit complete files."""
    path = Path(path)
    with path.open('xb') as stream:
        np.savez(stream, **arrays)
        stream.flush()
        os.fsync(stream.fileno())
    return file_binding(path)


def _advance(rhs, start, stop, initial, *, rtol, atol, stage, progress):
    """DOP853 with the original bounded seed and no retained adaptive history."""
    first = bounded_first_step(rhs, (start, stop), initial, rtol=rtol, atol=atol)
    solver = DOP853(rhs, start, initial, stop, rtol=rtol, atol=atol, first_step=first)
    steps, last_notice, started = 0, monotonic(), perf_counter()
    while solver.status == 'running':
        message = solver.step()
        steps += 1
        if solver.status == 'failed':
            raise RuntimeError(stage + ': ' + str(message))
        if not np.isfinite(solver.y).all():
            raise FloatingPointError(stage + ': nonfinite raw state; no clipping')
        if monotonic() - last_notice >= 5:
            progress(stage, time_inverse_ry=float(solver.t), accepted_steps=steps,
                     rhs_evaluations=int(solver.nfev))
            last_notice = monotonic()
    stats = dict(interval_inverse_ry=[float(start), float(stop)], accepted_steps=steps,
                 rhs_evaluations=int(solver.nfev), first_step=first,
                 elapsed_seconds=perf_counter()-started, method='DOP853')
    return solver.y.copy(), stats


def _difference(reconstructed, forward, nx):
    delta = reconstructed-forward
    result = {}
    for name, section in (('whole_state', slice(None)), ('mechanical', slice(0, nx)),
                          ('electronic', slice(nx, None))):
        value, reference = delta[section], forward[section]
        norm = float(np.linalg.norm(value))
        reference_norm = float(np.linalg.norm(reference))
        result[name] = dict(max_abs=float(np.max(np.abs(value))) if value.size else 0.,
                            l2=norm, forward_l2=reference_norm,
                            relative_l2=norm/reference_norm if reference_norm else None)
    return delta, result


@dataclass(frozen=True)
class HybridGammaEvolution:
    means: np.ndarray
    preparation_duals: np.ndarray
    ordered: np.ndarray
    coherence: np.ndarray
    artifacts: dict
    diagnostics: dict


def evolve_hybrid_gamma(model, coordinate, *, output, rtol, atol, progress,
                        publish_json):
    """One forward mean, then one continuous reverse mean with probe injections.

    Probe columns at preparation are chronological 0/1/10/30. The physical
    observable is Gamma Ni1/x, not a spatial Fourier sum. Ordered covariance is
    model.initial_ordered(z_i,z_j), with no extra spin/cell/root-N factor.
    """
    directory = Path(output)
    if not (np.isfinite(rtol) and rtol > 0 and np.isfinite(atol) and atol > 0):
        raise ValueError('finite positive original DOP853 tolerances required')
    initial = np.asarray(model.initial, dtype=complex)
    if initial.shape != (model.size,) or not np.isfinite(initial).all():
        raise ValueError('finite model initial state of declared size required')
    probe = np.asarray(model.probe(coordinate), dtype=complex)
    if probe.shape != initial.shape or not np.isfinite(probe).all():
        raise ValueError('finite original canonical displacement probe required')
    means = np.empty((len(TIMES), model.size), dtype=complex)
    means[0] = initial
    forward, checkpoints = [], []
    checkpoints.append(save_npz_exclusive(directory/'mean-000.npz',
                       time_inverse_ry=np.asarray(0.), mean=means[0]))
    progress('forward-start', time_inverse_ry=0., state_size=int(model.size))
    for index in range(1, len(TIMES)):
        value, stats = _advance(model.mean_rhs, TIMES[index-1], TIMES[index],
                                means[index-1], rtol=rtol, atol=atol,
                                stage='forward', progress=progress)
        means[index] = value
        forward.append(stats)
        checkpoints.append(save_npz_exclusive(directory/('mean-%03d.npz' % index),
                           time_inverse_ry=np.asarray(TIMES[index]), mean=value))
        progress('forward-checkpoint', time_inverse_ry=TIMES[index],
                 checkpoint=checkpoints[-1], solver=stats)
    baseline_archive = save_npz_exclusive(directory/'baseline-complete.npz',
                        times_inverse_ry=np.asarray(TIMES), means=means,
                        terminal_mean=means[-1])
    baseline_receipt = publish_json(directory/'baseline-complete.json', dict(
        schema='P254-A32-hybrid-baseline-complete-v1', status='baseline-complete',
        archive=baseline_archive, checkpoints=checkpoints, solver=forward,
        history='NEW hybrid forward mean, not original A29 adaptive history',
        dense_history_saved=False, recovery=RECOVERY))
    progress('baseline-complete', archive=baseline_archive, receipt=baseline_receipt)

    # Only already-injected probes are integrated. Mean follows f(t,y), duals
    # follow -Df(t,y).T z along decreasing physical time. No conjugation.
    mean = means[-1].copy()
    active = probe[:, None].copy()  # time 30
    column_times = [TIMES[-1]]
    reverse, reconstruction, raw_differences = [], [], []
    backward_checkpoints = []
    for index in range(len(TIMES)-1, 0, -1):
        count = active.shape[1]
        shape = (model.size, count)
        augmented = np.concatenate((mean, active.ravel()))

        def rhs(time, value):
            current = value[:model.size]
            dual = value[model.size:].reshape(shape)
            derivative = np.empty_like(value)
            derivative[:model.size] = model.mean_rhs(time, current)
            derivative[model.size:] = -model.transpose(time, current, dual).ravel()
            return derivative

        stop = TIMES[index-1]
        value, stats = _advance(rhs, TIMES[index], stop, augmented,
                                rtol=rtol, atol=atol, stage='backward', progress=progress)
        mean = value[:model.size].copy()
        active = value[model.size:].reshape(shape).copy()
        delta, defects = _difference(mean, means[index-1], model.nx)
        raw_differences.append(delta)
        reconstruction.append(dict(time_inverse_ry=stop, **defects))
        reverse.append(stats)
        boundary = save_npz_exclusive(directory/('backward-%03d.npz' % (index-1)),
                    time_inverse_ry=np.asarray(stop), reconstructed_mean=mean,
                    mean_minus_forward=delta, duals=active,
                    probe_times_inverse_ry=np.asarray(column_times))
        receipt = publish_json(directory/('backward-%03d.json' % (index-1)), dict(
            schema='P254-A32-backward-boundary-v1', status='backward-boundary-complete',
            archive=boundary, time_inverse_ry=stop, probe_times_inverse_ry=column_times,
            reconstruction=defects, solver=stats))
        backward_checkpoints.append(dict(archive=boundary, receipt=receipt))
        progress('backward-checkpoint', time_inverse_ry=stop,
                 checkpoint=boundary, receipt=receipt, reconstruction=defects, solver=stats)
        # At 10,1 inject next terminal observable; at preparation add t=0.
        active = np.column_stack((active, probe))
        column_times.append(stop)
    order = np.argsort(column_times)
    duals = active[:, order]
    ordered = np.asarray(model.initial_ordered(duals, duals), dtype=complex)
    if ordered.shape != (len(TIMES), len(TIMES)) or not np.isfinite(ordered).all():
        raise FloatingPointError('finite ordered 0/1/10/30 covariance required; no repair')
    variances = np.diag(ordered).copy()
    roots = np.sqrt(variances)
    denominators = roots[:, None]*roots[None, :]
    valid = np.empty(ordered.shape, dtype=bool)
    coherence = np.full_like(ordered, complex(np.nan, np.nan))
    normalization_records = []
    for left in range(len(TIMES)):
        for right in range(len(TIMES)):
            record = normalized_ordered(ordered[left, right], variances[left], variances[right])
            valid[left, right] = record['defined']
            if record['defined']:
                coherence[left, right] = record['complex_degree']
            normalization_records.append(dict(left_time=TIMES[left], right_time=TIMES[right], **record))
    # Defects are raw, not criteria used to Hermitize or relabel the state.
    diagnostics = dict(
        forward_solver=forward, backward_solver=reverse,
        reconstruction=reconstruction, reconstruction_repaired=False,
        ordered_hermiticity_max_abs=float(np.max(np.abs(ordered-ordered.conj().T))),
        ordered_variance_real=variances.real, ordered_variance_imag=variances.imag,
        variance_imaginary_max_abs=float(np.max(np.abs(variances.imag))),
        variance_negative_real_indices=np.flatnonzero(variances.real < 0),
        normalized_defined=valid, normalization='existing normalized_ordered; separate principal roots of raw complex variances',
        normalization_records=normalization_records,
        dual_equation='dz/dt=-L(t, reconstructed_mean).T z, true nonconjugating transpose',
        mean_history='NEW forward hybrid checkpoints; backward differences are raw reconstruction error, not repaired original history',
        dense_history_saved=False, recovery=RECOVERY)
    readout = save_npz_exclusive(directory/'hybrid-gamma-readout.npz',
        times_inverse_ry=np.asarray(TIMES), preparation_duals=duals,
        ordered_covariance=ordered, ordered_variances=variances,
        normalized_gamma_coherence=coherence,
        coherence_abs=np.abs(coherence), coherence_phase=np.angle(coherence),
        normalized_defined=valid, normalization_denominators=denominators,
        reconstructed_minus_forward=np.asarray(raw_differences),
        reconstruction_times_inverse_ry=np.asarray(TIMES[-2::-1]),
        mean_mechanical=means[:, :model.nx])
    progress('ordered-readout-complete', archive=readout)
    return HybridGammaEvolution(means, duals, ordered, coherence,
        dict(baseline=baseline_archive, baseline_receipt=baseline_receipt,
             forward_checkpoints=checkpoints, backward_checkpoints=backward_checkpoints,
             readout=readout), diagnostics)
