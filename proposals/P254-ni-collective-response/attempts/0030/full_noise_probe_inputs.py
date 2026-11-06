"""Original full-PW Cartesian/modal force factors, without H or state replay."""
from dataclasses import dataclass
import resource
from time import perf_counter, process_time

import numpy as np
from original_noise_source import restore_original_noise_source, _identity
from substrate_framework.physical_factor_preparation import ActionFactors


@dataclass(frozen=True)
class FullNoiseProbeInputs:
    psi: np.ndarray
    gram: np.ndarray
    density_coefficients: np.ndarray
    active_coordinates: np.ndarray
    selected: np.ndarray
    probes: dict
    factor_identity: str
    pw_basis_identity: str
    state_authority: str
    frame: dict
    source: dict
    quadrature_weight: float
    operator_work: dict

    def algebra_arguments(self):
        """Keyword arguments only; do not construct another full action bank."""
        return dict(psi=self.psi, gram=self.gram,
                    density_coefficients=self.density_coefficients,
                    active_coordinates=self.active_coordinates, probes=self.probes,
                    factor_identity=self.factor_identity,
                    pw_basis_identity=self.pw_basis_identity,
                    state_authority=self.state_authority)


def _force_actions(prepared, psi, block_size=8):
    action = np.empty(psi.shape, dtype=np.complex128)
    adjoint = np.empty(psi.shape, dtype=np.complex128)
    forward_seconds = adjoint_seconds = 0.
    blocks = 0
    for start in range(0, psi.shape[1], block_size):
        stop = min(start + block_size, psi.shape[1])
        t = perf_counter()
        action[:, start:stop] = prepared.apply(psi[:, start:stop])
        forward_seconds += perf_counter() - t
        t = perf_counter()
        adjoint[:, start:stop] = prepared.adjoint_apply(psi[:, start:stop])
        adjoint_seconds += perf_counter() - t
        blocks += 1
    # The operator API is the POSITIVE derivative. Force is its negative,
    # independently on both legs; no spin, cell, or quadrature scaling.
    action *= -1
    adjoint *= -1
    action.setflags(write=False)
    adjoint.setflags(write=False)
    return action, adjoint, dict(force_forward_actions=1, force_adjoint_actions=1,
        forward_block_calls=blocks, adjoint_block_calls=blocks,
        forward_action_columns=psi.shape[1], adjoint_action_columns=psi.shape[1],
        forward_seconds=forward_seconds, adjoint_seconds=adjoint_seconds,
        block_size=block_size)




def _force_result(context, directions, specifications):
    started, cpu_started = perf_counter(), process_time()
    probes, work_rows = {}, []
    for name, direction in directions.items():
        prepared = context.coupling.ionic.prepare(direction)
        action, adjoint, work = _force_actions(prepared, context.psi)
        identity = _identity(specifications[name])
        probe = ActionFactors(action, adjoint, identity,
            identity+':actual-Hilbert-adjoint', 'Ry/bohr',
            context.pw_basis_identity, context.state_authority)
        probe.validate(context.psi.shape, context.pw_basis_identity)
        probes[name] = probe
        work_rows.append(work)
    source = dict(context.source)
    if len(specifications) == 1:
        source['operator'] = next(iter(specifications.values()))
    else:
        source['operators'] = specifications
    work = dict(context.operator_work)
    work.update(force_forward_actions=sum(row['force_forward_actions'] for row in work_rows),
        force_adjoint_actions=sum(row['force_adjoint_actions'] for row in work_rows),
        forward_block_calls=sum(row['forward_block_calls'] for row in work_rows),
        adjoint_block_calls=sum(row['adjoint_block_calls'] for row in work_rows),
        forward_action_columns=sum(row['forward_action_columns'] for row in work_rows),
        adjoint_action_columns=sum(row['adjoint_action_columns'] for row in work_rows),
        forward_seconds=sum(row['forward_seconds'] for row in work_rows),
        adjoint_seconds=sum(row['adjoint_seconds'] for row in work_rows),
        elapsed_seconds=work['elapsed_seconds']+perf_counter()-started,
        process_cpu_seconds=work['process_cpu_seconds']+process_time()-cpu_started,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        block_size=8)
    return FullNoiseProbeInputs(context.psi, context.gram, context.density_coefficients,
        context.active_coordinates, context.selected, probes, context.factor_identity,
        context.pw_basis_identity, context.state_authority, context.frame, source,
        context.quadrature_weight, work)


def _operator_spec(context):
    return dict(kind='minus-current-physical-ionic-Hamiltonian-derivative',
        external_q=[0, 0, 0], time_inverse_ry=0, units='Ry/bohr',
        basis_identity=context.pw_basis_identity, source=context.source['source'],
        physical_sources=list(context.operator_sources), force_sign=-1,
        derivative_scale=1, quadrature_applied=False,
        adjoint_api='BlochDirectionalAction.adjoint_apply')


def prepare_full_noise_probe(initial_result, coefficient_manifest, original_record, atom, axis, bands=12):
    """Restore once and evaluate only the requested original Cartesian probe."""
    if type(atom) is not int or axis not in ('x', 'y', 'z'):
        raise ValueError('one exact original atom and Cartesian axis required')
    context = restore_original_noise_source(initial_result, coefficient_manifest, original_record, bands)
    if not 1 <= atom <= len(context.atom_names):
        raise ValueError('original atom outside the physical cell')
    direction = np.zeros((len(context.atom_names), 3))
    direction[atom-1, 'xyz'.index(axis)] = 1.
    name = '%s%d/%sGamma' % (context.atom_names[atom-1], atom, axis)
    specification = dict(_operator_spec(context), atom=atom,
        species=context.atom_names[atom-1], axis=axis,
        derivative_direction=direction.tolist())
    return _force_result(context, {name:direction}, {name:specification})


def prepare_full_noise_modes(initial_result, coefficient_manifest, original_record, directions):
    """Two complex physical modes; no unrelated Cartesian or H action."""
    context = restore_original_noise_source(initial_result, coefficient_manifest, original_record)
    directions = np.asarray(directions, dtype=complex)
    if directions.shape != (3*len(context.atom_names), 2) or not np.isfinite(directions).all():
        raise ValueError('two finite original Cartesian modal directions required')
    probes, specifications = {}, {}
    for mode in range(2):
        name = 'mode'+str(mode)
        direction = directions[:, mode].reshape(len(context.atom_names), 3)
        probes[name] = direction
        specifications[name] = dict(_operator_spec(context), mode=mode,
            derivative_direction_real=direction.real.tolist(),
            derivative_direction_imag=direction.imag.tolist(),
            convention='sum_i direction[i,mode] F_i; no conjugation; independent physical adjoint')
    return _force_result(context, probes, specifications)
