"""Original full-PW uniform electronic currents, with no H/force action.

The existing physical field vertex includes kinetic AND physical-UPF projector
terms. Its checked real-direction Hermiticity permits one owned action array
per axis to serve both legs, without changing the original frame, D or Gram.
"""
from pathlib import Path
import resource
import sys
from time import perf_counter, process_time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src')]+[str(HERE.parent/name) for name in ('0030', '0029', '0028')]


def prepare_uniform_current_inputs(initial_result, coefficient_manifest, original_record):
    """Return the existing FullNoiseProbeInputs for three actual J=-L/2 probes."""
    import numpy as np
    from full_noise_probe_inputs import FullNoiseProbeInputs
    from original_noise_source import restore_original_noise_source, _identity
    import matched_reference_io as material_io
    from substrate_framework.material_physical_operators import PhysicalFieldVertices
    from substrate_framework.physical_factor_preparation import ActionFactors
    from substrate_framework.qe_material_hamiltonian import MaterialHamiltonian
    from substrate_framework.streamed_material_state import file_binding

    started, cpu_started = perf_counter(), process_time()
    context = restore_original_noise_source(initial_result, coefficient_manifest, original_record)
    # Source restoration exposes its exact physical factors, but not background.
    # Reload the SAME saved background API, never reconstruct H/factors or replay
    # preparation. MaterialHamiltonian construction checks their actual binding;
    # neither its apply method nor any force/trajectory API is invoked here.
    background_started = perf_counter()
    _, registration = material_io.registration()
    control, background, ledger, rows, restored = material_io.load_case(context.source['case'], registration)
    if (control.provenance() != context.source['control_provenance']
            or restored['background_receipt'] != context.source['saved_background_receipt']
            or restored['population_ledger'] != context.source['population_ledger']):
        raise ValueError('field background/control differs from admitted original source')
    del control, ledger, rows
    hamiltonian = MaterialHamiltonian(background, context.coupling.factors, context.source['spin'], 8)
    vertices = PhysicalFieldVertices(hamiltonian, projector_block_size=64)
    background_seconds = perf_counter()-background_started

    probes, specifications = {}, {}
    preparation_seconds = action_seconds = 0.
    for axis, name in enumerate(('current_x', 'current_y', 'current_z')):
        direction = np.zeros(3)
        direction[axis] = 1.
        t = perf_counter()
        prepared = vertices.prepare(direction)
        preparation_seconds += perf_counter()-t
        t = perf_counter()
        # PhysicalFieldAction.apply allocates a new owned full-row output.
        # Scale ONLY that output in place, never the borrowed original frame.
        action = prepared.apply(context.psi)
        action *= -0.5
        action.setflags(write=False)
        action_seconds += perf_counter()-t
        specification = dict(kind='physical-uniform-electronic-current',
            axis='xyz'[axis], field_direction=direction.tolist(), external_q=[0, 0, 0],
            time_inverse_ry=0, units='electronic_current_hartree', current_scale=-0.5,
            field_api='PhysicalFieldVertices.prepare(real_direction).apply',
            adjoint_api='PhysicalFieldAction.adjoint_apply delegates apply',
            adjoint_reuse='real direction and checked original Hermitian UPF D; J_dagger=J',
            basis_identity=context.pw_basis_identity, source=context.source['source'],
            physical_sources=list(context.operator_sources), quadrature_applied=False)
        identity = _identity(specification)
        # Exactly the same array is correct on the equal Hilbert-adjoint leg;
        # this is three physical evaluations, NOT six independent evaluations.
        probe = ActionFactors(action, action, identity, identity,
            'electronic_current_hartree', context.pw_basis_identity, context.state_authority)
        probe.validate(context.psi.shape, context.pw_basis_identity)
        probes[name], specifications[name] = probe, specification
        del prepared

    source = dict(context.source)
    source['operators'] = specifications
    source['uniform_current_implementation'] = [file_binding(Path(__file__))]
    work = dict(context.operator_work)
    work.update(current_action_evaluations=3, adjoint_reuses_by_hermiticity=3,
        force_forward_actions=0, force_adjoint_actions=0, h0_actions=0,
        trajectory_actions=0, native_calls=0, gram_recomputations=0,
        population_refills=0, preparation_replays=0,
        field_background_restorations=1, field_background_seconds=background_seconds,
        current_preparation_seconds=preparation_seconds, current_action_seconds=action_seconds,
        current_action_columns=3*context.psi.shape[1],
        elapsed_seconds=perf_counter()-started, process_cpu_seconds=process_time()-cpu_started,
        peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*1024)
    return FullNoiseProbeInputs(context.psi, context.gram, context.density_coefficients,
        context.active_coordinates, context.selected, probes, context.factor_identity,
        context.pw_basis_identity, context.state_authority, context.frame, source,
        context.quadrature_weight, work)
