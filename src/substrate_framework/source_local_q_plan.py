"""Payload and replay-action accounting for exact context-only q pairs.

Physical moving frames are still required.  Cached-frame aliases count only
when the caller has separately authenticated their actual saved boundaries.
"""
import json

from .partitioned_q_response import frame_times, partition_plan, initial_partition_plan
from .streamed_material_state import authenticate


def _plan(base, times, second_times, aliases, edges, columns):
    required = frame_times(times) | frame_times(second_times)
    aliases = set(float(time) for time in aliases)
    new_frames = required-aliases
    intervals = [len(times)-1, len(second_times)-1]
    # For reverse interval c, each moment stage replays all 4c previous
    # derivatives, then its own prefix and current free/moment operation.
    # Final source contraction replays four derivatives per interval once.
    free_calls = sum(8*count*count+2*count for count in intervals)
    moment_calls = 4*sum(intervals)
    potentials = 4*sum(intervals)*base['complete_potential_bytes']
    moments = base['independent_moments_bytes'] if sum(intervals) else 0
    frames = len(new_frames)*base['actual_frame_bytes_per_time']
    source_scalars = 48*edges  # Two complex ordered values and two float64 norms.
    return dict(execution_scheme='source-local-discrete-rk4-v1', mode='pair',
        full_q_state_bytes_by_geometry=base['one_q_state_bytes'],
        retained_full_q_states=0, full_q_payload_bytes_allocated=0,
        all_required_frame_times=sorted(required),
        genuine_boundary_alias_times=sorted(required & aliases),
        owned_frame_times=sorted(new_frames), frame_payload_upper_bound=frames,
        moment_payload_upper_bound=moments,
        retained_kernel_context_payload_bytes=potentials,
        retained_source_scalar_payload_bytes=source_scalars,
        successful_peak_payload_upper_bound=frames+moments+potentials+source_scalars,
        final_payload_upper_bound=potentials+source_scalars,
        edge_generator_calls=(free_calls+moment_calls)*edges,
        free_generator_calls=free_calls*edges,
        moment_calls=moment_calls*edges,
        bare_hamiltonian_actions=2*free_calls*edges,
        bare_hamiltonian_columns=free_calls*columns,
        limitations='Payload only: excludes measured JSON/header/reservation slack and failed attempts. Frame reconstruction has additional charged actions. Native memory and runtime require a representative physical measurement; this is not a production grant or an ETA.',
        tradeoff='Recompute full original edge equations from durable global contexts; no Q-sized input, derivative, accumulator, seed or retained left output is written.')


def source_local_plan(response, times, second_times, workers, *,
                      aliased_frame_times=(), dynamics='screened-joint'):
    base = partition_plan(response, times, second_times, workers, mode='pair', dynamics=dynamics)
    columns = sum(response.sectors[e.source].shape[1]+response.sectors[e.target].shape[1]
                  for e in response.edges)
    return _plan(base, times, second_times, aliased_frame_times, len(response.edges), columns)


def source_local_initial_plan(binding, *, q, mesh, density_shape, times,
                              second_times, workers, dynamics='screened-joint'):
    base = initial_partition_plan(binding, q=q, mesh=mesh, density_shape=density_shape,
        times=times, second_times=second_times, workers=workers, mode='pair', dynamics=dynamics)
    initial = json.loads(authenticate(binding).read_text())
    columns = 2*sum(slot['shape'][1] for slot in initial['slots'])
    result = _plan(base, times, second_times, (0.,), len(initial['slots']), columns)
    result.update(initial_metadata=binding, model_identity=base['model_identity'],
        q=base['q'], mesh=base['mesh'], density_shape=base['density_shape'],
        numerical_work_executed=False, history_available=False,
        proposed_time_grid=list(times), proposed_second_time_grid=list(second_times))
    return result
