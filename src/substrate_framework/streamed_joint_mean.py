"""Self-consistent electronic/ionic/Maxwell exponential midpoint and history.

Only homogeneous mean evolution is implemented here, not local-q tangents or
quantum noise. A converged nonlinear defect is not a time-accuracy certificate.
"""
from pathlib import Path

import numpy as np

from .streamed_joint_reduction import JointReduction
from .streamed_material_evolution import exponential_action


class JointConvergenceError(RuntimeError):
    def __init__(self, residuals):
        self.residuals = residuals
        super().__init__('joint midpoint failed; no successor adopted: '+repr(residuals))


class StreamedJointMean(JointReduction):
    def __init__(self, *args, exponential=exponential_action, **kwargs):
        super().__init__(*args, exponential=exponential, **kwargs)

    @staticmethod
    def midpoint_coordinates(state, stage, dt):
        return tuple(initial+dt*rhs/2 for initial, rhs in zip(
            (state.u, state.p, state.coordinate, state.momentum),
            (stage.velocity, stage.force, stage.field_velocity, stage.field_force)))

    @staticmethod
    def residual(old, updated):
        names = ('potential', 'gradient', 'field_gradient_electronic', 'u', 'p',
                 'coordinate', 'momentum', 'velocity', 'force', 'field_velocity', 'field_force')
        return {name: float(np.max(np.abs(getattr(updated, name)-getattr(old, name)) /
                                  (1+np.abs(getattr(updated, name))), initial=0.)) for name in names}

    def step(self, state, directory, *, dt, max_iterations, tolerance):
        if (not np.isfinite(dt) or dt <= 0 or type(max_iterations) is not int
                or max_iterations < 1 or not 0 < tolerance < 1):
            raise ValueError('positive interval and bounded nonlinear solve required')
        directory = Path(directory).resolve()
        directory.mkdir(parents=True, exist_ok=False)
        stage = self.reduce(state)
        residuals = []
        for iteration in range(max_iterations):
            updated = self.reduce(state, frozen_stage=stage, interval=dt/2,
                                  coordinates=self.midpoint_coordinates(state, stage, dt))
            residuals.append(self.residual(stage, updated))
            if max(residuals[-1].values()) <= tolerance:
                break
            stage = updated
        else:
            self.budget.json(directory/'nonconvergence.json', dict(residuals=residuals,
                             prior_boundary=str(state.path), completed=False))
            raise JointConvergenceError(residuals)
        # Keep precisely the tested frozen stage; do not replace it with an
        # untested mixture of new coordinates and stale electronic reductions.
        end = tuple(2*middle-initial for initial, middle in zip(
            (state.u, state.p, state.coordinate, state.momentum),
            self.midpoint_coordinates(state, stage, dt)))
        from .streamed_material_state import file_binding
        record = stage.save(directory/'frozen-midpoint', self.budget, metadata=dict(
            kind='actual-frozen-midpoint', time_ry_inverse=state.time+dt/2,
            start_time_ry_inverse=state.time, dt_ry_inverse=dt,
            prior_boundary=str(state.path), midpoint_iterations=iteration+1,
            prior_boundary_binding=file_binding(state.path),
            midpoint_residuals=residuals, tolerance=tolerance,
            operator='physical H0 + K(delta n) + physical g(u) + field-only L(a)',
            frame_reconstruction='exp(-i dt/2 H_frozen) applied to prior original-column frames'))
        def frames():
            for slot in state.slots:
                psi, _ = state.arrays(slot)
                with self.operators(slot) as (coupling, h0, vertex):
                    yield self.propagate(coupling, h0, vertex, slot, stage, psi, dt)
        result = state.successor(directory/'completed', frames(), time=state.time+dt,
            u=end[0], p=end[1], coordinate=end[2], momentum=end[3], stage_record=record,
            evidence=dict(method='all-slot self-consistent joint exponential midpoint',
                exponential=getattr(self.exponential, '__name__', type(self.exponential).__name__),
                prior_boundary=str(state.path), complete_reductions=self.reductions,
                scope='homogeneous coupled mean; whole local-q dynamics separate'), budget=self.budget)
        return result

    def evolve(self, state, directory, *, physical_interval, steps, max_steps,
               max_iterations, tolerance, refinement_label):
        if (type(steps) is not int or type(max_steps) is not int or not 1 <= steps <= max_steps
                or not np.isfinite(physical_interval) or physical_interval <= 0 or not refinement_label):
            raise ValueError('bounded steps, positive physical interval and refinement label required')
        directory = Path(directory).resolve()
        directory.mkdir(parents=True, exist_ok=False)
        records = []
        for index in range(steps):
            state = self.step(state, directory/('step-%06d' % (index+1)),
                dt=physical_interval/steps, max_iterations=max_iterations, tolerance=tolerance)
            # Publish a durable completed-boundary pointer before diagnostic work
            # can exhaust a cap. It is not a claim the requested interval finished.
            self.budget.json(directory/('completed-%06d.json' % (index+1)), dict(
                boundary=str(state.path), time_ry_inverse=state.time, requested_steps=steps,
                completed_steps=index+1, refinement_label=refinement_label))
            endpoint = self.reduce(state)
            record = endpoint.save(directory/('endpoint-%06d' % (index+1)), self.budget,
                metadata=dict(kind='complete-endpoint-reduction', boundary=str(state.path),
                              time_ry_inverse=state.time, refinement_label=refinement_label))
            self.budget.json(directory/('diagnostics-%06d.json' % (index+1)), record)
            records.append(record)
        return state, records
