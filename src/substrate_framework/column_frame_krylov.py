"""Bounded independent Arnoldi spaces for every original full-PW frame column.

Columns share batched H calls, never a cross-column projected matrix. Original
normalization, phase, ordering and D are unchanged; even dependent/zero columns
remain present. Frozen-operator identity authentication is the caller's duty.
"""
from time import perf_counter

import numpy as np
from scipy.linalg import expm

from .column_krylov_bounds import column_envelope


def memory_bound_bytes(rows, columns, dimension):
    """Conservative numerical-array plan, not a process RSS guarantee.

    Two stored bases; twelve frame panels cover inputs, output, H batches and
    invariant/compare temporaries; eight single-column basis panels cover full
    defects/QR; small per-column matrices and expm work. Excludes operator caches,
    native allocator overhead and OS/disk caches. No growing (columns*k)^2 array.
    """
    return 16*(2*rows*columns*dimension+12*rows*columns
               +8*rows*dimension+32*columns*dimension*dimension)


class ColumnFrameKrylov:
    """Independent two-pass Arnoldi; all stopping limits are hard limits.

    initial_action may be reused only with matching nonempty operator identities.
    These labels do not authenticate a callable: the runtime must bind the full
    operator, frozen stage and exact original frame to the retained-action receipt.
    operator_norm_bound controls the seed's measured normalization defect; growth
    bounds both signed ambient semigroups. H-action evaluation errors themselves
    are not enclosed. max_dimension is PER ORIGINAL COLUMN, not total basis.
    """
    def __init__(self, action, frame, *, max_dimension, action_cap,
                 max_memory_bytes, ambient_growth_bound, operator_norm_bound,
                 operator_identity, initial_action=None, initial_action_identity=None,
                 breakdown_rtol=1e-13, shifts=None):
        frame = np.asarray(frame)
        if frame.ndim != 2 or min(frame.shape) < 1 or not np.isfinite(frame).all():
            raise ValueError('finite nonempty full raw frame required')
        n, m = frame.shape
        if (type(max_dimension) is not int or not 1 <= max_dimension <= n
                or type(action_cap) is not int or action_cap < 0
                or not np.isfinite([ambient_growth_bound, operator_norm_bound, breakdown_rtol]).all()
                or min(ambient_growth_bound, operator_norm_bound) < 0
                or not 0 <= breakdown_rtol < 1 or not operator_identity):
            raise ValueError('invalid dimension, action, growth, norm or identity limits')
        required = memory_bound_bytes(n, m, max_dimension)
        if required > max_memory_bytes:
            raise MemoryError(f'planned numerical arrays {required} exceed {max_memory_bytes}')
        if initial_action is not None:
            initial_action = np.asarray(initial_action)
            if (initial_action_identity != operator_identity or initial_action.shape != frame.shape
                    or not np.isfinite(initial_action).all()):
                raise ValueError('retained action must bind identical frozen operator and full frame')
        self.frame, self.action = frame, action
        self.max_dimension, self.action_cap = max_dimension, action_cap
        self.growth, self.operator_norm = float(ambient_growth_bound), float(operator_norm_bound)
        self.breakdown_rtol = breakdown_rtol
        self.operator_identity = operator_identity
        self.memory_plan = dict(array_workspace_bound_bytes=required, maximum_bytes=max_memory_bytes,
            excludes='operator caches/workspace, Python/native allocator and OS/disk caches')
        self.q = np.empty((m, max_dimension, n), dtype=complex)
        self.hq = np.empty_like(self.q)
        self.sizes = np.zeros(m, dtype=int)
        self.closed = np.zeros(m, dtype=bool)
        self.beta = np.linalg.norm(frame, axis=0)
        if not np.isfinite(self.beta).all():
            raise ValueError('column norm overflow')
        self.initial_errors = np.zeros(m)
        self.seed_defects = np.zeros(m)
        self.calls, self.requested_columns, self.generation = [], 0, 0
        active = np.flatnonzero(self.beta)
        for j in range(m):
            if self.beta[j] == 0:
                if np.any(frame[:, j] != 0):
                    raise ValueError('nonzero raw column norm underflow; no column dropped')
                self.closed[j] = True
            else:
                self.q[j, 0] = frame[:, j]/self.beta[j]
                self.sizes[j] = 1
                self.initial_errors[j] = np.linalg.norm(frame[:, j]-self.beta[j]*self.q[j, 0])
        if initial_action is None and active.size:
            applied = self._apply(self.q[active, 0].T.copy())
            for index, j in enumerate(active):
                self.hq[j, 0] = applied[:, index]
        elif initial_action is not None:
            for j in active:
                self.hq[j, 0] = initial_action[:, j]/self.beta[j]
                relation = np.linalg.norm(initial_action[:, j]-self.beta[j]*self.hq[j, 0])
                self.seed_defects[j] = (relation+self.operator_norm*self.initial_errors[j])/self.beta[j]
            if np.any(initial_action[:, self.closed] != 0):
                raise ValueError('zero raw column has nonzero retained H action')
        if shifts is None:
            self.shifts = np.array([float(np.vdot(self.q[j, 0], self.hq[j, 0]).real)
                                   if self.sizes[j] else 0. for j in range(m)])
        else:
            candidate = np.asarray(shifts)
            if candidate.shape != (m,) or not np.isrealobj(candidate) or not np.isfinite(candidate).all():
                raise ValueError('one finite real shift per original column required')
            self.shifts = candidate.astype(float, copy=True)
        self.seed_receipt = dict(raw_columns=m, nonzero_columns=int(active.size),
            retained_action_reused=initial_action is not None, operator_identity=operator_identity,
            reconstruction_column_norms=self.initial_errors.tolist(),
            seed_action_defect_column_norms=self.seed_defects.tolist(), shifts=self.shifts.tolist())
        self._stage = None

    def _apply(self, block):
        count = block.shape[1]
        if self.requested_columns+count > self.action_cap:
            raise RuntimeError('H-vector cap reached before action')
        self.requested_columns += count
        receipt = dict(columns=count, status='requested')
        self.calls.append(receipt)
        started = perf_counter()
        try:
            value = np.asarray(self.action(block))
            if value.shape != block.shape or not np.isfinite(value).all():
                raise ValueError('action returned invalid full-row panel')
            receipt['status'] = 'completed'
            return value
        finally:
            receipt['seconds'] = perf_counter()-started

    def _decomposition(self, j):
        k = int(self.sizes[j])
        q, hq = self.q[j, :k].T, self.hq[j, :k].T
        shifted = hq-self.shifts[j]*q
        projected = q.conj().T @ shifted
        residual = shifted-q @ projected
        return q, projected, residual

    def stage(self):
        if self._stage is None:
            small, diagnostics = [], []
            for j, k in enumerate(self.sizes):
                if not k:
                    small.append(np.empty((0, 0), dtype=complex))
                    diagnostics.append(dict(column=j, dimension=0, zero_original_column=True))
                    continue
                q, projected, residual = self._decomposition(j)
                small.append(projected)
                diagnostics.append(dict(column=j, dimension=int(k), breakdown=bool(self.closed[j]),
                    residual_frobenius=float(np.linalg.norm(residual)),
                    residual_column_norms=np.linalg.norm(residual, axis=0).tolist(),
                    initial_error=float(self.initial_errors[j]), seed_action_defect=float(self.seed_defects[j]),
                    orthogonality_frobenius=float(np.linalg.norm(q.conj().T @ q-np.eye(k))),
                    projected_antihermitian_frobenius=float(np.linalg.norm((projected-projected.conj().T)/(2j)))))
            self._stage = dict(projected=small, diagnostics=dict(generation=self.generation,
                original_columns=len(self.sizes), dimensions=self.sizes.tolist(),
                total_basis_vectors=int(self.sizes.sum()), columns=diagnostics,
                requested_h_vector_columns=self.requested_columns))
        return self._stage

    def propagate(self, time, *, bound_cells=4, bound_degree=4):
        if not np.isfinite(time) or bound_cells < 1 or bound_degree < 1:
            raise ValueError('finite time and positive envelope controls required')
        started = perf_counter()
        output = np.zeros(self.frame.shape, dtype=complex)
        records = []
        for j, k in enumerate(self.sizes):
            if not k:
                records.append(dict(column=j, absolute_error_envelope_estimate=0.))
                continue
            q, projected, residual = self._decomposition(j)
            small = self.beta[j]*expm(-1j*time*projected)[:, 0]
            output[:, j] = np.exp(-1j*time*self.shifts[j])*(q @ small)
            record = column_envelope(projected, residual, self.beta[j], self.initial_errors[j], time,
                ambient_growth=self.growth, seed_action_defect=self.seed_defects[j],
                cells=bound_cells, degree=bound_degree)
            records.append(dict(column=j, **record))
        # Frobenius aggregation across original columns is exact; no independence
        # of their physical states, orthogonality, or diagonal D is assumed.
        bounds = [r['absolute_error_envelope_estimate'] for r in records]
        envelope = float(np.hypot.reduce(bounds))
        return output, dict(time=float(time), absolute_error_envelope_estimate=envelope,
            columns=records, ambient_growth_bound=self.growth,
            bound_cells=bound_cells, bound_degree=bound_degree, seconds=perf_counter()-started,
            qualification='ordinary floating evaluation of full-row analytic envelopes; NOT interval certified')

    def expand(self):
        """One batched full-PW H call, at most one new direction per column.

        A small numerical breakdown stops expansion only, never asserts accuracy;
        its omitted direction remains in every measured full-row stage residual.
        """
        started = perf_counter()
        remaining = self.action_cap-self.requested_columns
        eligible = np.flatnonzero((self.sizes < self.max_dimension) & ~self.closed)
        if not eligible.size or remaining == 0:
            return dict(status='action_limit' if remaining == 0 and eligible.size else 'basis_or_breakdown_limit',
                        accepted_columns=0, requested_h_vector_columns=self.requested_columns)
        panel = np.empty((self.frame.shape[0], min(len(eligible), remaining)), dtype=complex)
        selected, candidates = [], []
        for j in eligible:
            if len(selected) == panel.shape[1]:
                break
            k = int(self.sizes[j])
            q = self.q[j, :k].T
            w = self.hq[j, k-1]-self.shifts[j]*self.q[j, k-1]
            scale = float(np.linalg.norm(w))
            for _ in range(2):
                w -= q @ (q.conj().T @ w)
            norm = float(np.linalg.norm(w))
            threshold = self.breakdown_rtol*scale
            closed = norm <= threshold or norm == 0.
            candidates.append(dict(column=int(j), residual_norm=norm, threshold=threshold, breakdown=closed))
            if closed:
                self.closed[j] = True
            else:
                panel[:, len(selected)] = w/norm
                selected.append(int(j))
        if selected:
            panel = panel[:, :len(selected)]
            applied = self._apply(panel)
            for index, j in enumerate(selected):
                k = int(self.sizes[j])
                self.q[j, k], self.hq[j, k] = panel[:, index], applied[:, index]
                self.sizes[j] += 1
            self.generation += 1
        self._stage = None
        return dict(status='expanded' if selected else 'numerical_breakdown',
                    accepted_columns=len(selected), selected_original_columns=selected,
                    candidates=candidates, dimensions=self.sizes.tolist(),
                    requested_h_vector_columns=self.requested_columns, seconds=perf_counter()-started)
