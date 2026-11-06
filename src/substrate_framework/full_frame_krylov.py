"""Full-row frozen-H block Krylov with raw-column reconstruction and residuals.

No PW-square operator is formed; no physical state or occupation is truncated.
The supplied action is used verbatim (including any caller's atomic counter).
Bounds below are exact-arithmetic formulas evaluated in ordinary floating point,
not interval certificates. See P254 A28/block-krylov-registration.json.
"""
from __future__ import annotations

from math import factorial
from time import perf_counter

import numpy as np
from scipy.linalg import expm, solve_triangular


def _norm(x):
    return float(np.linalg.norm(x))


def _exp(x):
    return float(np.exp(x)) if x < 709.0 else float('inf')


def memory_bound_bytes(rows, columns, basis):
    """Conservative array-workspace planning bound, excluding the H action.

    Includes stored Q/HQ, full defect, QR/SVD/reorthogonalization temporaries,
    original/output panels, and dense projected expm work. BLAS/LAPACK allocator
    overhead and the supplied operator's own caches/workspace are not bounded.
    """
    return 16 * (8 * rows * basis + 8 * rows * columns + 32 * basis * basis)


def physical_invariants(frame, raw_d):
    """Small Gram contractions of rho=frame D frame*, never a PW-square rho.

    D is neither diagonalized nor refilled. Complex traces expose raw defects.
    Density Hermiticity is measured with a thin QR, including nonorthogonal raw
    frames. These are diagnostics, not imposed repairs/conservation constraints.
    """
    gram = frame.conj().T @ frame
    dg = raw_d @ gram
    _, r = np.linalg.qr(frame, mode='reduced')
    density_small = r @ raw_d @ r.conj().T
    scalar = lambda z: dict(real=float(z.real), imag=float(z.imag))
    return dict(number=scalar(np.trace(dg)), purity=scalar(np.trace(dg @ dg)),
                frame_frobenius=_norm(frame), gram_frobenius=_norm(gram),
                density_frobenius=_norm(density_small),
                density_hermiticity_frobenius=_norm(density_small-density_small.conj().T),
                raw_D_hermiticity_frobenius=_norm(raw_d-raw_d.conj().T))


class FullFrameKrylov:
    """Bounded block Arnoldi workspace for one identical frozen operator.

    Initial QR retains every column; ill-conditioned reconstruction is rejected,
    not silently rank-truncated. Only expansion directions can be deflated. The
    original frame remains untouched, and its actual reconstruction residual is
    always part of the propagation estimate. `action_cap` counts newly requested
    H-vector columns, independent of an outer counter's call-based cap.
    """

    def __init__(self, action, frame, *, initial_action=None, max_blocks,
                 max_basis, action_cap, max_memory_bytes,
                 ambient_growth_bound, deflation_rtol=1e-13):
        frame = np.asarray(frame)
        if frame.ndim != 2 or not np.isfinite(frame).all():
            raise ValueError('finite full raw frame required')
        n, m = frame.shape
        if not (0 < m <= n and m <= max_basis <= n and max_blocks >= 1
                and action_cap >= 0 and np.isfinite(ambient_growth_bound)
                and ambient_growth_bound >= 0 and 0 <= deflation_rtol < 1):
            raise ValueError('invalid basis/action/growth/deflation limits')
        required = memory_bound_bytes(n, m, max_basis)
        if required > max_memory_bytes:
            raise MemoryError(f'planned numerical array workspace {required} exceeds cap')
        self.frame, self.action = frame, action
        self.max_blocks, self.max_basis, self.action_cap = max_blocks, max_basis, action_cap
        self.growth = float(ambient_growth_bound)
        self.deflation_rtol = deflation_rtol
        self.memory_plan = dict(array_workspace_bound_bytes=required,
                                maximum_bytes=max_memory_bytes,
                                excludes='operator caches/workspace, Python and native allocator overhead')
        self.calls = []
        self.requested_columns = 0
        self._q = np.empty((n, max_basis), dtype=complex, order='F')
        self._hq = np.empty_like(self._q)
        started = perf_counter()
        q, r = np.linalg.qr(frame, mode='reduced')
        condition = float(np.linalg.cond(r))
        if not np.isfinite(condition) or condition*np.finfo(float).eps > 1e-8:
            raise ValueError('initial raw reconstruction too ill-conditioned; no columns clipped')
        self.initial_r = r
        self.initial_qr_residual = frame-q @ r
        self._q[:, :m] = q
        if initial_action is None:
            hq = self._apply(q)
            reuse = dict(reused=False, columns=0)
        else:
            initial_action = np.asarray(initial_action)
            if initial_action.shape != frame.shape or not np.isfinite(initial_action).all():
                raise ValueError('retained H Psi must match every original row/column')
            # Transpose, not adjoint: HQ R=H Psi. No inverse and no seed H call.
            hq = solve_triangular(r.T, initial_action.T, lower=True).T
            reuse = dict(reused=True, columns=m,
                         triangular_solve_residual_frobenius=_norm(hq @ r-initial_action),
                         qualification='HQ inferred from retained H Psi; action/QR roundoff not interval enclosed')
        self._hq[:, :m] = hq
        self.size, self.blocks, self.last_start = m, 1, 0
        self.seed_receipt = dict(seconds=perf_counter()-started, raw_columns=m,
                                 qr_condition=condition,
                                 reconstruction_residual_frobenius=_norm(self.initial_qr_residual),
                                 retained_action=reuse)
        self.expansions = []
        self._stage = None

    @property
    def q(self):
        return self._q[:, :self.size]

    @property
    def hq(self):
        return self._hq[:, :self.size]

    def _apply(self, block):
        count = block.shape[1]
        if self.requested_columns+count > self.action_cap:
            raise RuntimeError('H-vector column cap reached before action')
        # Count requested work even if the caller rejects or fails the request.
        self.requested_columns += count
        receipt = dict(columns=count, status='requested')
        self.calls.append(receipt)
        start = perf_counter()
        try:
            result = np.asarray(self.action(block))
            if result.shape != block.shape or not np.isfinite(result).all():
                raise ValueError('action returned wrong shape or nonfinite full-row data')
            receipt['status'] = 'completed'
            return result
        finally:
            receipt['seconds'] = perf_counter()-start

    def stage(self):
        """Keep actual full-row defect and initial residual, not Ritz heuristics."""
        if self._stage is None:
            start = perf_counter()
            q = self.q
            t = q.conj().T @ self.hq
            c = q.conj().T @ self.frame
            residual = self.hq-q @ t
            initial = self.frame-q @ c
            self._stage = dict(projected=t, coefficients=c, residual=residual,
                               initial_residual=initial,
                               diagnostics=dict(basis=self.size, blocks=self.blocks,
                                   residual_frobenius=_norm(residual),
                                   residual_column_norms=np.linalg.norm(residual, axis=0).tolist(),
                                   initial_residual_frobenius=_norm(initial),
                                   orthogonality_frobenius=_norm(q.conj().T @ q-np.eye(self.size)),
                                   projected_antihermitian_frobenius=_norm((t-t.conj().T)/(2j)),
                                   requested_h_vector_columns=self.requested_columns,
                                   seconds=perf_counter()-start))
        return self._stage

    def propagate(self, time, *, bound_cells=4, bound_degree=4):
        """Return all raw columns plus an analytic residual-envelope estimate.

        Cell Taylor remainders are integrated analytically; no quadrature sample
        is labeled a bound. All non-Hermitian projected growth is included. The
        envelope is valid in exact arithmetic given the ambient growth bound.
        Floating-point action, QR, expm, eigensolver and norm errors are NOT
        enclosed, including cancellation in the residual Gram square root.
        """
        if not np.isfinite(time) or bound_cells < 1 or bound_degree < 1:
            raise ValueError('finite time and positive envelope controls required')
        start = perf_counter()
        s = self.stage()
        t, c, e = s['projected'], s['coefficients'], s['residual']
        tau, sign = abs(float(time)), (1.0 if time >= 0 else -1.0)
        b = -1j*sign*t
        small = expm(tau*b) @ c
        output = self.q @ small
        h = (b+b.conj().T)/2
        eigen = np.linalg.eigvalsh(h)
        projected_growth = max(0., float(eigen[-1]))
        two_sided_growth = float(np.max(np.abs(eigen)))
        initial_term = _exp(self.growth*tau)*_norm(s['initial_residual'])
        cheap = initial_term+tau*_exp(max(self.growth, projected_growth)*tau)*_norm(e)*_norm(c)
        # L*L = E*E in exact arithmetic; retain any negative roundoff spectrum
        # explicitly rather than presenting this square root as a certificate.
        eg = e.conj().T @ e
        values, vectors = np.linalg.eigh(eg)
        l = np.sqrt(np.maximum(values, 0.))[:, None]*vectors.conj().T
        bp = np.linalg.matrix_power(b, bound_degree)
        remainder_coefficient = _norm(l @ bp)
        local = initial_term
        half = tau/(2*bound_cells)
        for cell in range(bound_cells):
            left, midpoint = 2*cell*half, (2*cell+1)*half
            y = expm(midpoint*b) @ c
            power_y = y
            cell_integral = 0.
            for degree in range(bound_degree):
                cell_integral += 2*half**(degree+1)/factorial(degree+1)*_norm(l @ power_y)
                power_y = b @ power_y
            cell_integral += (2*half**(bound_degree+1)/factorial(bound_degree+1)
                              *remainder_coefficient*_norm(y)*_exp(two_sided_growth*half))
            local += _exp(self.growth*(tau-left))*cell_integral
        bound = min(cheap, local)
        return output, dict(time=float(time), absolute_error_envelope_estimate=float(bound),
                            cheap_envelope_estimate=float(cheap),
                            cell_taylor_envelope_estimate=float(local),
                            initial_error_term=float(initial_term),
                            ambient_growth_bound=self.growth,
                            projected_forward_growth=projected_growth,
                            projected_two_sided_growth=two_sided_growth,
                            residual_gram_minimum_eigenvalue=float(values[0]),
                            bound_cells=bound_cells, bound_degree=bound_degree,
                            seconds=perf_counter()-start,
                            qualification='exact-arithmetic bound formula, floating-point estimate; NOT an interval certificate')

    def expand(self):
        """Expand latest full-row H block; retain discarded numerical directions.

        Returns (receipt, discarded_residual). A caller can publish this before
        issuing another expansion. Basis/action limits stop rather than relax.
        """
        if self.blocks >= self.max_blocks or self.size >= self.max_basis:
            return dict(status='basis_limit', basis=self.size), None
        start = perf_counter()
        q = self.q
        w = self.hq[:, self.last_start:self.size].copy()
        incoming_scale = _norm(w)
        for _ in range(2):
            w -= q @ (q.conj().T @ w)
        u, singular, _ = np.linalg.svd(w, full_matrices=False)
        threshold = self.deflation_rtol*incoming_scale
        rank = int(np.count_nonzero(singular > threshold))
        available = min(self.max_basis-self.size, self.action_cap-self.requested_columns)
        accepted = min(rank, available)
        new = u[:, :accepted]
        # Account all omitted expansion directions, numerical AND budget-limited.
        discarded = w-new @ (new.conj().T @ w)
        receipt = dict(status='expanded' if accepted else ('numerical_closure' if not rank else 'action_limit'),
                       previous_basis=self.size, candidate_columns=w.shape[1],
                       accepted_columns=accepted, numerical_rank=rank,
                       threshold=threshold, incoming_action_frobenius=incoming_scale,
                       singular_values=singular.tolist(),
                       discarded_frobenius=_norm(discarded),
                       discarded_column_norms=np.linalg.norm(discarded, axis=0).tolist())
        if accepted:
            applied = self._apply(new)
            old = self.size
            self._q[:, old:old+accepted] = new
            self._hq[:, old:old+accepted] = applied
            self.last_start, self.size, self.blocks = old, old+accepted, self.blocks+1
            self._stage = None
        receipt.update(seconds=perf_counter()-start, basis=self.size,
                       requested_h_vector_columns=self.requested_columns)
        self.expansions.append(receipt)
        return receipt, discarded
