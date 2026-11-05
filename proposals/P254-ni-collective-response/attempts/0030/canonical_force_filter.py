"""Conditional active-model physical-force -> terminal displacement pullback.

The observable uses the bilinear pairing z.T @ state, never a Hermitian
pairing. A physical ionic force is added to p-prime, and a field-momentum
force to Pi-prime, both with plus sign. Thus their kernels are exactly the
p and Pi entries of z(tau) propagated with z'=-L(tau).T z backwards.

Only these small forcing entries of each accepted DOP853 dense interpolant
are retained. The full dual survives ONLY at preparation time zero, for
later preparation/source cross terms. This is not passive noise or loss.
"""
from __future__ import annotations

from dataclasses import dataclass

from pathlib import Path
import sys
import numpy as np
from scipy.integrate import DOP853

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / '0029'))
from reduced_ode_initial_step import bounded_first_step


_DEGREE = 7
# Installed SciPy rk.py:Dop853DenseOutput evaluates seven nested factors
# (dop853_coefficients.INTERPOLATOR_POWER=7): a degree <=7 polynomial.
# Eight Lobatto samples therefore change representation, not the physical
# approximation or its adaptive error tolerance. No private solver data used.
def chebyshev_lobatto_transform(degree):
    """Degree-preserving transform from descending Chebyshev Lobatto nodes."""
    nodes = np.cos(np.pi * np.arange(degree + 1) / degree)
    transform = (2.0 / degree) * np.cos(
        np.pi * np.outer(np.arange(degree + 1), np.arange(degree + 1)) / degree)
    transform[:, [0, -1]] *= 0.5
    transform[[0, -1], :] *= 0.5
    return nodes, transform


_NODES, _TRANSFORM = chebyshev_lobatto_transform(_DEGREE)


@dataclass(frozen=True)
class CanonicalForceFilter:
    """Small, ascending-time piecewise-Chebyshev force kernels.

    segment_bounds[s] = [left,right], and local x=(2*t-left-right)/(right-left).
    ionic_chebyshev[s,k,i,p] multiplies T_k(x) for ionic force coordinate i
    and terminal probe p; field_chebyshev has the same layout with field i.
    boundary_times and *_boundary contain every actual accepted endpoint,
    including zero and final_time. Exact-knot evaluation uses those values;
    interior evaluation uses the polynomial on that interval. Empty field
    axes for q != 0 are intentional: global field modes belong only to Gamma.
    """

    final_time: float
    initial_dual: np.ndarray
    segment_bounds: np.ndarray
    ionic_chebyshev: np.ndarray
    field_chebyshev: np.ndarray
    boundary_times: np.ndarray
    ionic_boundary: np.ndarray
    field_boundary: np.ndarray
    metadata: dict

    @property
    def arrays(self):
        """Numeric-only, NPZ-ready persistence payload (metadata is separate JSON)."""
        return {
            'initial_dual': self.initial_dual,
            'segment_bounds': self.segment_bounds,
            'ionic_chebyshev': self.ionic_chebyshev,
            'field_chebyshev': self.field_chebyshev,
            'boundary_times': self.boundary_times,
            'ionic_boundary': self.ionic_boundary,
            'field_boundary': self.field_boundary,
        }

    def evaluate(self, time):
        """Return (ionic_force[ndof,nprobes], field_force[nfield,nprobes]).

        These are response coefficients, not realizations of a stochastic
        force. No extrapolation, conjugation, mass or cell normalization.
        """
        if np.ndim(time) != 0:
            raise ValueError('one scalar force-insertion time required')
        time = float(time)
        if not np.isfinite(time) or not 0.0 <= time <= self.final_time:
            raise ValueError('force-insertion time must lie in [0, final_time]')
        right = int(np.searchsorted(self.boundary_times, time, side='left'))
        if right < len(self.boundary_times) and self.boundary_times[right] == time:
            return self.ionic_boundary[right].copy(), self.field_boundary[right].copy()
        segment = right - 1
        left_time, right_time = self.segment_bounds[segment]
        x = 2.0 * ((time - left_time) / (right_time - left_time)) - 1.0
        return (np.polynomial.chebyshev.chebval(x, self.ionic_chebyshev[segment]),
                np.polynomial.chebyshev.chebval(x, self.field_chebyshev[segment]))

    __call__ = evaluate


def build_canonical_force_filter(correlations, q, final_time, coordinates, *, capture=None):
    """Integrate the existing moving L_q(t).T using its original rtol/atol.

    The one q coefficient block is held just as in backward_linear_probes.
    DOP853.step(), .t_old/.t and .dense_output() are public SciPy interfaces;
    no OdeSolution/full-state history or transition matrix is retained.
    ``capture`` optionally receives ephemeral duals through boundary(time, dual)
    and segment(left, right, dense, shape). It must retain only its small
    physical insertions, never the full-state arrays. Existing force consumers
    need no callback and retain their unchanged array/evaluation contract.
    """
    if np.ndim(final_time) != 0:
        raise ValueError('scalar terminal time in inverse Ry required')
    final_time = float(final_time)
    if not np.isfinite(final_time) or final_time < 0.0:
        raise ValueError('finite nonnegative terminal time in inverse Ry required')
    if not isinstance(q, (int, np.integer)):
        raise ValueError('integer q index required')
    q = int(q)
    if not 0 <= q < correlations.c.nq:
        raise ValueError('q index outside coefficient mesh')
    coordinates = tuple(coordinates)
    terminal = correlations.ionic_probes(q, coordinates)
    shape = terminal.shape
    d = int(correlations.c.d)
    f = int(correlations.c.f) if q == 0 else 0
    p = slice(d, 2*d)
    pi = slice(2*d + f, 2*d + 2*f)

    def force_entries(value):
        dual = value.reshape(shape)
        return dual[p].copy(), dual[pi].copy()

    bounds, ionic_segments, field_segments = [], [], []
    boundary_times = [final_time]
    ionic_end, field_end = force_entries(terminal)
    ionic_boundary, field_boundary = [ionic_end], [field_end]
    if capture is not None:
        capture.boundary(final_time, terminal)
    nfev = 0
    first_step = None
    if final_time == 0.0:
        initial_dual = terminal.copy()
    else:
        block = correlations.c.q(q)

        def rhs(time, value):
            return -correlations._transpose(q, time, value.reshape(shape), block).ravel()

        interval = (final_time, 0.0)
        first_step = bounded_first_step(rhs, interval, terminal.ravel(),
                                        rtol=correlations.rtol, atol=correlations.atol)
        solver = DOP853(rhs, final_time, terminal.ravel(), 0.0,
                        rtol=correlations.rtol, atol=correlations.atol,
                        first_step=first_step)
        # terminal is no longer needed; DOP853 owns its initial state reference.
        while solver.status == 'running':
            message = solver.step()
            if solver.status == 'failed':
                raise RuntimeError(f'backward canonical force filter solve failed: {message}')
            left, right = float(solver.t), float(solver.t_old)
            if not left < right:
                raise RuntimeError('backward DOP853 accepted a nondecreasing interval')
            dense = solver.dense_output()
            ionic_samples = np.empty((_DEGREE+1, d, shape[1]), complex)
            field_samples = np.empty((_DEGREE+1, f, shape[1]), complex)
            for j, node in enumerate(_NODES):
                # Explicit endpoints avoid roundoff excursions outside a bounded
                # baseline. Dense evaluation itself never calls the baseline.
                time = right if j == 0 else left if j == _DEGREE else (
                    left + (right-left)*(node+1.0)/2.0)
                value = dense(time)
                ionic_samples[j], field_samples[j] = force_entries(value)
            ionic_coeff = np.tensordot(_TRANSFORM, ionic_samples, axes=(1, 0))
            field_coeff = np.tensordot(_TRANSFORM, field_samples, axes=(1, 0))
            if not (np.all(np.isfinite(ionic_coeff)) and np.all(np.isfinite(field_coeff))):
                raise FloatingPointError('nonfinite force polynomial; no clipping applied')
            bounds.append((left, right))
            ionic_segments.append(ionic_coeff)
            field_segments.append(field_coeff)
            boundary_times.append(left)
            ionic_end, field_end = force_entries(solver.y)
            ionic_boundary.append(ionic_end)
            field_boundary.append(field_end)
            if capture is not None:
                capture.segment(left, right, dense, shape)
                capture.boundary(left, solver.y.reshape(shape))
            # Do not retain a full-state interpolant between accepted steps.
            del dense, value
        initial_dual = solver.y.reshape(shape).copy()
        nfev = int(solver.nfev)
        del solver
    if not np.all(np.isfinite(initial_dual)):
        raise FloatingPointError('nonfinite initial dual; no clipping applied')
    count = len(bounds)
    segment_bounds = np.asarray(bounds[::-1], dtype=float).reshape(count, 2)
    ionic_chebyshev = (np.stack(ionic_segments[::-1]) if count else
                       np.empty((0, _DEGREE+1, d, shape[1]), complex))
    field_chebyshev = (np.stack(field_segments[::-1]) if count else
                       np.empty((0, _DEGREE+1, f, shape[1]), complex))
    metadata = {
        'schema': 'P254-canonical-force-filter-v1',
        'q_index': q,
        'coordinates': [int(i) for i in coordinates],
        'final_time_inverse_ry': final_time,
        'state_size': int(shape[0]),
        'ndof': d,
        'nfield': f,
        'nprobes': int(shape[1]),
        'method': 'DOP853',
        'rtol': np.asarray(correlations.rtol).tolist(),
        'atol': np.asarray(correlations.atol).tolist(),
        'first_step': first_step,
        'nfev': nfev,
        'accepted_intervals': count,
        'polynomial_degree': _DEGREE,
        'polynomial_basis': 'Chebyshev T_k, k=0..7; x=2*(t-left)/(right-left)-1',
        'polynomial_origin': 'eight public DOP853 dense-evaluator Lobatto samples per accepted interval',
        'approximation': 'representation of existing degree-seven DOP853 interpolant, not a new tolerance gate',
        'interval_order': 'ascending physical time, each row [left,right]',
        'boundary_policy': 'exact stored accepted-state values at every knot, including time zero and terminal time',
        'array_axes': {
            'initial_dual': ['state', 'probe'],
            'segment_bounds': ['segment', 'left_right'],
            'ionic_chebyshev': ['segment', 'degree', 'ionic_force', 'probe'],
            'field_chebyshev': ['segment', 'degree', 'field_momentum_force', 'probe'],
            'boundary_times': ['boundary'],
            'ionic_boundary': ['boundary', 'ionic_force', 'probe'],
            'field_boundary': ['boundary', 'field_momentum_force', 'probe'],
        },
        'state_layout': '[u(d), p(d), A(f), Pi(f), sigma(ns,m,m).ravel()] at Gamma; [u,p,sigma] otherwise',
        'ionic_force_state_slice': [d, 2*d],
        'field_force_state_slice': [2*d+f, 2*d+2*f],
        'pairing': 'nonconjugating z.T @ state; same q, no -q swap',
        'force_insertion': 'plus ionic force into p-prime, plus field-momentum force into Pi-prime',
        'extra_mass_weight_cell_factors': False,
        'history': 'only small canonical-force polynomials and boundaries; full dual only at preparation zero',
        'scope': 'conditional moving active-model causal output filter; not full passive noise, physical loss or lifetime',
    }
    return CanonicalForceFilter(
        final_time, initial_dual, segment_bounds, ionic_chebyshev, field_chebyshev,
        np.asarray(boundary_times[::-1], dtype=float),
        np.stack(ionic_boundary[::-1]), np.stack(field_boundary[::-1]), metadata,
    )
