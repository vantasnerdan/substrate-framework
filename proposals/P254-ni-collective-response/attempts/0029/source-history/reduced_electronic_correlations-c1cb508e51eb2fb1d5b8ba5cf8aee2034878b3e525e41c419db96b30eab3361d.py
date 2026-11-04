"""Connected Gaussian-linearized observables of the approved active-band model.

This consumer depends only on the numeric ReducedCoefficients interface and a
baseline(time)->(rho_b, h_b) callback.  It does not import a dynamics solver.
It retains the nonstationary preparation, explicit spins, q/source convention,
original product ionic/field covariance and ordered fermionic Wick noise.
Missing passive/complement electronic response AND noise are not supplied.
These correlations are not full-PW noise, entanglement or physical lifetimes.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

from reduced_correlation_transpose import (
    commutator_tangent, commutator_transpose, wick_edge_pair,
)
from reduced_ode_initial_step import bounded_first_step

Baseline = Callable[[float], tuple[np.ndarray, np.ndarray]]
NOISE_SCOPE = (
    'active-band Gaussian-linearized connected noise with original product '
    'ionic/field preparation; passive/complement electronic noise and '
    'active-passive initial coherences are not represented'
)


@dataclass(frozen=True)
class OrderedComponents:
    """AB=<A(t)B(s)>c and BA=<B(s)A(t)>c in the SAME output axis order."""

    ordered_ab: np.ndarray
    ordered_ba: np.ndarray
    symmetrized: np.ndarray
    commutator: np.ndarray
    noise_scope: str = NOISE_SCOPE

    @classmethod
    def from_ordered(cls, ab, ba):
        return cls(ab, ba, (ab + ba) / 2, ab - ba)


class ActiveBandCorrelations:
    """Matrix-free q tangents and backward nonconjugating dual propagation.

    State layout is [x_q, sigma_q.ravel()], with sigma[slot,a,b] the GLOBAL
    canonical density fluctuation (k+q,a;k,b), without a sector/spin weight.
    A column of dual coefficients z represents sum(z * state), not vdot.
    Only selected probe columns are propagated, never a full transition bank.
    """

    def __init__(self, coefficients, baseline: Baseline, *, rtol=2e-8,
                 atol=2e-10, method='DOP853'):
        self.c = coefficients
        self.baseline = baseline
        self.rtol, self.atol, self.method = rtol, atol, method
        self.root_n = np.sqrt(coefficients.n)

    def state_size(self, q):
        return len(self.c.mechanical.generators[q]) + self.c.ns * self.c.m**2

    def _split(self, q, value):
        value = np.asarray(value, dtype=complex)
        if value.ndim not in (1, 2) or value.shape[0] != self.state_size(q):
            raise ValueError('state/dual must have shape (q state size[, probes])')
        matrix = value[:, None] if value.ndim == 1 else value
        nx = len(self.c.mechanical.generators[q])
        return matrix[:nx], matrix[nx:].reshape(self.c.ns, self.c.m, self.c.m, -1), value.ndim

    @staticmethod
    def _join(x, sigma, ndim):
        result = np.concatenate((x, sigma.reshape(-1, x.shape[1])), axis=0)
        return result[:, 0] if ndim == 1 else result

    def _baseline(self, time):
        rho, h = self.baseline(float(time))
        shape = (self.c.ns, self.c.m, self.c.m)
        if np.shape(rho) != shape or np.shape(h) != shape:
            raise ValueError(f'baseline must supply rho,h with shape {shape}')
        return rho, h

    def _tangent(self, q, time, value, block):
        c = self.c
        x, sigma, ndim = self._split(q, value)
        rho, h = self._baseline(time)
        target = c.target[q]
        density = np.einsum('srab,sabp->rp', block.density, sigma,
                            optimize=True) / self.root_n
        potential = block.kernel @ density
        dh = np.einsum('siab,ip->sabp', block.ionic, x[:c.d],
                       optimize=True) / self.root_n
        dh += np.einsum('sabr,rp->sabp', block.potential, potential,
                        optimize=True) / self.root_n
        gradient = np.zeros_like(x)
        gradient[:c.d] = np.einsum('siab,sabp->ip', block.ionic_reverse,
                                    sigma, optimize=True) / self.root_n
        if q == 0:
            a = slice(2*c.d, 2*c.d+c.f)
            dh += np.einsum('siab,ip->sabp', c.field_modes, x[a], optimize=True)
            gradient[a] = np.einsum('siab,sabp->ip', c.field_reverse_modes,
                                     sigma, optimize=True)
        dx = c.mechanical.generators[q] @ x + c.mechanical.j[q] @ gradient
        dsigma = commutator_tangent(h[target], h, rho[target], rho, sigma, dh)
        return self._join(dx, dsigma, ndim)

    def tangent(self, q, time, state):
        """Forward L_q(time) state, including screening and electronic feedback."""
        return self._tangent(q, time, state, self.c.q(q))

    def _transpose(self, q, time, value, block):
        c = self.c
        zx, zs, ndim = self._split(q, value)
        rho, h = self._baseline(time)
        target = c.target[q]
        out_s, dh_dual = commutator_transpose(h[target], h, rho[target], rho, zs)
        out_x = c.mechanical.generators[q].T @ zx
        out_x[:c.d] += np.einsum('siab,sabp->ip', block.ionic, dh_dual,
                                  optimize=True) / self.root_n
        feature_dual = np.einsum('sabr,sabp->rp', block.potential, dh_dual,
                                 optimize=True) / self.root_n
        density_dual = block.kernel.T @ feature_dual
        out_s += np.einsum('srab,rp->sabp', block.density, density_dual,
                           optimize=True) / self.root_n
        gradient_dual = c.mechanical.j[q].T @ zx
        out_s += np.einsum('siab,ip->sabp', block.ionic_reverse,
                           gradient_dual[:c.d], optimize=True) / self.root_n
        if q == 0:
            a = slice(2*c.d, 2*c.d+c.f)
            out_x[a] += np.einsum('siab,sabp->ip', c.field_modes, dh_dual,
                                  optimize=True)
            out_s += np.einsum('siab,ip->sabp', c.field_reverse_modes,
                               gradient_dual[a], optimize=True)
        return self._join(out_x, out_s, ndim)

    def transpose(self, q, time, dual):
        """True L_q(time).T dual: deliberately no conjugation and no -q swap."""
        return self._transpose(q, time, dual, self.c.q(q))

    def ionic_probes(self, q, coordinates):
        coordinates = tuple(coordinates)
        if not coordinates or any(not isinstance(i, (int, np.integer)) or
                                  not 0 <= i < self.c.d for i in coordinates):
            raise ValueError('nonempty original Cartesian displacement coordinates required')
        probes = np.zeros((self.state_size(q), len(coordinates)), complex)
        probes[np.asarray(coordinates), np.arange(len(coordinates))] = 1
        return probes

    def backward_linear_probes(self, q, time, probes):
        """Pull arbitrary linear terminal probes at time back to preparation 0.

        z'(tau)=-L(tau).T z(tau), integrated from terminal time down to zero.
        Baseline is evaluated at the actual intermediate time, never frozen.
        A single lazy q block is held during this integration.
        """
        time = float(time)
        if not np.isfinite(time) or time < 0:
            raise ValueError('finite nonnegative time in inverse Ry required')
        probes = np.asarray(probes, dtype=complex)
        self._split(q, probes)
        if time == 0:
            return probes.copy()
        block = self.c.q(q)
        shape = probes.shape
        rhs = lambda tau, z: -self._transpose(q, tau, z.reshape(shape), block).ravel()
        interval, initial = (time, 0.0), probes.ravel()
        first = bounded_first_step(rhs, interval, initial, rtol=self.rtol, atol=self.atol)
        result = solve_ivp(rhs, interval, initial, method=self.method,
                           rtol=self.rtol, atol=self.atol, t_eval=[0.0], first_step=first)
        if not result.success:
            raise RuntimeError(f'backward active correlation solve failed: {result.message}')
        dual = result.y[:, -1].reshape(shape)
        if not np.all(np.isfinite(dual)):
            raise FloatingPointError('nonfinite active correlation dual; no clipping applied')
        return dual

    def backward_probes(self, q, time, ionic_coordinates):
        return self.backward_linear_probes(q, time, self.ionic_probes(q, ionic_coordinates))

    def initial_ordered(self, q, left, right):
        """<left_q right_-q>c in the exact specified PRODUCT preparation.

        right's source slot is target(q,left_source), NOT left_source.
        The electronic Wick sum has NO 1/N or extra spin factor.  The only
        1/N for spatial correlation occurs in the final Fourier readout.
        """
        minus = int(self.c.minus[q])
        lx, ls, _ = self._split(q, left)
        rx, rs, _ = self._split(minus, right)
        mechanical = lx.T @ self.c.mechanical.w0[q] @ rx
        rho = self.c.rho0
        target = self.c.target[q]
        electronic = wick_edge_pair(ls, rs[target], rho, rho[target])
        return mechanical + electronic

    def two_time_q(self, q, time_left, time_right, left_coordinates, right_coordinates=None):
        """Ordered/symmetric/commutator q,-q correlations of ionic displacements."""
        if right_coordinates is None:
            right_coordinates = left_coordinates
        minus = int(self.c.minus[q])
        left = self.backward_probes(q, time_left, left_coordinates)
        right = self.backward_probes(minus, time_right, right_coordinates)
        ab = self.initial_ordered(q, left, right)
        ba = self.initial_ordered(minus, right, left).T
        return OrderedComponents.from_ordered(ab, ba)

    def equal_time_q(self, q, time, coordinates):
        return self.two_time_q(q, time, time, coordinates)

    def spatial(self, time_left, time_right, displacements, left_coordinates,
                right_coordinates=None):
        """Real-space connected correlation for R_left-R_right (integer cells).

        Orthonormal Fourier convention: (1/N) sum_q exp(i q.R) C_q.
        No reality projection is performed: imaginary residuals remain visible.
        Each q is consumed independently; no giant transition/covariance bank.
        """
        cells = np.asarray(displacements)
        if cells.ndim == 1:
            cells = cells[None, :]
        if cells.ndim != 2 or cells.shape[1] != 3 or not np.all(np.isfinite(cells)):
            raise ValueError('finite cell displacements with shape (count,3) required')
        if not np.all(cells == np.rint(cells)):
            raise ValueError('spatial displacements must be integer cells')
        left_coordinates = tuple(left_coordinates)
        right_coordinates = left_coordinates if right_coordinates is None else tuple(right_coordinates)
        shape = (len(cells), len(left_coordinates), len(right_coordinates))
        ab, ba = np.zeros(shape, complex), np.zeros(shape, complex)
        for q in range(self.c.nq):
            components = self.two_time_q(q, time_left, time_right,
                                          left_coordinates, right_coordinates)
            phase = np.exp(2j*np.pi * ((cells / np.asarray(self.c.mesh)) @ self.c.q_indices[q]))
            ab += phase[:, None, None] * components.ordered_ab / self.c.n
            ba += phase[:, None, None] * components.ordered_ba / self.c.n
        return OrderedComponents.from_ordered(ab, ba)
