"""Original-basis trace dual of the finite Gamma linearized coupled action.

Pairing is O=u_probe.u+p_probe.p+sum_c Tr(A_c delta_rho_c), bilinear,
not sesquilinear and not elementwise dot. Forward action is the existing
+g*u, density K, and NEGATIVE ionic gradient in Ry time (hbar=1).
No exact interacting evolution, stationary FDT or general-q transfer is claimed.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class BilinearProbe:
    """Coefficients of an operator; arbitrary complex linear combinations allowed.

u/p coefficients carry observable/bohr and observable/(Ry-time momentum)
units, electronic matrices carry observable units. density/force helpers below
include the physical 1/N weight; raw constructor coefficients include no weight.
"""
    u: np.ndarray
    p: np.ndarray
    electronic: tuple[np.ndarray, ...]

    def adjoint(self):
        return BilinearProbe(self.u.conj(), self.p.conj(),
                             tuple(a.conj().T for a in self.electronic))

    def evaluate(self, u, p, matrices):
        return self.u @ u + self.p @ p + sum(
            np.trace(a @ r) for a, r in zip(self.electronic, matrices))


class PhysicalObservableTransfer:
    """Selected matrix probes: storage O(sum d_c^2), no all-observable covariance.

model is a validated InertialFermionDynamics (channels, interaction,
frozen_curvature, masses, nion, ndensity); its tangent/energy/preparation APIs
are NOT used. Only its original-basis action inputs are shared.
"""
    def __init__(self, model, modes):
        self.model, self.modes = model, modes
        modes.validate(model.channels)
        self._shapes = tuple(c.h0.shape for c in model.channels)
        self._sizes = tuple(c.h0.size for c in model.channels)

    def zero_probe(self):
        zero = np.zeros(self.model.nion, complex)
        return BilinearProbe(zero.copy(), zero.copy(),
                             tuple(np.zeros(s, complex) for s in self._shapes))

    def density_probe(self, coefficients):
        coefficients = np.asarray(coefficients, complex)
        if coefficients.shape != (self.model.ndensity,):
            raise ValueError("one density coefficient per density slot required")
        zero = self.zero_probe()
        return BilinearProbe(zero.u, zero.p, tuple(c.weight*np.einsum(
            "s,sij->ij", coefficients, c.densities) for c in self.model.channels))

    def force_probe(self, direction):
        """TOTAL linear ionic force -Cfrozen*u-sum_c w_c Tr(g_c delta_rho)."""
        direction = np.asarray(direction, complex)
        if direction.shape != (self.model.nion,):
            raise ValueError("one direction coefficient per ionic coordinate required")
        zero = self.zero_probe()
        return BilinearProbe(-self.model.frozen_curvature.T @ direction, zero.p,
                             tuple(-c.weight*np.einsum("a,aij->ij", direction, c.g)
                                   for c in self.model.channels))

    def dual(self, probe):
        """Derive dO/dt by cyclic trace, not by transposing a packed generator.

Tr(A*(-i[V,rho0]))=-i Tr(A[V,rho0]); density feedback then passes
through K.T and w*density. The free dual is i[H,A]. The force term
contributes -w*p_probe.g. No conjugation occurs in this bilinear dual.
"""
        self._flatten(probe)
        m = self.model
        gu = np.zeros(m.nion, complex)
        dn = np.zeros(m.ndensity, complex)
        for c, a in zip(m.channels, probe.electronic):
            rho = c.rho0
            gu += -1j*np.einsum("ij,aji->a", a, c.g @ rho-rho @ c.g)
            dn += -1j*np.einsum("ij,sji->s", a, c.densities @ rho-rho @ c.densities)
        feedback = m.interaction.T @ dn
        electronic = tuple(1j*(c.h0 @ a-a @ c.h0)
            - c.weight*np.einsum("a,aij->ij", probe.p, c.g)
            + c.weight*np.einsum("s,sij->ij", feedback, c.densities)
            for c, a in zip(m.channels, probe.electronic))
        return BilinearProbe(gu-m.frozen_curvature.T @ probe.p,
                             probe.u/m.masses, electronic)

    def _flatten(self, probe):
        if (np.shape(probe.u) != (self.model.nion,) or np.shape(probe.p) != (self.model.nion,)
                or len(probe.electronic) != len(self._shapes)
                or any(np.shape(a) != s for a, s in zip(probe.electronic, self._shapes))):
            raise ValueError("probe shape differs from original-basis model")
        value = np.concatenate((probe.u, probe.p, *(a.ravel() for a in probe.electronic)))
        if not np.all(np.isfinite(value)):
            raise ValueError("finite probe required")
        return value.astype(complex)

    def _unflatten(self, value):
        n = self.model.nion
        cursor, matrices = 2*n, []
        for size, shape in zip(self._sizes, self._shapes):
            matrices.append(value[cursor:cursor+size].reshape(shape))
            cursor += size
        return BilinearProbe(value[:n], value[n:2*n], tuple(matrices))

    def pullback(self, probe, time, *, rtol=1e-10, atol=1e-12):
        """O(time)=pulled coefficients contracted with INITIAL operators."""
        initial = self._flatten(probe)
        if not np.isfinite(time):
            raise ValueError("finite Ry time required")
        if time == 0:
            return self._unflatten(initial.copy())
        result = solve_ivp(lambda t, value: self._flatten(self.dual(self._unflatten(value))),
                           (0, time), initial, method="DOP853", t_eval=[time],
                           rtol=rtol, atol=atol)
        if not result.success:
            raise RuntimeError(result.message)
        return self._unflatten(result.y[:, -1])

    def two_time(self, preparation, left, left_time, right, right_time, **solver_options):
        """Prepared C(t,s), NOT stationary C(t-s); no equilibrium inference."""
        if (preparation.modes != self.modes or len(preparation.channels) != len(self.model.channels)
                or any(a is not b for a, b in zip(preparation.channels, self.model.channels))):
            raise ValueError("preparation and generator must share exact channel identities")
        result = preparation.correlations(self.pullback(left, left_time, **solver_options),
                                           self.pullback(right, right_time, **solver_options))
        result.update(left_time=float(left_time), right_time=float(right_time),
                      time_units="Ry^-1 with hbar=1", stationary_assumed=False,
                      scope="linearized mean-field prepared two-point transfer at Gamma")
        return result
