"""Same-frequency physical screening: (I - R0_nn(z) K_q) n = R0_ng(z) u.

The kernel is a callable on complex (nnr, 2) up/down density envelopes.
The supplied response fixes z for every action; no static-screened ionic vertex
is introduced. Storage, normalization and padding follow response.model.fft.
"""
from __future__ import annotations

from dataclasses import dataclass
from time import monotonic
from typing import Any, Callable

import numpy as np
from numpy.typing import NDArray
from scipy.sparse.linalg import LinearOperator, gmres


@dataclass(frozen=True)
class ScreenedResponseValue:
    density: NDArray[np.complex128]
    ionic_gradient: NDArray[np.complex128]
    induced_potential: NDArray[np.complex128]
    diagnostics: dict[str, Any]


class ScreenedResponseError(RuntimeError):
    """No completed column: reason and complete work counters are in diagnostics."""

    def __init__(self, reason: str, diagnostics: dict[str, Any]):
        self.reason = reason
        self.diagnostics = diagnostics
        super().__init__(f"screened response failed: {reason}")


def solve_screened_response(response, kernel, displacement, *, rtol: float,
                            max_matvec: int, wall_seconds: float,
                            progress: Callable[[dict[str, Any]], None] | None = None
                            ) -> ScreenedResponseValue:
    """Return one completed physical screened response, or raise explicitly.

    max_matvec is a HARD response-action budget, including the one bare action
    and the independent final residual/gradient action (minimum two). GMRES
    operator calls and inner iterations are recorded separately. Restart storage
    is at most 20 density vectors plus SciPy's fixed working arrays. Wall checks
    use monotonic time before/after each action and at iteration boundaries;
    an individual kernel/response action cannot be interrupted mid-call.

    The final residual is ||n - R0_nn K n - R0_ng u|| / ||R0_ng u||,
    independently evaluated with response.apply(K(n), u). For exactly zero RHS,
    the accepted solution is zero and its final residual must also be zero.
    Invalid arguments raise ValueError; execution failures raise
    ScreenedResponseError with diagnostics, never a partial result. progress
    receives diagnostic snapshots, not usable ionic columns.
    """
    if not np.isfinite(rtol) or rtol <= 0:
        raise ValueError("rtol must be finite and strictly positive")
    if isinstance(max_matvec, (bool, np.bool_)) or not isinstance(max_matvec, (int, np.integer)) or max_matvec < 2:
        raise ValueError("max_matvec must be an integer of at least two")
    if not np.isfinite(wall_seconds) or wall_seconds <= 0:
        raise ValueError("wall_seconds must be finite and strictly positive")
    u = np.asarray(displacement, dtype=np.complex128)
    if u.shape != (6,) or not np.all(np.isfinite(u)):
        raise ValueError("displacement must be a finite six-component vector")
    shape = (response.model.fft.nnr, 2)
    size = int(np.prod(shape))
    zero_u = np.zeros(6, dtype=np.complex128)
    started = monotonic()
    diagnostics: dict[str, Any] = {
        "status": "running", "z_ry": [float(response.z.real), float(response.z.imag)],
        "rtol": float(rtol), "max_matvec": int(max_matvec),
        "budget_scope": "all response.apply calls including bare and final",
        "wall_seconds": float(wall_seconds), "response_actions": 0,
        "response_seconds": 0.0, "kernel_actions": 0, "kernel_seconds": 0.0,
        "gmres_matvecs": 0, "iterations": 0, "completed_outputs": 0,
        "iteration_residuals": [], "elapsed_seconds": 0.0,
    }

    def snapshot():
        diagnostics["elapsed_seconds"] = monotonic() - started
        return {**diagnostics, "iteration_residuals": list(diagnostics["iteration_residuals"])}

    def fail(reason):
        diagnostics["status"] = "failed"
        diagnostics["reason"] = reason
        raise ScreenedResponseError(reason, snapshot())

    def check_wall():
        if monotonic() - started >= wall_seconds:
            fail("wall_limit")

    def notify():
        if progress is not None:
            progress(snapshot())
        check_wall()

    def checked_array(value, expected, name):
        array = np.asarray(value, dtype=np.complex128)
        if array.shape != expected or not np.all(np.isfinite(array)):
            diagnostics["invalid_output"] = name
            fail("invalid_output")
        return array

    def apply_response(potential, drive, *, reserve_final=False):
        check_wall()
        limit = max_matvec - int(reserve_final)
        if diagnostics["response_actions"] >= limit:
            fail("work_limit")
        diagnostics["response_actions"] += 1
        tick = monotonic()
        try:
            value = response.apply(potential, drive)
        finally:
            diagnostics["response_seconds"] += monotonic() - tick
        check_wall()
        density = checked_array(value.density, shape, "density")
        gradient = checked_array(value.ionic_gradient, (6,), "ionic_gradient")
        notify()
        return density, gradient

    def apply_kernel(density):
        check_wall()
        diagnostics["kernel_actions"] += 1
        tick = monotonic()
        try:
            potential = kernel(density)
        finally:
            diagnostics["kernel_seconds"] += monotonic() - tick
        check_wall()
        return checked_array(potential, shape, "induced_potential")

    def matvec(vector):
        check_wall()
        # Reserve the final independent residual/observable action before K work.
        if diagnostics["response_actions"] >= max_matvec - 1:
            fail("work_limit")
        diagnostics["gmres_matvecs"] += 1
        density = vector.reshape(shape)
        potential = apply_kernel(density)
        induced, _ = apply_response(potential, zero_u, reserve_final=True)
        return (density - induced).reshape(-1)

    def iteration(residual):
        diagnostics["iterations"] += 1
        residual = float(residual)
        if not np.isfinite(residual):
            fail("nonfinite_iteration_residual")
        diagnostics["iteration_residuals"].append(residual)
        notify()

    try:
        bare, _ = apply_response(np.zeros(shape, dtype=np.complex128), u,
                                 reserve_final=True)
        rhs_norm = float(np.linalg.norm(bare))
        if not np.isfinite(rhs_norm):
            fail("nonfinite_rhs_norm")
        diagnostics["bare_density_norm"] = rhs_norm
        if rhs_norm == 0:
            density = np.zeros(shape, dtype=np.complex128)
            info = 0
        else:
            operator = LinearOperator((size, size), matvec=matvec, dtype=np.complex128)
            solution, info = gmres(operator, bare.reshape(-1), rtol=rtol, atol=0.0,
                                   restart=min(20, size, max_matvec), maxiter=max_matvec,
                                   callback=iteration, callback_type="legacy")
            density = solution.reshape(shape)
        diagnostics["gmres_info"] = int(info)
        potential = apply_kernel(density)
        final_density, gradient = apply_response(potential, u)
        residual_norm = float(np.linalg.norm(density - final_density))
        diagnostics["final_residual_norm"] = residual_norm
        # Null RHS requires exact zero: there is no arbitrary unit-scale floor.
        relative = residual_norm / rhs_norm if rhs_norm else (0.0 if residual_norm == 0 else None)
        diagnostics["relative_residual"] = relative
        accepted = np.isfinite(residual_norm) and residual_norm <= rtol * rhs_norm
        if info != 0 or not accepted:
            fail("nonconvergence")
        check_wall()
        diagnostics["status"] = "complete"
        diagnostics["completed_outputs"] = 1
        return ScreenedResponseValue(density, gradient, potential, snapshot())
    except ScreenedResponseError:
        raise
    except Exception as error:
        diagnostics["cause"] = {"type": type(error).__name__, "message": str(error)}
        diagnostics["status"] = "failed"
        diagnostics["reason"] = "execution_failure"
        raise ScreenedResponseError("execution_failure", snapshot()) from error
