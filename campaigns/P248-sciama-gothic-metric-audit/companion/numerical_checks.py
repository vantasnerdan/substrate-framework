"""Numerical applicability and independent wave-source regressions for P248."""

from __future__ import annotations

import math

import numpy as np
import sympy as sp
from scipy import constants
from scipy.integrate import solve_ivp

from substrate_framework.optical_gothic import optical_adm_metric, recover_optical_adm
from substrate_framework.verification import CheckLedger


def _inside_printed_rounding(value: float, printed: float, last_place: float) -> bool:
    return abs(value - printed) <= 0.5 * last_place


def _wave_solution(points: int, final_time: float, source_sign: float = 1.0) -> tuple[float, float]:
    """Second-order periodic leapfrog for u_tt=u_xx+sign*sin(x)*cos(2t)."""

    domain = 2.0 * math.pi
    dx = domain / points
    spatial = np.arange(points, dtype=float) * dx
    steps = math.ceil(final_time / (0.2 * dx))
    dt = final_time / steps

    def laplacian(field: np.ndarray) -> np.ndarray:
        return (np.roll(field, -1) - 2.0 * field + np.roll(field, 1)) / dx**2

    previous = np.zeros(points)
    acceleration0 = source_sign * np.sin(spatial)
    current = previous + 0.5 * dt**2 * acceleration0
    for step in range(1, steps):
        time = step * dt
        acceleration = laplacian(current) + source_sign * np.sin(spatial) * np.cos(2.0 * time)
        following = 2.0 * current - previous + dt**2 * acceleration
        previous, current = current, following

    exact_amplitude = (math.cos(final_time) - math.cos(2.0 * final_time)) / 3.0
    exact = exact_amplitude * np.sin(spatial)
    error = float(np.sqrt(np.mean((current - exact) ** 2)))
    scale = float(np.sqrt(np.mean(exact**2)))
    return error, scale


def run() -> int:
    ledger = CheckLedger("P248-NUMERIC")

    # Independent random SPD round trips use fixed seed and direct NumPy
    # inversion, not the symbolic inverse returned by the package API.
    rng = np.random.default_rng(248)
    for sample in range(8):
        factor = rng.normal(size=(3, 3))
        spatial_np = factor.T @ factor + np.eye(3)
        flow_np = rng.normal(size=3)
        lapse = 0.5 + rng.random()
        speed = 2.0 + rng.random()
        # Floating samples are converted through decimal strings to exact
        # rationals before they cross the canonical symbolic API boundary.
        spatial_exact = [
            [sp.Rational(str(float(value))) for value in row]
            for row in spatial_np
        ]
        flow_exact = [sp.Rational(str(float(value))) for value in flow_np]
        lapse_exact = sp.Rational(str(lapse))
        speed_exact = sp.Rational(str(speed))
        metric = optical_adm_metric(
            lapse_exact, spatial_exact, flow_exact, speed_exact
        )
        covariant = np.array(metric.covariant, dtype=float)
        contravariant = np.array(metric.contravariant, dtype=float)
        residual = np.linalg.norm(covariant @ contravariant - np.eye(4), ord=np.inf)
        recovered = recover_optical_adm(metric.covariant, speed_exact)
        recovery_residual = np.linalg.norm(
            np.array(recovered.reconstructed_metric.covariant, dtype=float) - covariant,
            ord=np.inf,
        )
        ledger.check(
            f"SPD round trip {sample + 1}",
            residual < 2.0e-13 and recovery_residual < 2.0e-13,
            detail=f"inverse={residual:.3e}, recovery={recovery_residual:.3e}",
        )

    c = constants.c
    newton = constants.G
    hbar = constants.hbar
    planck_length = math.sqrt(hbar * newton / c**3)
    planck_density = c**5 / (hbar * newton**2)
    density = c**2 / (16.0 * math.pi * newton * planck_length**2)
    modulus = density * c**2
    impedance = density * c
    ledger.check(
        "printed Planck length rounding interval",
        _inside_printed_rounding(planck_length, 1.616255e-35, 1.0e-41),
        detail=f"{planck_length:.12e}",
    )
    ledger.check(
        "printed Planck density rounding interval",
        _inside_printed_rounding(planck_density, 5.155e96, 1.0e93),
        detail=f"{planck_density:.12e}",
    )
    ledger.check(
        "printed substrate density rounding interval",
        _inside_printed_rounding(density, 1.026e95, 1.0e92),
        detail=f"{density:.12e}",
    )
    ledger.check(
        "printed modulus rounding interval",
        _inside_printed_rounding(modulus, 9.22e111, 1.0e109),
        detail=f"{modulus:.12e}",
    )
    ledger.check(
        "printed impedance rounding interval",
        _inside_printed_rounding(impedance, 3.07e103, 1.0e101),
        detail=f"{impedance:.12e}",
    )
    ledger.check(
        "Planck-cutoff relation is algebraically circular",
        abs(density / (planck_density / (16.0 * math.pi)) - 1.0) < 2.0e-15,
    )

    final_time = 1.0
    errors = [_wave_solution(points, final_time)[0] for points in (64, 128, 256)]
    orders = [math.log(errors[index] / errors[index + 1], 2.0) for index in (0, 1)]
    ledger.check(
        "retarded-wave regression refines at second order",
        min(orders) > 1.8 and errors[-1] < 2.0e-5,
        detail=f"errors={errors}, orders={orders}",
    )
    wrong_error, exact_scale = _wave_solution(256, final_time, source_sign=-1.0)
    ledger.check(
        "wave-source sign mutation fails",
        wrong_error > exact_scale,
        detail=f"wrong_error={wrong_error:.3e}, scale={exact_scale:.3e}",
    )

    modal = solve_ivp(
        lambda time, state: [state[1], -state[0] + math.cos(2.0 * time)],
        (0.0, final_time),
        [0.0, 0.0],
        method="DOP853",
        rtol=1.0e-12,
        atol=1.0e-14,
    )
    exact_amplitude = (math.cos(final_time) - math.cos(2.0 * final_time)) / 3.0
    ledger.check("independent DOP853 solve succeeds", modal.success, modal.message)
    ledger.check(
        "independent modal amplitude agrees",
        abs(modal.y[0, -1] - exact_amplitude) < 2.0e-12,
        detail=f"error={abs(modal.y[0, -1] - exact_amplitude):.3e}",
    )
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(run())
