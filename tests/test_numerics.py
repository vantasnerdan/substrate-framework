from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import diags

from substrate_framework.numerics import (
    SolverTolerances,
    interpolating_spline_time_derivative,
    local_polynomial_time_derivative,
    refinement_study,
    solve_bvp_evidence,
    solve_ivp_evidence,
    solve_method_of_lines,
    trapezoid_integral,
    trapezoid_integral_axis,
)


TIGHT = SolverTolerances(rtol=1.0e-11, atol=1.0e-13)


def test_ivp_records_conserved_quantity_drift() -> None:
    result = solve_ivp_evidence(
        lambda _time, state: np.array([state[1], -state[0]]),
        (0.0, 2.0 * np.pi),
        [1.0, 0.0],
        sample_times=np.linspace(0.0, 2.0 * np.pi, 101),
        tolerances=TIGHT,
        invariant=lambda state: float(state @ state),
    )

    np.testing.assert_allclose(result.state[:, -1], [1.0, 0.0], atol=2.0e-10)
    assert result.max_abs_invariant_drift is not None
    assert result.max_abs_invariant_drift < 2.0e-10
    assert result.function_evaluations > 0


def test_bvp_records_collocation_residuals() -> None:
    coordinate = np.linspace(0.0, np.pi / 2.0, 21)
    guess = np.vstack((coordinate / coordinate[-1], np.ones_like(coordinate)))
    result = solve_bvp_evidence(
        lambda _x, state: np.vstack((state[1], -state[0])),
        lambda left, right: np.array([left[0], right[0] - 1.0]),
        coordinate,
        guess,
        tolerance=1.0e-8,
    )

    np.testing.assert_allclose(result.state[0], np.sin(result.coordinate), atol=2.0e-8)
    assert result.max_rms_residual < 1.0e-8
    assert result.iterations > 0
    assert result.parameters is None


def test_bvp_records_fitted_parameters() -> None:
    coordinate = np.linspace(0.0, 1.0, 11)
    guess = coordinate[None, :]
    result = solve_bvp_evidence(
        lambda _x, state, parameter: np.full_like(state, parameter[0]),
        lambda left, right, _parameter: np.array([left[0], right[0] - 1.0]),
        coordinate,
        guess,
        initial_parameters=[0.8],
        tolerance=1.0e-10,
    )

    np.testing.assert_allclose(result.state[0], result.coordinate, atol=1.0e-12)
    assert result.parameters is not None
    np.testing.assert_allclose(result.parameters, [1.0], atol=1.0e-12)


def test_method_of_lines_heat_equation_has_second_order_spatial_convergence() -> None:
    final_time = 0.02

    def solve(intervals: int) -> tuple[np.ndarray, np.ndarray]:
        spacing = 1.0 / intervals
        coordinate = np.linspace(0.0, 1.0, intervals + 1)[1:-1]
        count = coordinate.size
        laplacian = diags(
            [np.ones(count - 1), -2.0 * np.ones(count), np.ones(count - 1)],
            offsets=[-1, 0, 1],
            format="csr",
        ) / spacing**2
        result = solve_method_of_lines(
            lambda _time, state: laplacian @ state,
            (0.0, final_time),
            np.sin(np.pi * coordinate),
            sample_times=[final_time],
            tolerances=TIGHT,
        )
        return coordinate, result.state[:, -1]

    def rms_error(solution: tuple[np.ndarray, np.ndarray]) -> float:
        coordinate, numeric = solution
        exact = np.exp(-(np.pi**2) * final_time) * np.sin(np.pi * coordinate)
        return float(np.linalg.norm(numeric - exact) / np.sqrt(coordinate.size))

    study = refinement_study((16, 32, 64), solve, rms_error)

    assert study.errors_strictly_decrease
    assert study.final_error < 3.0e-5
    assert all(order is not None and order > 1.9 for order in study.observed_orders)


def test_trapezoid_integral_uses_current_numpy_api() -> None:
    coordinate = np.linspace(0.0, 1.0, 1001)
    assert trapezoid_integral(coordinate**2, coordinate) == pytest.approx(
        1.0 / 3.0, rel=2.0e-6
    )


def test_trapezoid_integral_axis_reuses_dispatch_and_validates_shape() -> None:
    coordinate = np.linspace(0.0, 1.0, 101)
    values = np.vstack((coordinate, coordinate**2))
    result = trapezoid_integral_axis(values, coordinate, axis=1)
    assert result == pytest.approx([0.5, 1.0 / 3.0], abs=2.0e-5)
    with pytest.raises(ValueError, match="coordinate length"):
        trapezoid_integral_axis(values, coordinate[:-1], axis=1)
    with pytest.raises(ValueError, match="axis"):
        trapezoid_integral_axis(values, coordinate, axis=3)


def test_trapezoid_integral_supports_legacy_numpy_name(monkeypatch) -> None:
    current = np.trapezoid
    monkeypatch.setattr(np, "trapezoid", None)
    monkeypatch.setattr(np, "trapz", current, raising=False)
    coordinate = np.linspace(0.0, 1.0, 101)
    assert trapezoid_integral(coordinate, coordinate) == 0.5


def test_independent_time_derivatives_recover_an_exact_polynomial_interior() -> None:
    time = np.linspace(-2.0, 2.0, 81)
    values = time**5 - 2.0 * time**3 + 3.0 * time
    expected_third = 60.0 * time**2 - 12.0
    local = local_polynomial_time_derivative(
        time,
        values,
        3,
        window_duration=0.6,
        polynomial_order=5,
    )
    spline = interpolating_spline_time_derivative(
        time,
        values,
        3,
        spline_degree=5,
    )
    interior = np.abs(time) <= 1.5
    np.testing.assert_allclose(local[interior], expected_third[interior], atol=2.0e-9)
    np.testing.assert_allclose(spline[interior], expected_third[interior], atol=2.0e-9)


def test_local_polynomial_sine_derivative_converges_when_sampling_is_halved() -> None:
    errors = []
    for interval in (0.08, 0.04, 0.02):
        time = np.arange(0.0, 8.0 + interval / 2.0, interval)
        numeric = local_polynomial_time_derivative(
            time,
            np.sin(1.7 * time),
            3,
            window_duration=8.0 * interval,
            polynomial_order=5,
        )
        interior = (time >= 1.0) & (time <= 7.0)
        exact = -(1.7**3) * np.cos(1.7 * time[interior])
        errors.append(float(np.sqrt(np.mean((numeric[interior] - exact) ** 2))))
    assert errors[1] < errors[0]
    assert errors[2] < errors[1]


@pytest.mark.parametrize(
    ("call", "message"),
    [
        (
            lambda: local_polynomial_time_derivative(
                [0.0, 0.5, 1.1],
                [0.0, 1.0, 2.0],
                1,
                window_duration=1.0,
            ),
            "uniformly sampled",
        ),
        (
            lambda: local_polynomial_time_derivative(
                np.linspace(0.0, 1.0, 11),
                np.ones(11),
                3,
                window_duration=0.2,
                polynomial_order=2,
            ),
            "at least derivative_order",
        ),
        (
            lambda: interpolating_spline_time_derivative(
                np.linspace(0.0, 1.0, 11),
                np.ones(11),
                6,
                spline_degree=5,
            ),
            "cannot exceed",
        ),
    ],
)
def test_invalid_sampled_time_derivative_inputs_are_rejected(call, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        call()
