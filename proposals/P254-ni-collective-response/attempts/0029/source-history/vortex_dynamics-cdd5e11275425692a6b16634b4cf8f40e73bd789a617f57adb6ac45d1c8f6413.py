"""Periodic point-vortex dynamics: exact solutions of 2-D Euler.

Conventions (all arrays float64):

- N point vortices, circulations ``Gamma_i``, periodic square box ``L``;
- streamfunction per vortex
  ``psi = -Gamma_j/(4 pi) ln[cosh(2 pi dx/L) - cos(2 pi dy/L)]``;
- velocity ``u = (d psi/dy, -d psi/dx)``, self-interaction excluded;
- Hamiltonian ``H = -1/(4 pi) sum_{i<j} Gamma_i Gamma_j ln D_ij``;
- controls: circulation, angular impulse modulo the box, relative H drift;
- deviatoric kinetic-stress observable: deviatoric part of the vortex-set
  velocity covariance, normalized by its initial magnitude; sustained
  strain-coupled restoring stress requires it not to decay.

Importing this module executes no simulation.
"""

from __future__ import annotations

from typing import Any

import numpy as np


def pair_separations(positions: np.ndarray, box: float) -> tuple[np.ndarray, np.ndarray]:
    """Minimum-image dx, dy matrices between all vortices."""

    delta = positions[:, None, :] - positions[None, :, :]
    delta -= box * np.round(delta / box)
    return delta[:, :, 0], delta[:, :, 1]


def velocities(
    positions: np.ndarray, circulations: np.ndarray, box: float
) -> np.ndarray:
    """Induced velocities from the closed-form periodic kernel."""

    dx, dy = pair_separations(positions, box)
    denom = np.cosh(2.0 * np.pi * dx / box) - np.cos(2.0 * np.pi * dy / box)
    np.fill_diagonal(denom, np.inf)
    common = circulations[None, :] / (2.0 * box * denom)
    u_x = -np.sum(common * np.sin(2.0 * np.pi * dy / box), axis=1)
    u_y = np.sum(common * np.sinh(2.0 * np.pi * dx / box), axis=1)
    return np.stack([u_x, u_y], axis=1)


def hamiltonian(
    positions: np.ndarray, circulations: np.ndarray, box: float
) -> float:
    """Periodic point-vortex Hamiltonian (self terms excluded)."""

    dx, dy = pair_separations(positions, box)
    denom = np.cosh(2.0 * np.pi * dx / box) - np.cos(2.0 * np.pi * dy / box)
    np.fill_diagonal(denom, np.inf)
    log_d = np.log(denom)
    iu = np.triu_indices(len(circulations), k=1)
    weights = circulations[iu[0]] * circulations[iu[1]]
    return float(-np.sum(weights * log_d[iu]) / (4.0 * np.pi))


def rk4_step(
    positions: np.ndarray,
    circulations: np.ndarray,
    box: float,
    dt: float,
) -> np.ndarray:
    """One classical RK4 step of the point-vortex equations."""

    def flow(pos: np.ndarray) -> np.ndarray:
        return velocities(pos, circulations, box)

    k1 = flow(positions)
    k2 = flow(positions + 0.5 * dt * k1)
    k3 = flow(positions + 0.5 * dt * k2)
    k4 = flow(positions + dt * k3)
    return positions + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def impulse(positions: np.ndarray, circulations: np.ndarray) -> np.ndarray:
    """Angular impulse P = sum Gamma_i r_i (conserved modulo the box)."""

    return circulations[:, None] * positions

def evolve(
    positions: np.ndarray,
    circulations: np.ndarray,
    box: float,
    dt: float,
    steps: int,
    sample_every: int = 10,
) -> dict[str, Any]:
    """Integrate; return trajectory samples and conservation diagnostics."""

    trajectory = [positions.copy()]
    h_history = [hamiltonian(positions, circulations, box)]
    for step in range(steps):
        positions = np.mod(rk4_step(positions, circulations, box, dt), box)
        if (step + 1) % sample_every == 0:
            trajectory.append(positions.copy())
            h_history.append(hamiltonian(positions, circulations, box))
    reference = h_history[0]
    drift = max(abs(h - reference) for h in h_history) / max(abs(reference), 1e-300)
    return {
        "trajectory": trajectory,
        "hamiltonian_history": h_history,
        "relative_h_drift": drift,
        "final_positions": trajectory[-1],
    }


def deviatoric_kinetic_stress(
    positions: np.ndarray, circulations: np.ndarray, box: float
) -> np.ndarray:
    """Deviatoric velocity covariance of the vortex set (2x2 matrix)."""

    vel = velocities(positions, circulations, box)
    centered = vel - vel.mean(axis=0)
    covariance = centered.T @ centered / len(vel)
    trace = float(np.trace(covariance))
    return covariance - 0.5 * trace * np.eye(2)


def stress_projection_record(
    positions: np.ndarray,
    circulations: np.ndarray,
    box: float,
    dt: float,
    steps: int,
    sample_every: int = 10,
) -> dict[str, Any]:
    """Persistence of initially-imposed anisotropy along one run.

    ``c(t) = <Pi(t), Pi(0)> / |Pi(0)|**2`` projects the current deviatoric
    kinetic stress on its initial direction. An elastic medium with
    memory keeps ``c`` near 1 (oscillating about it under a restoring
    stress); an ergodic fluid decorrelates, flipping sign without a
    preferred restored state.
    """

    result = evolve(positions, circulations, box, dt, steps, sample_every)
    initial = deviatoric_kinetic_stress(positions, circulations, box)
    norm_sq = float(np.sum(initial * initial))
    projections = [
        float(np.sum(deviatoric_kinetic_stress(pos, circulations, box) * initial))
        / max(norm_sq, 1e-300)
        for pos in result["trajectory"]
    ]
    return {
        "projections": projections,
        "initial_projection": projections[0],
        "final_projection": projections[-1],
        "relative_h_drift": result["relative_h_drift"],
    }


def ensemble_projection_record(
    box: float,
    vortex_count: int,
    seeds: int,
    dt: float,
    steps: int,
    sample_every: int = 10,
) -> dict[str, Any]:
    """Ensemble statistics of anisotropy persistence for random clouds."""

    finals, sign_flips = [], 0
    drifts = []
    for seed in range(seeds):
        rng = np.random.default_rng(1000 + seed)
        positions = rng.uniform(0.0, box, size=(vortex_count, 2))
        circulations = np.ones(vortex_count)
        record = stress_projection_record(
            positions, circulations, box, dt, steps, sample_every
        )
        finals.append(record["final_projection"])
        drifts.append(record["relative_h_drift"])
        flips = [
            a * b < 0.0
            for a, b in zip(record["projections"], record["projections"][1:])
        ]
        if any(flips):
            sign_flips += 1
    return {
        "final_projections": finals,
        "median_abs_final": float(np.median(np.abs(finals))),
        "sign_flip_fraction": sign_flips / seeds,
        "max_relative_h_drift": float(np.max(drifts)),
    }


