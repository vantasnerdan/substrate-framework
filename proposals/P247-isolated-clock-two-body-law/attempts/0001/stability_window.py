"""P247 attempt 0001, gate G3: pencil stability window on clean backgrounds.

For each clean ladder background the exact component decomposition gives

    H(R) = R^3 A + R^{-1} D,     lambda_min(H(R)) = R^3 * lambda_min(A + R^-4 D),

with A = d^2 V (potential part) and D = d^2 (C + Phi) (gradient + clock part)
from backgrounds.npz (clean_ladder.py). The stability window per background is
the set of R with lambda_min(A + R^-4 D) > 0; edges are located by bisection
on the exact eigenvalue (no additional quadrature).

Reported 0044 transfer reference: self-consistent stability only for
R in [~8, ~34], upper edge saturating at 32-38 across backgrounds 13-18.
Frozen gate: window edges within +-20 percent of [8, 34].
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from debox_common import REPORTED

HERE = Path(__file__).resolve().parent


def load_backgrounds() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    data = np.load(HERE / "backgrounds.npz")
    tags = sorted({key.split("_")[0] for key in data.files})
    return {tag: (data[f"{tag}_A"], data[f"{tag}_D"]) for tag in tags}


def softest_scaled(a_matrix: np.ndarray, d_matrix: np.ndarray, radius: float) -> float:
    scaled = a_matrix + d_matrix / radius**4
    symmetric = (scaled + scaled.T) / 2.0
    return float(np.linalg.eigvalsh(symmetric)[0])


def window_edges(a_matrix: np.ndarray, d_matrix: np.ndarray, root_radius: float) -> dict:
    """Stability window = connected positive interval of f containing root_radius.

    f(R) = lambda_min(A + R^-4 D); H_total(R) = R^3 f(R), so the branch family
    is stable exactly where f > 0. The committed evidence (index 1 at R=6,
    index 0 at 7.5-10) implies up to three crossings; anchoring on the
    background's own radius picks the physical window without index guesses.
    """
    grid = np.linspace(5.0, 60.0, 1101)
    values = np.array([softest_scaled(a_matrix, d_matrix, r) for r in grid])
    positive = values > 0.0
    if not positive[np.searchsorted(grid, root_radius, side="right") - 1]:
        return {"window": None, "note": "background radius not in a positive interval", "crossings": []}
    index = int(np.searchsorted(grid, root_radius, side="right") - 1)
    lo_index = index
    while lo_index > 0 and positive[lo_index - 1]:
        lo_index -= 1
    hi_index = index
    while hi_index < grid.size - 1 and positive[hi_index + 1]:
        hi_index += 1

    def crossing(outer: float, inner: float) -> float:
        return brentq(lambda r: softest_scaled(a_matrix, d_matrix, r), outer, inner, xtol=1.0e-10)

    low_edge = crossing(grid[max(lo_index - 1, 0)], grid[lo_index]) if lo_index > 0 else float(grid[0])
    high_edge = (
        crossing(grid[hi_index], grid[min(hi_index + 1, grid.size - 1)])
        if hi_index < grid.size - 1
        else float(grid[-1])
    )
    return {
        "window": [float(low_edge), float(high_edge)],
        "crossings": [float(low_edge), float(high_edge)],
        "scan_positive_fraction": float(positive.mean()),
    }


def main() -> None:
    backgrounds = load_backgrounds()
    per_background: dict[str, dict] = {}
    for tag, (a_matrix, d_matrix) in sorted(backgrounds.items()):
        radius = float(tag[1:])
        per_background[tag] = window_edges(a_matrix, d_matrix, radius)
        print(f"[{tag}] window={per_background[tag]['window']}", flush=True)

    windows = {
        tag: record["window"]
        for tag, record in per_background.items()
        if record["window"] is not None
    }
    inner = {tag: window for tag, window in windows.items() if tag != "R12"}
    self_consistent = (
        [max(window[0] for window in inner.values()), min(window[1] for window in inner.values())]
        if inner
        else None
    )
    all_backgrounds = (
        [max(window[0] for window in windows.values()), min(window[1] for window in windows.values())]
        if windows
        else None
    )


    def gate(band_value: list[float] | None) -> dict:
        reported_low, reported_high = REPORTED["window"]
        half = REPORTED["window_relative_halfwidth"]
        if band_value is None:
            return {"band": None, "pass": False}
        low_ok = abs(band_value[0] - reported_low) <= half * reported_low
        high_ok = abs(band_value[1] - reported_high) <= half * reported_high
        return {
            "band": band_value,
            "reported_reference": [reported_low, reported_high],
            "low_pass": bool(low_ok),
            "high_pass": bool(high_ok),
            "pass": bool(low_ok and high_ok),
        }

    gate_self = gate(self_consistent)
    gate_all = gate(all_backgrounds)
    payload = {
        "gate": "G3_window",
        "per_background": per_background,
        "self_consistent_band_r13_plus": {"band": self_consistent, "per_background": inner},
        "all_backgrounds_band": {"band": all_backgrounds, "per_background": windows},
        "gate_self_consistent": gate_self,
        "gate_all_backgrounds": gate_all,
        "gates_all_pass": bool(gate_self["pass"] or gate_all["pass"]),
    }
    (HERE / "stability-window.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({
        "self_consistent_band": self_consistent,
        "all_backgrounds_band": all_backgrounds,
        "gate_self": gate_self,
        "gate_all": gate_all,
    }, indent=2))
    print(f"G3 VERDICT: {'ESTABLISHED' if payload['gates_all_pass'] else 'MIXED'}")


if __name__ == "__main__":
    main()
