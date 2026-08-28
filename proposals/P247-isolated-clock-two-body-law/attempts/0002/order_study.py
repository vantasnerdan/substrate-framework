"""P247 attempt 0002: basis-order study 16/20/24 for the floor eigenvalues.

At R = 12, 14, 16, 18 the root is re-solved at basis orders 20 and 24 (seed =
least-squares Chebyshev reprojection of the order-16 root, exact for order
>= 16), certified at scaled quadrature, and lambda_min(A) is extrapolated in
basis order. Order-16 rows are reused verbatim from the attempt-0001 g5
audit (identical protocol: residuals vs raw operator, jitter scales 1e-13 and
1e-11, 3 draws, seed 20260828, grids 72x24 / 144x48).

Extrapolation model: for the Chebyshev modal basis the truncation error is
exponential in the order N, so when the successive differences grow
(ratio >= 1) the three-order data are fit exactly to
lambda(N) = lambda_inf + A*exp(-b*N) with b = ln(ratio)/4, and the power-law
observed order p is solved for comparison. When the differences shrink
(0 < ratio < 1) no single-exponential fit exists; the conservative geometric
extension bounds the remaining shift by d23/(1-ratio), so lambda_inf lies in
a sign-aware interval between the order-24 value and that bound. The
extrapolation is quoted only when the two differences share a sign and the
newer one exceeds the jitter spread (asymptotic-regime requirement).

Verdicts follow manifest floor_verdict_rule, applied to the conservative
interval: positive above budget upgrades the floor claim to positivity within
the radial family; negative above budget voids it and re-adjudicates G2/G3;
otherwise flatness-to-floor stands.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
PREV = HERE.parent / "0001"
sys.path.insert(0, str(PREV))

import debox_common as base  # noqa: E402  (installs P240 0041/0042 paths, pins threads)

import solve_radial_1d  # noqa: E402

RADII = [12.0, 14.0, 16.0, 18.0]
NEW_ORDERS = (20, 24)
JITTER_SCALES = (1.0e-13, 1.0e-11)
JITTER_DRAWS = 3
JITTER_SEED = 20260828
REL_GRAD_GATE = 1.0e-10
ALIAS_GATE = 1.0e-6
SOLVE_NODES = {20: (60, 30), 24: (72, 36)}
ALIAS_NODES = {20: (120, 60), 24: (144, 72)}
CERT_PRIMARY = {16: (72, 24), 20: (90, 30), 24: (108, 36)}
CERT_DOUBLED = {16: (144, 48), 20: (180, 60), 24: (216, 72)}


def reproject(flat16: np.ndarray, new_order: int) -> np.ndarray:
    """Least-squares reprojection of the modal coordinate to a larger basis.

    The modal coordinate m_k(x) = sum_i c_i T_i(2x^2 - 1) is a Cheb-16
    polynomial in the modal argument, so the fit to a basis of order
    >= 16 is exact up to roundoff; the asserted residual gates that.
    """
    old = np.asarray(flat16, dtype=np.float64).reshape(3, 16)
    x = np.linspace(0.02, 0.98, 4001)
    coord = 2.0 * x**2 - 1.0

    def basis(order: int) -> np.ndarray:
        angle = np.arccos(np.clip(coord, -1.0, 1.0))
        return np.cos(np.multiply.outer(angle, np.arange(order)))

    b_old = basis(16)
    b_new = basis(new_order)
    out = np.zeros((3, new_order), dtype=np.float64)
    worst = 0.0
    for channel in range(3):
        sampled = b_old @ old[channel]
        out[channel], *_ = np.linalg.lstsq(b_new, sampled, rcond=None)
        residual = float(np.max(np.abs(b_new @ out[channel] - sampled)))
        worst = max(worst, residual)
    if worst > 1.0e-12:
        raise RuntimeError(f"reprojection residual {worst:.2e} exceeds roundoff gate")
    return out.ravel()


def solve_at(order: int, radius: float, seed: np.ndarray):
    nodes, angular = SOLVE_NODES[order]
    settings = dict(
        radial_order=order, radial_nodes=nodes, angular_nodes=angular, radius=radius
    )
    return solve_radial_1d.solve_order(order, seed, settings)


def alias_energy(values: np.ndarray, order: int, radius: float) -> float:
    nodes, angular = ALIAS_NODES[order]
    oracle = solve_radial_1d.Oracle(
        dict(radial_order=order, radial_nodes=nodes, angular_nodes=angular, radius=radius)
    )
    return oracle.evaluate(np.asarray(values, dtype=np.float64))[0]


def soft_metrics(matrix: np.ndarray) -> dict:
    raw = np.asarray(matrix, dtype=np.float64)
    symmetric = (raw + raw.T) / 2.0
    asymmetry = float(np.max(np.abs(raw - symmetric)) / max(1.0, float(np.max(np.abs(raw)))))
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    lam_min = float(eigenvalues[0])
    lam_2 = float(eigenvalues[1])
    vector = eigenvectors[:, 0]
    residual = float(
        np.linalg.norm(raw @ vector - lam_min * vector) / max(abs(lam_min), 1e-300)
    )
    return {
        "lambda_min": lam_min,
        "lambda_2": lam_2,
        "ratio_lambda_min_over_lambda_2": lam_min / lam_2,
        "softest_relative_residual": residual,
        "asymmetry_relative": asymmetry,
        "morse_gate": base.committed_morse_gate(eigenvalues),
    }


def jitter_spread(values: np.ndarray, order: int, nodes: tuple[int, int]) -> dict:
    rng = np.random.default_rng(JITTER_SEED)
    flats = [
        np.asarray(values, dtype=np.float64)
        + scale
        * np.linalg.norm(values)
        * rng.standard_normal(values.size)
        / np.sqrt(values.size)
        for scale in JITTER_SCALES
        for _ in range(JITTER_DRAWS)
    ]
    mins = []
    for flat in flats:
        a_matrix, _ = base.component_hessians(flat, nodes[0], nodes[1], order=order)
        mins.append(base.lambda_min(a_matrix)[0])
    return {
        "nodes": list(nodes),
        "min": float(np.min(mins)),
        "max": float(np.max(mins)),
        "spread": float(np.max(mins) - np.min(mins)),
    }


def extrapolate(order_list: list[int], lam_list: list[float], jitter: float) -> dict:
    """Three-order extrapolation of lambda_min in basis order N (see module doc)."""
    d12 = lam_list[0] - lam_list[1]
    d23 = lam_list[1] - lam_list[2]
    out = {
        "orders": order_list,
        "lambda_by_order": lam_list,
        "d_16_20": d12,
        "d_20_24": d23,
    }
    if d12 == 0.0 or d23 == 0.0 or (d12 > 0.0) != (d23 > 0.0):
        out["asymptotic_regime"] = False
        out["continuum_estimate"] = None
        return out
    ratio = d23 / d12
    out["asymptotic_regime"] = bool(abs(d23) > jitter)
    if ratio >= 1.0:
        b = float(np.log(ratio) / 4.0)
        weight = float(np.exp(-b * order_list[2])) / (
            float(np.exp(-b * order_list[1])) - float(np.exp(-b * order_list[2]))
        )
        lam_inf = float(lam_list[2] - d23 * weight)

        def power_residual(p: float) -> float:
            return (order_list[0] ** (-p) - order_list[1] ** (-p)) / (
                order_list[1] ** (-p) - order_list[2] ** (-p)
            ) - ratio

        try:
            p = float(brentq(power_residual, 1e-6, 60.0))
        except Exception:
            p = float("nan")
        out.update(
            model="three_point_exponential",
            exponential_rate_b=b,
            power_law_observed_order=p,
            continuum_estimate=lam_inf,
        )
    else:
        out.update(
            model="geometric_bound_differences_shrinking",
            continuum_estimate=float(lam_list[2] + d23 / (1.0 - ratio)),
        )
    out["continuum_shift_from_order24"] = out["continuum_estimate"] - lam_list[2]
    return out


def verdict_for(per_radius: dict) -> dict:
    lam16 = per_radius["order_16"]["lambda_min"]
    rows = [lam16]
    for order in NEW_ORDERS:
        block = per_radius.get(f"order_{order}", {})
        if block.get("converged"):
            rows.append(block["lambda_min"])
    if len(rows) < 3:
        return {
            "band": "OTHER",
            "mechanism": f"only {len(rows)} converged orders; extrapolation not attemptable",
        }
    jitter20 = per_radius["order_20"]["jitter"]["spread"]
    extrap = extrapolate([16, 20, 24], rows, jitter20)
    budget_parts = []
    for order in NEW_ORDERS:
        block = per_radius.get(f"order_{order}", {})
        if block.get("converged"):
            budget_parts.append(block["quadrature_gauge"])
            budget_parts.append(block["jitter"]["spread"])
    budget = max(budget_parts + [1e-18])
    estimate = extrap.get("continuum_estimate")
    lam24 = rows[2]
    if not extrap.get("asymptotic_regime", False) or estimate is None:
        band, mechanism = (
            "UNRESOLVED_BELOW_BUDGET",
            "successive order differences change sign or sit below jitter; flatness-to-floor stands",
        )
    else:
        lo, hi = min(lam24, estimate), max(lam24, estimate)
        if lo > 10.0 * budget:
            band, mechanism = (
                "POSITIVE_ABOVE_BUDGET",
                "conservative extrapolation interval positive above the error budget; "
                "floor claim upgraded to positivity within the radial family",
            )
        elif hi < -10.0 * budget:
            band, mechanism = (
                "NEGATIVE_ABOVE_GAUGE",
                "conservative extrapolation interval negative above the error budget; "
                "floor claim voided, G2/G3 re-adjudication required",
            )
        else:
            band, mechanism = (
                "UNRESOLVED_BELOW_BUDGET",
                "extrapolation interval within a factor 10 of the budget; flatness-to-floor stands",
            )
    return {"extrapolation": extrap, "error_budget": budget, "band": band, "mechanism": mechanism}


def main() -> None:
    started = time.time()
    with open(PREV / "clean-ladder.json") as handle:
        previous = json.load(handle)
    prior = {row["radius"]: row for row in previous["rungs"] if row.get("accepted")}
    with open(PREV / "g5-small-ratio-audit.json") as handle:
        g5 = json.load(handle)
    g5_by_radius = {row["radius"]: row for row in g5["records"]}

    records: list[dict] = []
    for radius in RADII:
        row16 = prior[radius]
        root16 = np.asarray(row16["values"], dtype=np.float64)
        per_radius: dict = {"radius": radius}
        for order in NEW_ORDERS:
            step_started = time.time()
            seed = reproject(root16, order)
            solution = solve_at(order, radius, seed)
            rel_grad = solution["relative_gradient"]
            if not rel_grad < REL_GRAD_GATE:
                per_radius[f"order_{order}"] = {
                    "converged": False,
                    "relative_gradient": rel_grad,
                }
                print(f"[R={radius} order={order}] FAILED relgrad={rel_grad:.3e}", flush=True)
                continue
            values = np.asarray(solution["values"], dtype=np.float64)
            energy = float(solution["energy"])
            alias = alias_energy(values, order, radius)
            alias_relative = abs(alias - energy) / abs(energy)
            primary = CERT_PRIMARY[order]
            doubled = CERT_DOUBLED[order]
            a_primary, _ = base.component_hessians(values, primary[0], primary[1], order=order)
            a_doubled, _ = base.component_hessians(values, doubled[0], doubled[1], order=order)
            metrics = soft_metrics(a_primary)
            lam_doubled = base.lambda_min(a_doubled)[0]
            quadrature_gauge = abs(lam_doubled - metrics["lambda_min"])
            jit = jitter_spread(values, order, primary)
            per_radius[f"order_{order}"] = {
                "converged": True,
                "relative_gradient": rel_grad,
                "energy": energy,
                "alias_relative_gap": alias_relative,
                "alias_pass": bool(alias_relative <= ALIAS_GATE),
                "cert_nodes": list(primary),
                "doubled_nodes": list(doubled),
                **metrics,
                "lambda_min_doubled": lam_doubled,
                "quadrature_gauge": quadrature_gauge,
                "jitter": jit,
                "minutes": round((time.time() - step_started) / 60.0, 2),
            }
            print(
                f"[R={radius} order={order}] lamA={metrics['lambda_min']:+.6e} "
                f"lam2={metrics['lambda_2']:.3e} "
                f"res={metrics['softest_relative_residual']:.2e} "
                f"quad_gauge={quadrature_gauge:.2e} jit={jit['spread']:.2e} "
                f"alias={alias_relative:.2e} ({per_radius[f'order_{order}']['minutes']} min)",
                flush=True,
            )
        # order-16 row reused from the g5 audit (identical protocol)
        g = g5_by_radius[radius]
        per_radius["order_16"] = {
            "source": "attempts/0001/g5-small-ratio-audit.json",
            "cert_nodes": [72, 24],
            "doubled_nodes": [144, 48],
            "lambda_min": g["lambda_min_A_coarse_72x24"],
            "lambda_2": g["lambda_2_A_coarse_72x24"],
            "ratio_lambda_min_over_lambda_2": g["ratio_lambda_min_over_lambda_2_coarse_72x24"],
            "softest_relative_residual": g["residual_softest_A_coarse_72x24"],
            "asymmetry_relative": g["asym_max_A_coarse_72x24"],
            "quadrature_gauge": g["quadrature_gauge_lambda_min_A"],
            "jitter_spread_1e_11": g["jitter_spread_A_scale_1e-11"],
            "lambda_min_doubled": g["lambda_min_A_fine_144x48"],
        }
        records.append(per_radius)

    # projector-background restricted spectrum across orders (operator truncation gauge)
    projector_gauge: dict = {}
    for order in (16, 20, 24):
        nodes = CERT_PRIMARY[order]
        zeros = np.zeros(3 * order, dtype=np.float64)
        a_projector_ansatz, _ = base.component_hessians(
            zeros, nodes[0], nodes[1], order=order
        )
        spectrum = np.linalg.eigvalsh((a_projector_ansatz + a_projector_ansatz.T) / 2.0)
        projector_gauge[str(order)] = {
            "nodes": list(nodes),
            "background": "zero modal coefficients (committed ansatz background, NOT the rank-1 projector)",
            "lambda_min": float(spectrum[0]),
            "lambda_2": float(spectrum[1]),
            "lambda_max": float(spectrum[-1]),
        }
        print(
            f"[ansatz background order={order}] restricted spectrum "
            f"lam_min={spectrum[0]:.6e} lam2={spectrum[1]:.6e} lam_max={spectrum[-1]:.6e}",
            flush=True,
        )

    verdicts = {str(per_radius["radius"]): verdict_for(per_radius) for per_radius in records}
    for radius, verdict in verdicts.items():
        print(f"[verdict R={radius}] {verdict['band']}: {verdict['mechanism']}", flush=True)

    payload = {
        "records": records,
        "ansatz_background_gauge": projector_gauge,
        "verdicts": verdicts,
        "minutes_total": round((time.time() - started) / 60.0, 2),
    }
    (HERE / "g6-order-study.json").write_text(json.dumps(payload, indent=2))
    print("WROTE g6-order-study.json")


if __name__ == "__main__":
    main()
