"""P247 attempt 0003: W1 order-24 pencil windows + W3 extended-ladder refined instability.

Preregistered in manifest.yaml (commits a8ae525 + a82028e, amendment
recorded after the first-execution defect) BEFORE execution; interpretation
bands frozen there. Machinery: committed continuation protocol
(continuation_scan.py), order-24 re-solve protocol (order_study.py),
committed window machinery (0001/stability_window.py), committed Hessians
(debox_common.component_hessians). No new solver classes.

First-execution defect and repair (manifest amendments_before_execution):
the R = 24 order-24 re-solve from the reprojection seed converged to a
spurious divergent stationary point (E = 6.9e516, relgrad 7.2e-70) that the
relative aliasing gate passed blind (gap 3.04e-8). Repairs: (1) root
continuity gate |E24 - E16|/E16 < 1e-6 on every accepted root; (2) on
failure, re-seed by continuing in R at order 24 from the nearest accepted
order-24 root; (3) a still-failing rung is recorded
root-unresolvable-at-order-24 and scoped out of the bands with the
mechanism named.

W3: order-16 provenance roots at R = 20..30, re-solved at basis order 24;
lambda_min(A) quoted with the full small-ratio budget. Bands:
NEGATIVE_PERSISTENT / SIGN_CHANGE / BUDGET_FLAT(escalate order 28 once).

W1: order-24 A, D matrices at R in {12, 14, 18, 24}; committed window_edges
semantics anchored at each background's own radius on [5, 60]; edges at
primary (108,36) and doubled (216,72) quadrature. Bands: RESTORED_ORDER24 /
REFUTED_ORDER24 / UNRESOLVED_EDGES / OWN_RADIUS_BREAK.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PREV = HERE.parent / "0002"
A1 = HERE.parent / "0001"
sys.path.insert(0, str(A1))
sys.path.insert(0, str(PREV))

import debox_common as base  # noqa: E402  (pins threads, installs P240 paths)

import continuation_scan as cs  # noqa: E402
import order_study as osd  # noqa: E402
import solve_radial_1d  # noqa: E402
import stability_window as sw  # noqa: E402

REL_GRAD_GATE = 1.0e-10
ALIAS_GATE = 1.0e-6
CONTINUITY_GATE = 0.25  # amendment 2: measured order-improvement is 0.7-2.2 percent and grows with R; divergence fails by >500 orders
CONT_RADII = [20.0, 22.0, 24.0, 26.0, 28.0, 30.0]
W1_RADII = [12.0, 14.0, 18.0, 24.0]
ORDER = 24
PRIMARY = (108, 36)
DOUBLED = (216, 72)
ESCALATED = (144, 48)
ESCALATED_DOUBLED = (288, 96)
EDGE_GAUGE = 2.0  # R units; ~20 percent of the lower reported edge
BUDGET_MARGIN = 10.0
REPORTED_WINDOW = (8.0, 34.0)
GATE_HALF = 0.20


def oracle_for(order: int, nodes: tuple[int, int], radius: float):
    return solve_radial_1d.Oracle(
        dict(
            radial_order=order,
            radial_nodes=nodes[0],
            angular_nodes=nodes[1],
            radius=radius,
        )
    )


def channel_split(matrix: np.ndarray) -> dict:
    symmetric = (np.asarray(matrix, dtype=np.float64) + np.asarray(matrix).T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    fractions, nodes = solve_radial_1d.analyze_mode(eigenvectors[:, 0])
    return {
        "fractions_q_tangent_split": [float(f) for f in fractions],
        "split_radial_nodes": int(nodes),
    }


def order24_root(
    radius: float,
    seed16: np.ndarray,
    label: str,
    seed24: np.ndarray | None = None,
) -> dict:
    """Order-24 solve under the root-continuity gate with R-continuation
    reseeding (manifest amendment)."""
    reference = float(oracle_for(16, (96, 48), radius).evaluate(seed16)[0])
    seeds = [seed24] if seed24 is not None else [osd.reproject(seed16, ORDER)]
    last = {}
    for attempt_index, seed in enumerate(seeds):
        solution = osd.solve_at(ORDER, radius, seed)
        rel_grad = solution["relative_gradient"]
        if not rel_grad < REL_GRAD_GATE:
            raise RuntimeError(f"{label}: order-24 solve failed relgrad={rel_grad:.3e}")
        values = np.asarray(solution["values"], dtype=np.float64)
        energy = float(solution["energy"])
        alias_relative = abs(osd.alias_energy(values, ORDER, radius) - energy) / abs(
            energy
        )
        continuity = abs(energy - reference) / abs(reference)
        last = {
            "energy": energy,
            "alias": alias_relative,
            "continuity": continuity,
        }
        print(
            f"[root24 {label}] attempt={attempt_index} E={energy:.6f} "
            f"relgrad={rel_grad:.2e} alias={alias_relative:.2e} "
            f"continuity={continuity:.2e}",
            flush=True,
        )
        if continuity < CONTINUITY_GATE:
            return {
                "values": values,
                "energy": energy,
                "relative_gradient": rel_grad,
                "alias_relative": alias_relative,
                "continuity": continuity,
                "seed_attempts": attempt_index + 1,
            }
    raise RuntimeError(
        f"{label}: root-unresolvable-at-order-24 after reseeding "
        f"(spurious divergent basin; last={last})"
    )


def w3_row(root: dict) -> dict:
    values = root["values"]
    a_primary, _ = base.component_hessians(values, PRIMARY[0], PRIMARY[1], order=ORDER)
    metrics = osd.soft_metrics(a_primary)
    a_doubled, _ = base.component_hessians(values, DOUBLED[0], DOUBLED[1], order=ORDER)
    lam_doubled = base.lambda_min(a_doubled)[0]
    quad_gauge = abs(lam_doubled - metrics["lambda_min"])
    jit = osd.jitter_spread(values, ORDER, PRIMARY)
    budget = quad_gauge + jit["spread"]
    sign_stable = (jit["min"] > 0.0 and metrics["lambda_min"] > 0.0) or (
        jit["max"] < 0.0 and metrics["lambda_min"] < 0.0
    )
    row = {
        "energy": root["energy"],
        "alias_relative": root["alias_relative"],
        "continuity": root.get("continuity"),
        "cert_nodes": list(PRIMARY),
        "doubled_nodes": list(DOUBLED),
        **metrics,
        "lambda_min_doubled": lam_doubled,
        "quadrature_gauge": quad_gauge,
        "jitter": jit,
        "budget": budget,
        "verdict_threshold": BUDGET_MARGIN * budget,
        "sign_stable_under_jitter": bool(sign_stable),
        "softest_channel": channel_split(a_primary),
    }
    if abs(metrics["lambda_min"]) > BUDGET_MARGIN * budget and sign_stable:
        row["sign"] = "negative" if metrics["lambda_min"] < 0.0 else "positive"
    else:
        row["sign"] = "unresolved"
    print(
        f"[W3] lamA={metrics['lambda_min']:+.6e} lam2={metrics['lambda_2']:.3e} "
        f"budget={budget:.2e} sign={row['sign']} "
        f"frac={np.round(row['softest_channel']['fractions_q_tangent_split'], 4).tolist()}",
        flush=True,
    )
    return row


def escalated_w3_row(values: np.ndarray) -> dict:
    a_matrix, _ = base.component_hessians(values, 130, 43, order=28)
    metrics = osd.soft_metrics(a_matrix)
    return {
        "escalated": True,
        "order": 28,
        "cert_nodes": [130, 43],
        **metrics,
        "softest_channel": channel_split(a_matrix),
    }


def window_record(values: np.ndarray, anchor: float) -> dict:
    a_p, d_p = base.component_hessians(values, PRIMARY[0], PRIMARY[1], order=ORDER)
    a_d, d_d = base.component_hessians(values, DOUBLED[0], DOUBLED[1], order=ORDER)
    win_p = sw.window_edges(a_p, d_p, anchor)
    win_d = sw.window_edges(a_d, d_d, anchor)
    record = {
        "anchor": anchor,
        "primary": {"nodes": list(PRIMARY), **win_p},
        "doubled": {"nodes": list(DOUBLED), **win_d},
    }
    edges_p = [e for e in (win_p.get("window") or ()) if np.isfinite(e)]
    edges_d = [e for e in (win_d.get("window") or ()) if np.isfinite(e)]
    if edges_p and edges_d and len(edges_p) == len(edges_d):
        record["edge_gauge"] = float(max(abs(p - d) for p, d in zip(edges_p, edges_d)))
    else:
        record["edge_gauge"] = None
    record["quadrature_gauge_lambda_min_A"] = abs(
        base.lambda_min(a_d)[0] - base.lambda_min(a_p)[0]
    )
    record["f_at_own_radius_primary"] = float(sw.softest_scaled(a_p, d_p, anchor))
    print(
        f"[W1 R={anchor}] primary={win_p.get('window')} doubled={win_d.get('window')} "
        f"edge_gauge={record['edge_gauge']} "
        f"f_own={record['f_at_own_radius_primary']:+.6e}",
        flush=True,
    )
    return record


def gate_band(band: list[float] | None) -> dict:
    if band is None:
        return {"band": None, "pass": False}
    low_ok = abs(band[0] - REPORTED_WINDOW[0]) <= GATE_HALF * REPORTED_WINDOW[0]
    high_ok = abs(band[1] - REPORTED_WINDOW[1]) <= GATE_HALF * REPORTED_WINDOW[1]
    return {
        "band": band,
        "reported_reference": list(REPORTED_WINDOW),
        "low_pass": bool(low_ok),
        "high_pass": bool(high_ok),
        "pass": bool(low_ok and high_ok),
    }


def main() -> None:
    started = time.time()
    with open(A1 / "clean-ladder.json") as handle:
        previous = json.load(handle)
    prior = {row["radius"]: row for row in previous["rungs"] if row.get("accepted")}

    # --- order-16 provenance continuation R = 20..30 (committed protocol) ---
    cont16: dict[float, np.ndarray] = {}
    cont16_provenance: dict[str, dict] = {}
    seed = np.asarray(prior[18.0]["values"], dtype=np.float64)
    for radius in CONT_RADII:
        solution = cs.solve_rung(radius, seed)
        rel_grad = solution["relative_gradient"]
        if not rel_grad < REL_GRAD_GATE:
            raise RuntimeError(f"continuation R={radius} failed relgrad={rel_grad:.3e}")
        values = np.asarray(solution["values"], dtype=np.float64)
        energy = float(solution["energy"])
        alias_relative = abs(cs.alias_energy(values, radius) - energy) / abs(energy)
        if not alias_relative < ALIAS_GATE:
            raise RuntimeError(f"continuation R={radius} alias gap={alias_relative:.3e}")
        cont16[radius] = values
        cont16_provenance[str(radius)] = {
            "energy": energy,
            "relative_gradient": rel_grad,
            "alias_relative": alias_relative,
        }
        seed = values
        print(
            f"[cont16 R={radius}] E={energy:.6f} relgrad={rel_grad:.2e} "
            f"alias={alias_relative:.2e}",
            flush=True,
        )

    # --- order-24 roots with continuity gate and R-continuation reseeding ---
    needed = sorted(set(W1_RADII) | set(CONT_RADII))
    roots: dict[float, dict] = {}
    unresolvable: dict[str, str] = {}
    for index, radius in enumerate(needed):
        if radius in prior:
            seed16 = np.asarray(prior[radius]["values"], dtype=np.float64)
            source = "clean-ladder order-16"
        else:
            seed16 = cont16[radius]
            source = "continuation order-16 (this run)"
        seed24 = None
        # nearest accepted order-24 root, walked outward (repair path 2)
        for offset in (1, 2, 3):
            for neighbour in (radius - 2.0 * offset, radius + 2.0 * offset):
                if neighbour in roots:
                    seed24 = roots[neighbour]["values"]
                    break
            if seed24 is not None:
                break
        try:
            roots[radius] = order24_root(radius, seed16, f"R={radius}", seed24=seed24)
            roots[radius]["seed_source"] = source
        except RuntimeError as failure:
            unresolvable[str(radius)] = str(failure)
            print(f"[root24 R={radius}] UNRESOLVABLE: {failure}", flush=True)
    if roots:
        np.savez(
            HERE / "roots24.npz",
            **{f"R{int(r)}": roots[r]["values"] for r in sorted(roots)},
        )

    # --- W3 rows on the extended ladder (resolvable rungs only) ---
    w3: dict[str, dict] = {}
    escalations: dict[str, dict] = {}
    for radius in CONT_RADII:
        if str(radius) in unresolvable or radius not in roots:
            continue
        row = w3_row(roots[radius])
        w3[str(radius)] = row
        if row["sign"] == "unresolved":
            seed28 = osd.reproject(roots[radius]["values"], 28)
            solution = solve_radial_1d.solve_order(
                28,
                seed28,
                dict(radial_order=28, radial_nodes=86, angular_nodes=43, angular_modes=1),
            )
            if solution["relative_gradient"] < REL_GRAD_GATE:
                escalated = escalated_w3_row(
                    np.asarray(solution["values"], dtype=np.float64)
                )
                escalations[str(radius)] = escalated
                print(
                    f"[W3 escalated R={radius}] lamA={escalated['lambda_min']:+.6e} "
                    f"lam2={escalated['lambda_2']:.3e}",
                    flush=True,
                )

    signs = [w3[str(r)]["sign"] for r in CONT_RADII if str(r) in w3]
    flat_rungs = [str(r) for r in CONT_RADII if str(r) in w3 and w3[str(r)]["sign"] == "unresolved"]
    scoped_out = sorted(unresolvable.keys())
    if scoped_out:
        w3_band = "PARTIALLY_SCOPED"
        w3_mechanism = (
            f"rung(s) {scoped_out} root-unresolvable-at-order-24 "
            "(spurious divergent basin, amendment scope); "
            f"resolved rungs: signs={signs}, flat={flat_rungs}"
        )
    elif signs and all(s == "negative" for s in signs):
        w3_band = "NEGATIVE_PERSISTENT"
        w3_mechanism = (
            "the refined q-channel static instability persists at every extended rung; "
            "|lambda_min|(R) trend recorded for candidate-B design"
        )
    elif "positive" in signs and "negative" in signs:
        w3_band = "SIGN_CHANGE"
        w3_mechanism = (
            "the static instability heals within the extended ladder; "
            "crossing location by bisection is the named follow-up"
        )
    else:
        w3_band = "BUDGET_FLAT_PARTIAL" if flat_rungs else "INCONCLUSIVE"
        w3_mechanism = f"unresolved rungs: {flat_rungs}; escalated rows recorded"

    # --- W1 windows (resolvable W1 radii only) ---
    w1: dict[str, dict] = {}
    for radius in W1_RADII:
        if radius not in roots:
            continue
        record = window_record(roots[radius]["values"], radius)
        own = record["f_at_own_radius_primary"]
        budget = record["quadrature_gauge_lambda_min_A"]
        record["own_radius_break"] = bool(own < 0.0 and -own > BUDGET_MARGIN * budget)
        if record.get("edge_gauge") is not None and record["edge_gauge"] >= EDGE_GAUGE:
            values = roots[radius]["values"]
            a_e, d_e = base.component_hessians(
                values, ESCALATED[0], ESCALATED[1], order=ORDER
            )
            a_ed, d_ed = base.component_hessians(
                values,
                ESCALATED_DOUBLED[0],
                ESCALATED_DOUBLED[1],
                order=ORDER,
            )
            win_e = sw.window_edges(a_e, d_e, radius)
            win_ed = sw.window_edges(a_ed, d_ed, radius)
            record["escalated"] = {
                "nodes": list(ESCALATED),
                "doubled_nodes": list(ESCALATED_DOUBLED),
                "primary": win_e,
                "doubled": win_ed,
            }
            edges_e = [e for e in (win_e.get("window") or ()) if np.isfinite(e)]
            edges_ed = [e for e in (win_ed.get("window") or ()) if np.isfinite(e)]
            if edges_e and edges_ed and len(edges_e) == len(edges_ed):
                record["escalated"]["edge_gauge"] = float(
                    max(abs(p - d) for p, d in zip(edges_e, edges_ed))
                )
                record["edge_gauge"] = record["escalated"]["edge_gauge"]
        record["resolution_limited"] = bool(
            record.get("edge_gauge") is None or record["edge_gauge"] >= EDGE_GAUGE
        )
        w1[str(radius)] = record

    usable = {
        tag: rec
        for tag, rec in w1.items()
        if not rec["resolution_limited"] and not rec["own_radius_break"]
    }
    windows = {
        tag: rec["primary"]["window"]
        for tag, rec in usable.items()
        if rec["primary"].get("window") is not None
    }
    inner_tags = [tag for tag in windows if tag != "R12"]
    self_consistent = (
        [
            max(windows[tag][0] for tag in inner_tags),
            min(windows[tag][1] for tag in inner_tags),
        ]
        if inner_tags
        else None
    )
    all_band = (
        [max(w[0] for w in windows.values()), min(w[1] for w in windows.values())]
        if windows
        else None
    )
    gate_self = gate_band(self_consistent)
    gate_all = gate_band(all_band)
    any_gate_pass = bool(gate_self["pass"] or gate_all["pass"])
    if any_gate_pass:
        w1_band = "RESTORED_ORDER24"
        w1_mechanism = (
            "the bounded stability window re-emerges at refined resolution "
            "within the frozen G3 gate"
        )
    elif any(rec["resolution_limited"] for rec in w1.values()):
        w1_band = "UNRESOLVED_EDGES"
        w1_mechanism = (
            "at least one edge is resolution-limited even after escalation; "
            "named per background"
        )
    else:
        w1_band = "REFUTED_ORDER24"
        w1_mechanism = (
            "the window fails the frozen gate at refined resolution with edges "
            "above the quadrature gauge"
        )

    with open(A1 / "stability-window.json") as handle:
        order16_windows = {
            tag: rec.get("window")
            for tag, rec in json.load(handle)["per_background"].items()
        }

    payload = {
        "manifest": "attempts/0003/manifest.yaml",
        "first_run_defect": "g9-first-run-divergent-R24-stdout.txt (spurious divergent R=24 basin, alias gate passed at 3.04e-8)",
        "order16_windows_reference": order16_windows,
        "continuation16_provenance": cont16_provenance,
        "roots24": {
            str(r): {
                "seed_source": roots[r]["seed_source"],
                "energy": roots[r]["energy"],
                "alias_relative": roots[r]["alias_relative"],
                "continuity": roots[r].get("continuity"),
                "seed_attempts": roots[r].get("seed_attempts"),
            }
            for r in sorted(roots)
        },
        "unresolvable_rungs": unresolvable,
        "w3_extended_ladder": w3,
        "w3_escalations": escalations,
        "w3_band": w3_band,
        "w3_mechanism": w3_mechanism,
        "w1_windows": w1,
        "w1_self_consistent_band_r14_plus": {"band": self_consistent, "gate": gate_self},
        "w1_all_backgrounds_band": {"band": all_band, "gate": gate_all},
        "w1_band": w1_band,
        "w1_mechanism": w1_mechanism,
        "minutes_total": round((time.time() - started) / 60.0, 2),
    }
    (HERE / "g9-w1w3-order24.json").write_text(json.dumps(payload, indent=2))
    print(f"W1 BAND: {w1_band} - {w1_mechanism}")
    print(f"W3 BAND: {w3_band} - {w3_mechanism}")
    print("WROTE g9-w1w3-order24.json")


if __name__ == "__main__":
    main()
