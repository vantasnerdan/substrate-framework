"""Attempt 0010 -- de-boxing ladder Delta_E(R) about window-family roots.

PREREGISTRATION (frozen in proposals/P244-clock-full-band-spectrum/
proposal.yaml `next_loop` BEFORE execution; mirrored here):

QUESTION
  Is delta-E_bare = 72.58859645998888 an intrinsic property of the clock or
  a confining-window artifact? Attempt 0009 showed the empty-window
  subtraction is structurally undefined, so the comparison keeps the clock
  present and varies the box: Delta_E(R) on the window-family root ladder
  R in {10, 12, 14}.

CONSTRUCTION
  Verbatim committed machinery per rung: coefficients read read-only from
  P240 attempts/0042 largeR-roots.json; kinetic_stage2.RADIUS and
  route_a_corrected.RADIUS set to the rung radius so build_cache's gauss
  grid and corrected_kinetic_value's x=r/R basis map the identical code
  path; corrected_kinetic_hessian for the kinetic metric; pencil_float64
  for the spectrum. Gauss-Legendre 96x48 primary at every rung plus a
  48x24 coarse member of the same family as the rung cross-check (the
  Chebyshev 160x80 cross-family stays an R12-only certified result and is
  NOT recomputed here).

  Budget substitution disclosed: C-M5S-011 composes
  sigma_i = max(float64 entry jitter, cross-family spread). At non-R12
  rungs the cross-family member does not exist, so the independent-
  discretization spread |omega(96x48) - omega(48x24)| replaces it inside
  the same max(); the entry-jitter half is unchanged. Certification
  margin stays omega >= 10*sigma; zero-point sum stays fsum over
  certified positive modes with RSS and linear bounds.

FROZEN INTERPRETATION BANDS (fixed before any Delta_E(R>12) existed):
  PLATEAU_INTRINSIC : |D10/D12 - 1| <= 0.05 AND |D14/D12 - 1| <= 0.05
  DECAY_BOX_ARTIFACT: monotone decrease, ratio <= 2/3 per outward rung
  GROWTH_VOLUME     : monotone increase, ratio >= 3/2 per outward rung
  MIXED / certification failure at any rung: INCONCLUSIVE with named
  mechanism; no consumer-mass verdict may be drawn.

GATES PER RUNG
  G0 transfer  : Oracle energy at 96x48 equals the rung's pre-recorded
                 frozen-background reference within 1e-6 relative
                 (references measured in this session before preregistration
                 was written; R12 equals the committed E to 4e-13).
  G0R regression: R12 zero-point reproduces C-M5S-011's 72.58859645998888
                 to <= 1e-9 absolute (same code path; expected exact).
  G3 mutation  : +1e-4 on coefficient index 7 moves some certified
                 frequency relatively more than 1e-7.
  G5 certification: every kept mode certified at 96x48 on every rung;
                 otherwise the rung fails and the verdict is INCONCLUSIVE.

SCOPE
  Order-16 truncated model class throughout; a de-boxing TREND
  measurement, not a continuum renormalization theorem. Consumer-mass
  consequences are applied only through separately reviewed edits.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
for p in (
    "campaigns/P243-clock-sourced-induced-coupling/attempts/0008",
    "proposals/P240-m5-kinetic-axis/attempts/0041",
    "proposals/P244-clock-full-band-spectrum/attempts/0002",
    "proposals/P244-clock-full-band-spectrum/attempts/0003",
):
    sys.path.insert(0, str(REPO / p))

import kinetic_stage2  # noqa: E402
import route_a_corrected  # noqa: E402
import solve_radial_1d  # noqa: E402
from route_a_corrected import CERT_MARGIN, ORDER  # noqa: E402
from route_final import pencil_float64  # noqa: E402

torch.set_num_threads(1)

RUNGS = (
    ("R10", 10.0, 53.193754983),
    ("R12", 12.0, 55.10418278043526),
    ("R14", 14.0, 56.628217935),
)
DE12_COMMITTED = 72.58859645998888
MUT_INDEX = 7
VERDICT_TOL = 0.05


def evaluate_rung(tag: str, radius: float, e_ref: float, roots: dict):
    background = np.asarray(roots[tag]["values"], dtype=float)
    kinetic_stage2.RADIUS = radius
    route_a_corrected.RADIUS = radius

    out = {"rung": tag, "radius": radius}

    oracle = solve_radial_1d.Oracle(dict(radial_order=ORDER,
                                         radial_nodes=96, angular_nodes=48,
                                         radius=radius))
    e_cl, _, h_raw, _ = oracle.evaluate(background)
    e_cl = float(e_cl)
    out["E_classical"] = e_cl
    out["g0_transfer_rel"] = abs(e_cl - e_ref) / e_ref
    out["g0_pass"] = bool(out["g0_transfer_rel"] <= 1e-6)

    h_fin = (np.asarray(h_raw) + np.asarray(h_raw).T) / 2
    m_fin = route_a_corrected.corrected_kinetic_hessian(
        kinetic_stage2.build_cache(background, 96, 48))
    m_fin = (m_fin + m_fin.T) / 2

    oracle_c = solve_radial_1d.Oracle(dict(radial_order=ORDER,
                                           radial_nodes=48, angular_nodes=24,
                                           radius=radius))
    _, _, h_raw_c, _ = oracle_c.evaluate(background)
    h_coarse = (np.asarray(h_raw_c) + np.asarray(h_raw_c).T) / 2
    m_coarse = route_a_corrected.corrected_kinetic_hessian(
        kinetic_stage2.build_cache(background, 48, 24))
    m_coarse = (m_coarse + m_coarse.T) / 2

    omega_f, vecs, stiff_ray, kept_n, m_evals, keep_mask = \
        pencil_float64(h_fin, m_fin)
    omega_c = pencil_float64(h_coarse, m_coarse)[0]
    n_cmp = min(len(omega_f), len(omega_c))
    wa = np.sort(np.abs(np.asarray(omega_f[:n_cmp], dtype=float)))
    wc = np.sort(np.abs(np.asarray(omega_c[:n_cmp], dtype=float)))

    rng = np.random.default_rng(20260825)
    seeds = []
    for _ in range(8):
        hj = h_fin * (1.0 + 2e-15 * rng.standard_normal(h_fin.shape))
        mj = m_fin * (1.0 + 2e-15 * rng.standard_normal(m_fin.shape))
        seeds.append(pencil_float64((hj + hj.T) / 2, (mj + mj.T) / 2)[0][:n_cmp])
    sigma_entry = np.array(seeds).std(axis=0, ddof=1)

    rows = []
    for i in range(n_cmp):
        w = float(math.sqrt(max(abs(wa[i]), 0.0)))
        ladder_sigma_w = float(abs(wc[i] - wa[i]))
        sigma_w = max(float(sigma_entry[i]) / max(2.0 * w, 1e-300),
                      ladder_sigma_w)
        ok = bool(w > CERT_MARGIN * sigma_w)
        rows.append({"mode": i,
                     "stiffness_rayleigh": float(stiff_ray[i]),
                     "omega": w,
                     "sigma_omega": sigma_w,
                     "ladder_abs_gap": ladder_sigma_w,
                     "jitter_sigma_over_2w": float(
                         sigma_entry[i] / max(2.0 * w, 1e-300)),
                     "certified_margin_ok": ok})

    mutated = background.copy()
    mutated[MUT_INDEX] += 1e-4
    oracle_m = solve_radial_1d.Oracle(dict(radial_order=ORDER,
                                           radial_nodes=96, angular_nodes=48,
                                           radius=radius))
    _, _, h_mut, _ = oracle_m.evaluate(mutated)
    m_mut = route_a_corrected.corrected_kinetic_hessian(
        kinetic_stage2.build_cache(mutated, 96, 48))
    w_mut = pencil_float64((np.asarray(h_mut) + np.asarray(h_mut).T) / 2,
                           (m_mut + m_mut.T) / 2)[0]
    n_mut = min(len(w_mut), n_cmp)
    shifts = np.abs(np.sort(np.abs(w_mut[:n_mut])) - wa[:n_mut])
    rel_shifts = shifts / np.maximum(wa[:n_mut], 1e-300)
    out["mutation_max_relative_shift"] = float(rel_shifts.max())
    out["mutation_pass"] = bool(rel_shifts.max() > 1e-7)

    uncertified = [r["mode"] for r in rows if not r["certified_margin_ok"]]
    certified = [r for r in rows if r["certified_margin_ok"] and r["omega"] > 0]
    delta_e = 0.5 * math.fsum(r["omega"] for r in certified)
    sig_rss = 0.5 * math.sqrt(math.fsum(r["sigma_omega"] ** 2 for r in certified))
    sig_lin = 0.5 * math.fsum(r["sigma_omega"] for r in certified)
    # Cross-resolution Delta_E comparison is deliberately NOT computed:
    # frequency crossings reorder sorted positions between quadratures, so
    # position-wise sums differ without implying physical disagreement.
    # This mirrors the pairing-free stance of attempt 0009's preregistration;
    # per-mode ladder gaps still enter every budget via sigma_w.
    cf_note = ("not computed: crossings prevent position-wise "
               "cross-resolution mode pairing")

    out.update({
        "kept_mode_count": int(keep_mask.sum()),
        "compared_modes": n_cmp,
        "certified_mode_count": len(certified),
        "uncertified_modes": uncertified,
        "delta_E": delta_e,
        "sigma_rss_half": sig_rss,
        "sigma_linear_half": sig_lin,
        "omega_min_certified": (min(r["omega"] for r in certified)
                                if certified else None),
        "coarse_fine_deltaE_agreement": cf_note,
        "rows": rows,
    })
    return out


def main():
    started = time.time()
    roots = json.loads((REPO / "proposals/P240-m5-kinetic-axis/"
                        "attempts/0042/largeR-roots.json").read_text())

    rungs = [evaluate_rung(tag, R, e_ref, roots) for tag, R, e_ref in RUNGS]

    checks = []
    for r in rungs:
        checks.append({"name": f"G0_transfer_{r['rung']}",
                       "rel": r["g0_transfer_rel"],
                       "passed": bool(r["g0_pass"])})
        checks.append({"name": f"G3_mutation_{r['rung']}",
                       "max_relative_shift": r["mutation_max_relative_shift"],
                       "passed": bool(r["mutation_pass"])})
        checks.append({"name": f"G5_all_kept_certified_{r['rung']}",
                       "uncertified_modes": r["uncertified_modes"],
                       "passed": len(r["uncertified_modes"]) == 0})
    g0r = abs(rungs[1]["delta_E"] - DE12_COMMITTED)
    checks.append({"name": "G0R_regression_CM5S_011",
                   "abs_diff": g0r,
                   "passed": bool(g0r <= 1e-9)})

    d = {r["rung"]: r["delta_E"] for r in rungs}
    ratio_out = d["R10"] / d["R12"]
    ratio_in = d["R14"] / d["R12"]
    p_fit = math.log(d["R14"] / d["R10"]) / math.log(14.0 / 10.0)

    # Outward order is R10 -> R12 -> R14; "per outward rung" multiplies
    # Delta_E by the band factor on each step away from R10.
    plateau = (abs(ratio_out - 1.0) <= VERDICT_TOL
               and abs(ratio_in - 1.0) <= VERDICT_TOL)
    decay = (d["R10"] > d["R12"] > d["R14"]
             and ratio_out >= 1.0 / (2.0 / 3.0) and ratio_in <= 2.0 / 3.0)
    growth = (d["R10"] < d["R12"] < d["R14"]
              and ratio_out <= 2.0 / 3.0 and ratio_in >= 1.5)

    gates_ok = all(c["passed"] for c in checks)
    if not gates_ok:
        verdict = "INCONCLUSIVE: gate failure -- " + "; ".join(
            c["name"] for c in checks if not c["passed"])
    elif plateau:
        verdict = "PLATEAU_INTRINSIC"
    elif decay:
        verdict = "DECAY_BOX_ARTIFACT"
    elif growth:
        verdict = "GROWTH_VOLUME"
    else:
        verdict = "MIXED"

    report = {
        "attempt": "0010-debox-ladder",
        "preregistration": ("proposal.yaml next_loop (attempt '0010'), "
                            "frozen before execution; mirrored in module "
                            "docstring"),
        "thread_pin": "torch.set_num_threads(1)",
        "budget_substitution": ("cross-family spread replaced by GL "
                                "48x24-vs-96x48 ladder gap at every rung; "
                                "entry jitter unchanged"),
        "rungs": rungs,
        "ratios": {"R10_over_R12": ratio_out, "R14_over_R12": ratio_in},
        "power_law_exponent_p": p_fit,
        "interpretation_bands_frozen": {
            "plateau_tol": VERDICT_TOL,
            "decay_ratio_per_rung": 2.0 / 3.0,
            "growth_ratio_per_rung": 1.5,
        },
        "verdict": verdict,
        "checks": checks,
    }
    print(f"[ladder] D10={d['R10']:.9f} D12={d['R12']:.9f} "
          f"D14={d['R14']:.9f}", flush=True)
    print(f"[ladder] ratios out={ratio_out:.6f} in={ratio_in:.6f} "
          f"p={p_fit:.4f}", flush=True)
    tally = sum(1 for c in checks if c.get("passed"))
    report["tally"] = f"{tally}/{len(checks)} CHECKS PASS"
    report["runtime_seconds"] = round(time.time() - started, 1)
    (HERE / "debox-verdict.json").write_text(json.dumps(report, indent=1))
    print(f"VERDICT: {verdict}", flush=True)
    print(report["tally"], flush=True)


if __name__ == "__main__":
    main()
