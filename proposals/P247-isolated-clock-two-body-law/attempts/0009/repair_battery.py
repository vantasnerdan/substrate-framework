"""P247 attempt 0009: chi-sector kinetic repair and full re-verification.

Defect (attempts 0005/0008): the extended functional's sigma term used
-eta(du, du), which for a static boost profile is
-(|grad chi|^2 + sinh^2 chi |grad n_hat|^2) <= 0 - a ghost kinetic
norm, contradicting the accepted census (C-M5S-001: the boost orbit is
a POSITIVE propagating species, metric entry 1/16).

Repair: sigma_density = +KAPPA * eta(du, du) - the positive Legendre
density. It reduces to the census linearization (positive |grad chi|^2)
at small chi and stays positive at finite chi.

Battery:
  R1  S-sector regression: at chi = 0 the repaired functional's total
      energy equals the original's exactly (sigma vanishes at chi = 0,
      so the committed roots' stationarity and energies are untouched).
  R2  sigma positivity at finite chi (several chi amplitudes/channels).
  R3  Z2 evenness on the REPAIRED functional (root + randomized S).
  R4  full fixed-R Hessian spectrum at the R=24 root on the repaired
      functional: is the ghost gone?  Morse index?
  R5  all-channel linear chi source on the repaired functional.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "0004"))
sys.path.insert(0, str(HERE.parent / "0002"))
sys.path.insert(0, str(HERE.parent / "0001"))

import c3_functional as c3  # noqa: E402

ORDER = 24

# Build the repaired functional: same source, sigma sign flipped.
src = (HERE.parent / "0004" / "c3_functional.py").read_text()
marker = """    sigma_density = KAPPA * (
        -torch.sum(du_r * (du_r @ ETA), dim=-1)
        - torch.sum(du_theta * (du_theta @ ETA), dim=-1)"""
repaired = """    sigma_density = KAPPA * (
        torch.sum(du_r * (du_r @ ETA), dim=-1)
        + torch.sum(du_theta * (du_theta @ ETA), dim=-1)"""
assert src.count(marker) == 1, "sigma block not found"
src2 = src.replace(marker, repaired)
src2 = src2.replace("def extended_static(", "def extended_static_repaired(", 1)
src2 = src2.replace("__file__", repr(str(HERE.parent / "0004" / "c3_functional.py")))
ns = {"__name__": "repaired_c3"}
exec(compile(src2, "repaired", "exec"), ns)
ext_rep = ns["extended_static_repaired"]


def run(flat, radius, nodes=(48, 24), raw=False):
    fn = ext_rep if raw else c3.extended_static
    total, comp = fn(
        flat, radial_order=ORDER, radial_nodes=nodes[0],
        angular_nodes=nodes[1], radius=radius, chi_angular_degrees=(0,),
    )
    return total, comp


def r1_regression(roots) -> dict:
    out = {}
    for name in ["R20", "R24", "R30"]:
        R = float(name[1:])
        flat = torch.zeros((4, ORDER), dtype=torch.float64)
        flat[:3] = torch.tensor(np.asarray(roots[name], dtype=np.float64).reshape(3, ORDER))
        t0, _ = c3.extended_static(flat, radial_order=ORDER, radial_nodes=48,
                                   angular_nodes=24, radius=R, chi_angular_degrees=(0,))
        t1, _ = ext_rep(flat, radial_order=ORDER, radial_nodes=48,
                        angular_nodes=24, radius=R, chi_angular_degrees=(0,))
        out[name] = {"orig": float(t0), "repaired": float(t1),
                     "abs_diff": abs(float(t0) - float(t1))}
        print(f"[R1 {name}] orig={float(t0):.10f} repaired={float(t1):.10f}", flush=True)
    return out


def r2_sigma_positivity(roots) -> dict:
    root = np.asarray(roots["R24"], dtype=np.float64).reshape(3, ORDER)
    out = {}
    for deg, label in [((0,), "even-radial"), ((1,), "odd-axial")]:
        for amp in [0.05, 0.2, 0.5]:
            flat = torch.zeros((4, ORDER), dtype=torch.float64)
            flat[:3] = torch.tensor(root)
            flat[3] = torch.tensor(amp * np.ones(ORDER) / np.sqrt(ORDER))
            _t, comp = ext_rep(flat, radial_order=ORDER, radial_nodes=48,
                               angular_nodes=24, radius=24.0, chi_angular_degrees=deg)
            s = float(comp["sigma"])
            out[f"{label}_amp{amp}"] = s
            print(f"[R2 {label} amp={amp}] sigma = {s:.6e} {'OK' if s > 0 else 'NEGATIVE'}", flush=True)
    return out


def r3_z2(roots) -> dict:
    root = np.asarray(roots["R24"], dtype=np.float64).reshape(3, ORDER)
    rng = np.random.default_rng(20260828)
    gaps = {}
    for label, s_bg in [("root", root),
                        ("randomized", root + 0.05 * np.linalg.norm(root) * rng.standard_normal(root.size).reshape(root.shape) / np.sqrt(root.size))]:
        for deg in [(0,), (1,), (1, 3)]:
            worst = 0.0
            for _ in range(3):
                c = rng.standard_normal((len(deg), ORDER))
                e = {}
                for sgn in [1.0, -1.0]:
                    flat = torch.zeros((3 + len(deg), ORDER), dtype=torch.float64)
                    flat[:3] = torch.tensor(s_bg.reshape(3, ORDER))
                    flat[3:] = torch.tensor(sgn * c)
                    t, _ = ext_rep(flat, radial_order=ORDER, radial_nodes=48,
                                   angular_nodes=24, radius=24.0, chi_angular_degrees=deg)
                    e[sgn] = float(t)
                worst = max(worst, abs(e[1.0] - e[-1.0]))
            gaps[f"{label}_deg{deg}"] = worst
            print(f"[R3 {label} deg={deg}] max|E(+c)-E(-c)| = {worst:.3e}", flush=True)
    return gaps


def r4_hessian(roots) -> dict:
    root = np.asarray(roots["R24"], dtype=np.float64)
    flat = torch.zeros((4, ORDER), dtype=torch.float64)
    flat[:3] = torch.tensor(root.reshape(3, ORDER))
    flat.requires_grad_(True)
    total, comp = ext_rep(flat, radial_order=ORDER, radial_nodes=48,
                          angular_nodes=24, radius=24.0, chi_angular_degrees=(0,))
    total = total + comp["fixed_j"]
    g = torch.autograd.grad(total, flat, create_graph=True)[0].reshape(-1)
    rows = [torch.autograd.grad(g[i], flat, retain_graph=True)[0].reshape(-1)
            for i in range(g.numel())]
    H = torch.stack(rows).detach().numpy()
    H = (H + H.T) / 2
    ev = np.linalg.eigvalsh(H)
    neg = int((ev < 0).sum())
    print(f"[R4] min={ev[0]:.6e} num_negative={neg}/96 bottom6={np.round(ev[:6], 6)}", flush=True)
    return {"min": float(ev[0]), "num_negative": neg, "bottom6": ev[:6].tolist()}


def r5_source(roots) -> dict:
    root = np.asarray(roots["R24"], dtype=np.float64).reshape(3, ORDER)
    out = {}
    for deg in [(0,), (1,), (1, 3)]:
        flat = torch.zeros((3 + len(deg), ORDER), dtype=torch.float64)
        flat[:3] = torch.tensor(root)
        flat.requires_grad_(True)
        t, _ = ext_rep(flat, radial_order=ORDER, radial_nodes=48,
                       angular_nodes=24, radius=24.0, chi_angular_degrees=deg)
        (g,) = torch.autograd.grad(t, flat)
        out[str(deg)] = float(g[3:].norm())
        print(f"[R5 deg={deg}] |dE/dchi| = {float(g[3:].norm()):.6e}", flush=True)
    return out


def main() -> None:
    started = time.time()
    roots = np.load(HERE.parent / "0003" / "roots24.npz")
    payload = {}
    payload["r1_regression"] = r1_regression(roots)
    payload["r2_sigma_positivity"] = r2_sigma_positivity(roots)
    payload["r3_z2_evenness"] = r3_z2(roots)
    payload["r4_hessian"] = r4_hessian(roots)
    payload["r5_source"] = r5_source(roots)
    payload["minutes_total"] = round((time.time() - started) / 60.0, 2)
    (HERE / "repair-battery.json").write_text(json.dumps(payload, indent=2, default=float))
    print("WROTE repair-battery.json")


if __name__ == "__main__":
    main()
