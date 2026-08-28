"""P247 attempt 0006 continuation, gate F1-routing: is the extended
functional even under global chi -> -chi?

If E[+chi] == E[-chi] exactly (machine precision) on GENERIC S
backgrounds, the linear chi channel is dead by an exact Z2 symmetry for
EVERY configuration, and the pair interaction law must come from the
SECOND-order chi channel (still massive, still e^{-Lambda d}-type).
If instead E[+chi] != E[-chi] on generic backgrounds, the linear source
is alive off the symmetric class and the dipole-current pairing proceeds
as preregistered.

Discriminator: random modal chi profiles on (a) the committed symmetric
root and (b) randomized S backgrounds. Multiple draws, several angular
channels.
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
NODES = (72, 36)


def energy(root_s, chi_rows, degrees, radius=24.0):
    flat = torch.zeros((3 + len(degrees), ORDER), dtype=torch.float64)
    flat[:3] = torch.tensor(root_s.reshape(3, ORDER))
    flat[3:] = torch.tensor(chi_rows.reshape(len(degrees), ORDER))
    total, _ = c3.extended_static(
        flat,
        radial_order=ORDER,
        radial_nodes=NODES[0],
        angular_nodes=NODES[1],
        radius=radius,
        chi_angular_degrees=degrees,
    )
    return float(total)


def main() -> None:
    started = time.time()
    roots = np.load(HERE.parent / "0003" / "roots24.npz")
    root_s = np.asarray(roots["R24"], dtype=np.float64)
    rng = np.random.default_rng(20260828)

    payload = {"cases": {}}
    for label, s_background in [
        ("symmetric_root", root_s),
        ("randomized_s", root_s + 0.05 * np.linalg.norm(root_s) * rng.standard_normal(root_s.size) / np.sqrt(root_s.size)),
    ]:
        for degrees in [(0,), (1,), (1, 3)]:
            gaps = []
            for draw in range(4):
                c = rng.standard_normal((len(degrees), ORDER))
                e_plus = energy(s_background, c, degrees)
                e_minus = energy(s_background, -c, degrees)
                gaps.append(abs(e_plus - e_minus))
            payload["cases"][f"{label}_deg{degrees}"] = {
                "gaps": gaps,
                "max_gap": max(gaps),
                "scale": abs(energy(s_background, np.zeros((len(degrees), ORDER)), degrees)),
            }
            print(f"[{label} deg={degrees}] max|E(+c)-E(-c)| = {max(gaps):.3e}", flush=True)

    generic_max = max(
        v["max_gap"] / max(v["scale"], 1e-300)
        for k, v in payload["cases"].items()
        if k.startswith("randomized")
    )
    payload["verdict"] = (
        "EVEN_FUNCTIONAL_Z2_EXACT" if generic_max < 1e-12 else "LINEAR_CHANNEL_ALIVE_OFF_SYMMETRIC_CLASS"
    )
    payload["routing"] = {
        "EVEN_FUNCTIONAL_Z2_EXACT": "F1 targets the SECOND-order chi channel (quadratic-density overlap pairing)",
        "LINEAR_CHANNEL_ALIVE_OFF_SYMMETRIC_CLASS": "F1 proceeds with the preregistered dipole-current pairing",
    }[payload["verdict"]]
    payload["minutes_total"] = round((time.time() - started) / 60.0, 2)
    (HERE / "f1-routing-z2.json").write_text(json.dumps(payload, indent=2))
    print(f"F1-ROUTING VERDICT: {payload['verdict']}")
    print("WROTE f1-routing-z2.json")


if __name__ == "__main__":
    main()
