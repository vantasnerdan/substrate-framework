"""P247 attempt 0006, gate D1': does the symmetric clock source an AXIAL
(parity-odd, mu = cos theta) boost tail?

Attempt 0005 measured the chi-source only with radially symmetric test
fields (even angular channels). The clock's own velocity response is
axial-structured, so the parity argument forbids only the EVEN
(spherically symmetric) tail; the ODD (mu, mu^3, ...) channels were never
probed. This measurement decides the one-clock route:

  gradient of the fixed-J extended energy in the axial chi channel at
  chi = 0 on the committed order-24 R = 24 root, for angular degrees
  (1,) and (1, 3), with jitter cross-checks.

SURVIVES (nonzero above roundoff) -> the symmetric clock sources an axial
massive tail: the localized one-clock construction reopens inside the
single-centered machinery (build D2' Newton). VANISHES -> the one-clock
no-go extends to all channels and D3/route-G stand.
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
JITTER_DRAWS = 3


def source_for(degrees, root_s, jitter_seed: int | None = None):
    flat = torch.zeros((3 + len(degrees), ORDER), dtype=torch.float64)
    values = root_s
    if jitter_seed is not None:
        rng = np.random.default_rng(jitter_seed)
        values = root_s + 1.0e-11 * np.linalg.norm(root_s) * rng.standard_normal(root_s.size) / np.sqrt(root_s.size)
    flat[:3] = torch.tensor(values.reshape(3, ORDER))
    flat.requires_grad_(True)
    total, _comp = c3.extended_static(
        flat,
        radial_order=ORDER,
        radial_nodes=NODES[0],
        angular_nodes=NODES[1],
        radius=24.0,
        chi_angular_degrees=degrees,
    )
    (grad,) = torch.autograd.grad(total, flat)
    return float(grad[3:].norm()), float(grad[3:].abs().max())


def main() -> None:
    started = time.time()
    roots = np.load(HERE.parent / "0003" / "roots24.npz")
    root_s = np.asarray(roots["R24"], dtype=np.float64)

    payload = {"probes": {}}
    for degrees in [(1,), (3,), (1, 3), (0,), (2,)]:
        norm, max_entry = source_for(tuple(degrees), root_s)
        jitter = [
            source_for(tuple(degrees), root_s, jitter_seed=20260828 + draw)[0]
            for draw in range(JITTER_DRAWS)
        ]
        payload["probes"][str(degrees)] = {
            "source_norm": norm,
            "max_abs_entry": max_entry,
            "jitter_norms": jitter,
        }
        print(
            f"[D1' degrees={degrees}] source_norm={norm:.6e} "
            f"max={max_entry:.3e} jitter={['%.2e' % j for j in jitter]}",
            flush=True,
        )

    axial = payload["probes"]["(1, 3)"]
    survives = bool(axial["source_norm"] > 0.0)
    payload["verdict"] = "SURVIVES_AXIAL" if survives else "VANISHES_ALL_CHANNELS"
    payload["routing"] = (
        "D2' axial-tail Newton (one-clock reopens in single-centered machinery)"
        if survives
        else "D3 / route-G stand (one-clock no-go extends to all channels)"
    )
    payload["minutes_total"] = round((time.time() - started) / 60.0, 2)
    (HERE / "d1p-axial-source.json").write_text(json.dumps(payload, indent=2))
    print(f"D1' VERDICT: {payload['verdict']}")
    print("WROTE d1p-axial-source.json")


if __name__ == "__main__":
    main()
