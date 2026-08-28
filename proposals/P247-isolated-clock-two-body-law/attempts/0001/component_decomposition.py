"""Attempt-0001 addendum: x-space component decomposition of the clean ladder roots.

E(R)[c] = R^3 V[c] + (C[c] + Phi[c])/R, Phi = 1/(4 I). V, C, I are radius-free
functionals of the root coefficients, so this isolates which scaling class
carries the de-boxed energy growth (the competing scaling classes preregistered
in P240 attempt 0042).
"""

import json
import sys

import numpy as np
import torch

sys.path.insert(0, ".")
from debox_common import ORDER  # noqa: E402  (installs attempt sys.path entries)

import xspace_energy  # noqa: E402

d = json.load(open("clean-ladder.json"))
print(f"{'R':>5} {'E_total':>12} {'V*R^3':>12} {'C/R':>10} {'Phi/R':>10} {'sum':>12} {'V':>12} {'I*R':>10}")
for row in d["rungs"]:
    if not row.get("accepted"):
        continue
    radius = row["radius"]
    flat = torch.tensor(np.asarray(row["values"]), dtype=torch.float64, requires_grad=False)
    curvature, potential, inertia = xspace_energy.xspace_components(flat, ORDER, 48, 16)
    phi = float(1.0 / (4.0 * inertia))
    r3v = float(radius**3 * potential)
    cover = float(curvature / radius)
    phir = phi / radius
    total = r3v + cover + phir
    print(
        f"{radius:5.0f} {row['energy_total']:12.6f} {r3v:12.6f} {cover:10.6f} {phir:10.6f} "
        f"{total:12.6f} {float(potential):12.6e} {float(inertia * radius):10.6f}"
    )
