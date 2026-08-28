"""P247 attempt 0010: fixed-width existence probe (physical construction).

Question: does the model admit a FIXED-WIDTH stationary clock in a box
larger than the de-boxed family's continuation?  Construction: embed
the R=24 root's physical eigenvalue profiles into a larger box, continue
each row QUADRATICALLY to the boundary (the ansatz-native pinning form,
which keeps every modal function finite at n=1), reproject into the
modal parametrization at order 24, gate on seed fidelity, then relax
with the committed Newton solver.

  CONVERGES NARROW (width ~= R24's, E < family) -> fixed-width branch
      member exists; the isolated de-boxed clock exists.
  SLIDES TO SPREAD (width/E -> family values)   -> the spread branch is
      the only attractor; the width quasi-moduli obstruction is
      confirmed at solver level.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "0001"))
sys.path.insert(0, str(HERE.parent / "0002"))

import debox_common as base  # noqa: F401,E402  (installs P240 paths)
import solve_radial_1d  # noqa: E402

ORDER = 24


def profiles_from_coef(coef: np.ndarray, R: float, r: np.ndarray):
    """Physical eigenvalue profiles (lam_n, lam_p, lam_a) on the equator."""
    n2 = np.clip(r / R, 0.0, None) ** 2
    x = np.clip(2.0 * n2 - 1.0, -1.0, 1.0)
    m0 = chebval(x, coef[0])
    m1 = chebval(x, coef[1])
    m2 = chebval(x, coef[2])
    q = n2 + n2 * (1.0 - n2) * m0
    tangent = (1.0 - n2) * (1.0 / 3.0 + m1)
    delta = n2**2 * (1.0 - n2) * m2
    return tangent + q, tangent + delta, tangent - delta


def embed_and_reproject(coef24: np.ndarray, R_from: float, R_to: float) -> np.ndarray:
    """Embed R_from root into an R_to box; quadratic tail beyond r=R_from."""
    n_fine = np.linspace(1e-4, 1.0 - 1e-6, 4000)
    r_fine = n_fine * R_to
    n2 = n_fine**2
    inner = r_fine <= R_from
    lam_in = profiles_from_coef(coef24, R_from, r_fine[inner])
    lam = [np.empty_like(r_fine) for _ in range(3)]
    kappa = (1.0 - n2) / (1.0 - (R_from / R_to) ** 2)
    vac = np.array([1.0, 0.0, 0.0])  # projector eigenvalues (lam_n, lam_p, lam_a)
    for k in range(3):
        lam[k][inner] = lam_in[k]
        a = float(lam_in[k][-1])  # value at the junction r = R_from
        lam[k][~inner] = vac[k] + (a - vac[k]) * kappa[~inner]
    tangent_t = (lam[1] + lam[2]) / 2.0
    delta_t = (lam[1] - lam[2]) / 2.0
    q_t = lam[0] - tangent_t
    denom02 = n2 * (1.0 - n2)
    m0_t = (q_t - n2) / denom02
    m1_t = tangent_t / (1.0 - n2) - 1.0 / 3.0
    m2_t = delta_t / denom02
    x_fine = 2.0 * n2 - 1.0
    seed = np.zeros((3, ORDER))
    for row, m_t in enumerate([m0_t, m1_t, m2_t]):
        seed[row] = chebfit(x_fine, m_t, ORDER - 1)
    return seed


def width_proxy(coef: np.ndarray, R: float) -> float:
    n = np.linspace(0.02, 1.0, 400)
    lam = profiles_from_coef(coef, R, n * R)
    dev = np.zeros_like(n)
    for lam_row in lam:
        dev += np.abs(lam_row - lam_row[-1])
    half = dev.max() / 2.0
    idx = int(np.argmax(dev))
    beyond = np.where(dev[idx:] < half)[0]
    return float(n[idx + beyond[0]] * R) if len(beyond) else float(R)


def main() -> None:
    started = time.time()
    roots = np.load(HERE.parent / "0003" / "roots24.npz")
    coef24 = np.asarray(roots["R24"], dtype=np.float64).reshape(3, ORDER)
    w24 = width_proxy(coef24, 24.0)

    payload = {"width_R24": w24, "runs": {}}
    for R_box in [26.0, 30.0, 36.0]:
        seed = embed_and_reproject(coef24, 24.0, R_box)
        oracle = solve_radial_1d.Oracle(dict(radial_order=ORDER, radial_nodes=48,
                                             angular_nodes=24, radius=R_box))
        e_seed = float(oracle.evaluate(seed.reshape(-1))[0])
        w_seed = width_proxy(seed, R_box)
        rec = {"E_seed": e_seed, "width_seed": w_seed}
        print(f"[R{R_box:.0f}] E(seed)={e_seed:.4f} width(seed)={w_seed:.3f} "
              f"(R24 physical width {w24:.3f})", flush=True)
        if e_seed > 100.0:
            rec["verdict"] = "SEED_UNFAITHFUL_ABORT"
            payload["runs"][f"R{int(R_box)}"] = rec
            continue
        sol = solve_radial_1d.solve_order(ORDER, seed.reshape(-1),
                                          dict(radial_nodes=48, angular_nodes=24,
                                               radius=R_box))
        values = np.asarray(sol["values"], dtype=np.float64).reshape(3, ORDER)
        energy = float(sol["energy"])
        w_conv = width_proxy(values, R_box)
        rec.update({"E_converged": energy, "width_converged": w_conv,
                    "relgrad": sol.get("rel_grad", sol.get("relgrad"))})
        payload["runs"][f"R{int(R_box)}"] = rec
        print(f"[R{R_box:.0f}] E={energy:.6f} width_conv={w_conv:.3f} "
              f"relgrad={rec['relgrad']}", flush=True)

    payload["minutes_total"] = round((time.time() - started) / 60.0, 2)
    (HERE / "fixed-width-probe.json").write_text(json.dumps(payload, indent=2, default=float))
    print("WROTE fixed-width-probe.json")


if __name__ == "__main__":
    main()
