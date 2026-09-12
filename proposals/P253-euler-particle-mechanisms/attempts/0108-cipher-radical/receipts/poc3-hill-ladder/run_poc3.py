"""PoC-3 archive: horn-1 (B=eps*w) Hill-ball material flux + H_c + m* ladder.

Reproduces receipts/poc3-hill-ladder verdicts (final fix3 generation; see README
for the append-only debug trail B1-B3). Run: python3 run_poc3.py
Scope: frozen analytic Hill background (self-consistency OUT); horn-1 only.
"""
import time
import numpy as np

V, A = 1.0, 1.0


def uvw(sig, z, Vv):
    return 1.5 * Vv * z * sig / (A * A), 1.5 * Vv * (1 - (2 * sig * sig + z * z) / (A * A))


def w_of(s, Vv):
    return 7.5 * Vv * s / (A * A)


def signed_areas(S, Z):
    out = np.zeros((S.shape[0] - 1, S.shape[1] - 1))
    for i in range(out.shape[0]):
        for j in range(out.shape[1]):
            x = [S[i, j], S[i + 1, j], S[i + 1, j + 1], S[i, j + 1]]
            y = [Z[i, j], Z[i + 1, j], Z[i + 1, j + 1], Z[i, j + 1]]
            out[i, j] = 0.5 * sum(x[k] * y[(k + 1) % 4] - x[(k + 1) % 4] * y[k] for k in range(4))
    return out


def cell_center(B):
    return 0.25 * (B[:-1, :-1] + B[1:, :-1] + B[:-1, 1:] + B[1:, 1:])


def run_flux(Vv, eps, nr=24, nph=48, dt=0.001):
    rr = np.linspace(0.04, 0.5, nr)
    ph = np.linspace(-np.pi / 2, np.pi / 2, nph)
    RR, PH = np.meshgrid(rr, ph)
    S0 = (RR * np.cos(PH)) * 0.9
    Z0 = (RR * np.sin(PH)) * 0.9
    B = eps * w_of(S0, Vv)
    S1, Z1, B1 = S0.copy(), Z0.copy(), B.copy()

    def vel(s, z):
        return uvw(s, z, Vv)

    Phi0 = np.sum(cell_center(B) * signed_areas(S0, Z0))
    T = 1.0 / Vv
    for _ in range(int(T / dt)):
        u1, v1 = vel(S1, Z1)
        g1 = np.where(S1 > 1e-9, u1 / S1, 0.0)
        Sm = S1 + 0.5 * dt * u1
        Zm = Z1 + 0.5 * dt * v1
        Bm = B1 + 0.5 * dt * B1 * g1
        u2, v2 = vel(Sm, Zm)
        g2 = np.where(Sm > 1e-9, u2 / Sm, 0.0)
        S1 += dt * u2
        Z1 += dt * v2
        B1 += dt * Bm * g2
    A1 = signed_areas(S1, Z1)
    Phi1 = np.sum(cell_center(B1) * A1)
    return Phi0, Phi1, A1.min()


def main():
    t0 = time.time()
    th = np.linspace(0.05, np.pi - 0.05, 9)
    sg_, zz_ = A * np.sin(th), A * np.cos(th)
    usi, uzi = uvw(sg_, zz_, V)
    uti = usi * np.cos(th) - uzi * np.sin(th)
    print("boundary: max|u_interior_t - 1.5V sin| =", np.max(np.abs(uti - 1.5 * V * np.sin(th))))
    for Vv in [0.5, 1.0, 2.0]:
        P0, P1, q = run_flux(Vv, 0.05)
        print(f"Vv={Vv}: Phi0={P0:.6f} Phi1={P1:.6f} "
              f"drift={abs(P1-P0)/abs(P0):.2e} Phi/(eps*V)={P1/0.05/Vv:.4f}")
    print("H_c identically 0 for swirl-free Hill (u poloidal, w toroidal)")
    for Gc in [0.5, 1.0, 2.0]:
        R, a0 = 1.0, 0.05
        L = np.log(8 * R / a0)
        Iv = np.pi * Gc * R ** 2
        Vv = Gc / (4 * np.pi * R) * (L - 0.25)
        print(f"Gamma={Gc}: m*=I/V={Iv/Vv:.4f}")
    print(f"elapsed {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
