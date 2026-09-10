#!/usr/bin/env python3
"""P253/0117 stage-1 member build (beacon, frozen design.md).

Solves -D*u = r^2 f(u), f = e^-2 (u - c r^2/2 - mu)_+^p on the even half
slice with decaying Dirichlet outer BC, via damped fixed-point. Monitors:
solver residual, increment, circulation reproduction, boundary decay,
crossed h x R refinement with OBSERVED order. Exploratory until design
gates pass; then production record below.
"""

from __future__ import annotations

import math
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np
from scipy.sparse.linalg import spsolve

from skfem import Basis, ElementTriP1, MeshTri, asm
from skfem.helpers import ddot, grad
from skfem.models.poisson import unit_load  # noqa: F401 (documents API family)


P = 6
EPS = 0.15
R_RING = 1.0
KAPPA = 1.0
LLOG = math.log(1.0 / EPS)
C_SPEED = KAPPA * LLOG / (4 * math.pi * R_RING)
MU = 3 * KAPPA * R_RING * LLOG / (8 * math.pi)


def build_mesh(rmax: float, zmax: float, nr: int, nz: int):
    return MeshTri.init_tensor(np.linspace(0.0, rmax, nr + 1),
                               np.linspace(0.0, zmax, nz + 1))



def dirichlet_dofs(basis, mesh, rmax, zmax, tol=1e-9):
    locs = basis.doflocs
    mask = (locs[0] > rmax - tol) | (locs[1] > zmax - tol)
    return np.where(mask)[0]

def solve_stage1(rmax=6.0, zmax=3.0, nr=80, nz=40, itmax=300,
                 tol=1e-8, alpha=0.1, verbose=True, mu=MU, u0=None,
                 res_every=10):
    mesh = build_mesh(rmax, zmax, nr, nz)
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0]
    zn = basis.doflocs[1]
    u = (u0 if u0 is not None
         else 1.5 * np.exp(-(((rn - R_RING) ** 2 + zn**2) / 0.25**2)))
    dd = dirichlet_dofs(basis, mesh, rmax, zmax)
    free = np.setdiff1d(np.arange(basis.N), dd)
    from skfem import BilinearForm
    from skfem.helpers import dot, grad

    @BilinearForm
    def stiff(u, v, w):
        return w.x[0] * dot(grad(u), grad(v))

    A = asm(stiff, basis).tocsr()
    Af = A[free][:, free].tocsc()
    hist = []
    for it in range(itmax):
        src = u - C_SPEED * rn**2 / 2 - mu
        f = EPS**-2 * np.where(src > 0, src, 0.0) ** P
        from skfem import LinearForm
        finterp = basis.interpolator(f)

        @LinearForm
        def load(v, w):
            return (w.x[0] ** 3) * finterp(w.x) * v

        b = asm(load, basis)
        unew = np.zeros_like(u)
        unew[free] = spsolve(Af, b[free])
        unew[dd] = 0.0
        du = np.max(np.abs(unew - u))
        if not np.isfinite(du):
            print(f"it={it} DIVERGED (non-finite step); mu={mu}", flush=True)
            hist.append((np.inf, np.inf))
            break
        u = (1 - alpha) * u + alpha * unew
        # true residual of damped iterate (periodic cadence for speed)
        if it % res_every == 0 or du < tol:
            src2 = u - C_SPEED * rn**2 / 2 - mu
            with np.errstate(over="ignore", invalid="ignore"):
                f2 = EPS**-2 * np.where(src2 > 0, src2, 0.0) ** P
            finterp2 = basis.interpolator(f2)

            @LinearForm
            def load2(v, w):
                return (w.x[0] ** 3) * finterp2(w.x) * v

            res = A @ u - asm(load2, basis)
            rnorm = float(np.sqrt(res[free] @ res[free] / max(1, free.size)))
        else:
            rnorm = hist[-1][1] if hist else np.inf
        hist.append((du, rnorm))
        if verbose and (it % 25 == 0 or du < tol):
            print(f"it={it} maxdu={du:.3e} res={rnorm:.3e}", flush=True)
        if du < tol:
            break
    # monitors
    src = u - C_SPEED * rn**2 / 2 - mu
    from skfem import BilinearForm as BF
    from skfem.helpers import inner

    @BF
    def massr(u, v, w):
        return w.x[0] * inner(u, v)

    M = asm(massr, basis)
    kappa_hat = float(math.fsum((M @ f).tolist()))
    umax = float(np.max(u))
    rcore = float(rn[np.argmax(u)])
    bnd = float(np.max(np.abs(u[dd]))) if dd.size else 0.0
    return {"u": u, "mesh": mesh, "basis": basis, "hist": hist,
            "kappa_hat": kappa_hat, "umax": umax, "rcore": rcore,
            "bnd": bnd, "rmax": rmax, "zmax": zmax, "nr": nr, "nz": nz,
            "res": hist[-1][1], "du": hist[-1][0], "iters": len(hist)}

def main(argv=None) -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--bordered", action="store_true")
    ap.add_argument("--newton", action="store_true")
    ap.add_argument("--nr", type=int, default=80)
    ap.add_argument("--nz", type=int, default=40)
    ap.add_argument("--reg", type=float, default=0.0)
    args = ap.parse_args(argv)
    if args.bordered:
        print(f"MESH nr={args.nr} nz={args.nz} box=6x3 p-chain", flush=True)
        uu, mm, cc = None, None, None
        out = None
        for pp in (2, 3, 4, 5, 6):
            out = solve_bordered(nr=args.nr, nz=args.nz, reg=1e-3,
                                 p_pw=pp, u0=uu, mu0=mm, c0=cc)
            print(f"P {pp}: res={out['res']:.3e} "
                  f"iters={out['iters']} umax={out['u'].max():.4f} "
                  f"kap={out['kap']:.4f} rbar={out['rbar']:.4f} "
                  f"iz={out['iz']:.4f} mu={out['mu']:.4f} c={out['c']:.4f}",
                  flush=True)
            uu, mm, cc = out["u"], out["mu"], out["c"]
            np.savez("proposals/P253-euler-particle-mechanisms/attempts/0117-beacon-member/"
                     f"member-p{pp}-exploratory.npz",
                     u=out["u"], mu=out["mu"], c=out["c"],
                     kap=out["kap"], rbar=out["rbar"], iz=out["iz"],
                     res=out["res"])
        return
    if args.newton:
        print(f"MESH nr={args.nr} nz={args.nz} box=6x3", flush=True)
        out = solve_newton(nr=args.nr, nz=args.nz, reg=args.reg)
        print(f"NEWTON res={out['res']:.3e} iters={out['iters']} "
              f"umax={out['u'].max():.4f} r0={out['r0']:.3e} "
              f"kappa={out['kappa_hat']:.4f} reg={args.reg}", flush=True)
        np.savez("proposals/P253-euler-particle-mechanisms/attempts/0117-beacon-member/"
                 "member-newton-exploratory.npz",
                 u=out["u"], res=out["res"])
        return
    u0 = None
    out = None
    for mu in (0.8, 0.55, 0.35, MU):
        out = solve_stage1(mu=mu, u0=u0)
        print(f"mu={mu:.4f} kappa_hat={out['kappa_hat']:.6f} "
              f"umax={out['umax']:.4f} res={out['res']:.3e} "
              f"iters={out['iters']}", flush=True)
        if not np.isfinite(out["du"]):
            break
        u0 = out["u"]
    print(f"kappa_hat={out['kappa_hat']:.6f} (target 1)")
    print(f"umax={out['umax']:.4f} rcore={out['rcore']:.4f} bnd={out['bnd']:.2e}")
    print(f"res={out['res']:.3e} iters={out['iters']}")
    np.savez("proposals/P253-euler-particle-mechanisms/attempts/0117-beacon-member/"
             "member-stage1-exploratory.npz",
             u=out["u"], kappa_hat=out["kappa_hat"], res=out["res"])

def solve_newton(rmax=6.0, zmax=3.0, nr=80, nz=40, itmax=30, tol=1e-9,
                 verbose=True, mu=MU, u0=None, reg=0.0):
    from skfem import BilinearForm, LinearForm
    from skfem.helpers import dot, grad, inner
    mesh = build_mesh(rmax, zmax, nr, nz)
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0]
    zn = basis.doflocs[1]
    u = (u0 if u0 is not None
         else 1.5 * np.exp(-(((rn - R_RING) ** 2 + zn**2) / 0.25**2)))
    dd = dirichlet_dofs(basis, mesh, rmax, zmax)
    free = np.setdiff1d(np.arange(basis.N), dd)

    @BilinearForm
    def stiff(u, v, w):
        return w.x[0] * dot(grad(u), grad(v))

    A = asm(stiff, basis).tocsr()

    def smax(src):
        if reg <= 0:
            return np.where(src > 0, src, 0.0), (src > 0).astype(float)
        root = np.sqrt(src**2 + reg**2)
        return (src + root) / 2, 0.5 * (1 + src / root)

    def residual(u):
        src = u - C_SPEED * rn**2 / 2 - mu
        with np.errstate(over="ignore", invalid="ignore"):
            s, _ = smax(src)
            f = EPS**-2 * s ** P
        fi = basis.interpolator(f)
        @LinearForm
        def load(v, w):
            return (w.x[0] ** 3) * fi(w.x) * v

        r = A @ u - asm(load, basis)
        r[dd] = 0.0
        return r, f, src

    r, f, src = residual(u)
    n0 = float(np.sqrt(r[free] @ r[free] / free.size))
    nrm = n0
    for it in range(itmax):
        s, ds = smax(src)
        jac_w = 6 * EPS**-2 * s ** 5 * ds
        ji = basis.interpolator(jac_w)

        @BilinearForm
        def jacform(a, b, w):
            return (w.x[0] * dot(grad(a), grad(b))
                    - (w.x[0] ** 3) * ji(w.x) * inner(a, b))

        J = asm(jacform, basis).tocsr()
        Jf = J[free][:, free].tocsc()
        du = np.zeros_like(u)
        du[free] = spsolve(Jf, -r[free])
        step = 1.0
        for _ in range(30):
            u_try = u + step * du
            u_try[dd] = 0.0
            r_try, _, _ = residual(u_try)
            n_try = float(np.sqrt(r_try[free] @ r_try[free] / free.size))
            n_cur = float(np.sqrt(r[free] @ r[free] / free.size))
            if n_try < n_cur:
                break
            step *= 0.5
        u = u + step * du
        u[dd] = 0.0
        r, f, src = residual(u)
        nrm = float(np.sqrt(r[free] @ r[free] / free.size))
        if verbose:
            print(f"newt it={it} res={nrm:.3e} step={step:.2f} "
                  f"umax={u.max():.4f} active={(src > 0).sum()}", flush=True)
        if nrm < tol:
            break
    from skfem import BilinearForm as BF2
    from skfem.helpers import inner as inner2

    @BF2
    def massr2(u, v, w):
        return w.x[0] * inner2(u, v)

    import math as _math
    kappa_hat = float(_math.fsum((asm(massr2, basis) @ f).tolist()))
    return {"u": u, "res": nrm, "iters": it + 1, "r0": n0,
            "mesh": mesh, "basis": basis, "mu": mu,
            "kappa_hat": kappa_hat, "reg": reg}

def solve_bordered(rmax=6.0, zmax=3.0, nr=40, nz=20, itmax=20, tol=1e-9,
                   verbose=True, reg=1e-3, rbar_target=None,
                   u0=None, mu0=None, c0=None, p_pw=None):
    """Bordered Newton: unknowns (u, mu, c), rows (kappa, mean-radius).
    0080 BR-border numerically: (kappa, I_z) are near-parallel at frozen
    geometry (r~=r^2 over a thin core at R=1, measured det ratio 3.23 vs
    3.25); mean radius breaks the degeneracy and frees radial motion.
    I_z becomes a post-hoc output check against the leading-jet pi."""
    import math as _math
    from skfem import BilinearForm, LinearForm
    from skfem.helpers import dot, grad, inner
    if p_pw is None:
        p_pw = P
    if rbar_target is None:
        rbar_target = 1.0
    mesh = build_mesh(rmax, zmax, nr, nz)
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0]
    zn = basis.doflocs[1]
    u = (u0 if u0 is not None
         else 1.5 * np.exp(-(((rn - R_RING) ** 2 + zn**2) / 0.25**2)))
    mu = MU if mu0 is None else mu0
    c = C_SPEED if c0 is None else c0
    dd = dirichlet_dofs(basis, mesh, rmax, zmax)
    free = np.setdiff1d(np.arange(basis.N), dd)

    @BilinearForm
    def stiff(a, b, w):
        return w.x[0] * dot(grad(a), grad(b))

    @BilinearForm
    def massr(a, b, w):
        return w.x[0] * inner(a, b)

    @BilinearForm
    def massr2(a, b, w):
        return (w.x[0] ** 2) * inner(a, b)

    @BilinearForm
    def massr3(a, b, w):
        return (w.x[0] ** 3) * inner(a, b)

    A = asm(stiff, basis).tocsr()
    Mr = asm(massr, basis)
    Mr2 = asm(massr2, basis)
    M3 = asm(massr3, basis)
    m1 = np.asarray(Mr.sum(axis=0)).ravel()
    m2 = np.asarray(Mr2.sum(axis=0)).ravel()

    def smax(src):
        root = np.sqrt(src**2 + reg**2)
        return (src + root) / 2, 0.5 * (1 + src / root)

    def full(u, mu, c):
        src = u - c * rn**2 / 2 - mu
        with np.errstate(over="ignore", invalid="ignore"):
            s, ds = smax(src)
            f = EPS**-2 * s ** p_pw
        fi = basis.interpolator(f)

        @LinearForm
        def load(v, w):
            return (w.x[0] ** 3) * fi(w.x) * v

        F = A @ u - asm(load, basis)
        F[dd] = 0.0
        kap = float(_math.fsum((Mr @ f).tolist()))
        num = float(_math.fsum((Mr @ (f * rn)).tolist()))
        rbar = num / kap if kap > 0 else float("nan")
        iz = float(_math.pi * _math.fsum((Mr2 @ f).tolist()))
        return F, f, s, ds, kap, rbar, iz

    for it in range(itmax):
        F, f, s, ds, kap, rbar, iz = full(u, mu, c)
        rows = np.array([kap - 1.0, rbar - rbar_target])
        nrm = float(np.sqrt(F[free] @ F[free] / free.size + rows @ rows))
        if verbose:
            print(f"bord it={it} res={nrm:.3e} umax={u.max():.4f} "
                  f"kap={kap:.4f} rbar={rbar:.4f} iz={iz:.4f} "
                  f"mu={mu:.4f} c={c:.4f}", flush=True)
        jf = p_pw * EPS**-2 * s ** (p_pw - 1) * ds
        ji = basis.interpolator(jf)

        @BilinearForm
        def jacform(a, b, w):
            return (w.x[0] * dot(grad(a), grad(b))
                    - (w.x[0] ** 3) * ji(w.x) * inner(a, b))

        J = asm(jacform, basis).tocsr()
        B = np.zeros((basis.N, 2))
        B[:, 0] = (M3 @ jf)
        B[:, 1] = (M3 @ (jf * rn**2 / 2))
        Jf = J[free][:, free].tocsc()
        t0 = np.zeros(basis.N)
        t0[free] = spsolve(Jf, -F[free])
        T = np.zeros((basis.N, 2))
        T[free, :] = spsolve(Jf, -B[free, :])
        dN = m1 * rn * jf
        dD = m1 * jf
        # quotient rule: d(N/D) = (dN*D - N*dD)/D^2, N = kap*rbar
        Crow = np.stack([dD, (dN * kap - kap * rbar * dD) / max(kap, 1e-300)**2])
        S = Crow @ T
        try:
            dp = np.linalg.solve(S, -rows - Crow @ t0)
        except np.linalg.LinAlgError:
            print("bord: singular Schur complement", flush=True)
            break
        t0n = float(np.max(np.abs(t0)))
        dpn = float(np.max(np.abs(dp)))
        step = 1.0
        cur = nrm
        improved = False
        ray = []
        for _ in range(64):
            u_t = u + step * (t0 + T @ dp)
            u_t[dd] = 0.0
            F_t, _, _, _, k_t, rb_t, _ = full(u_t, mu + step * dp[0],
                                             c + step * dp[1])
            r_t = np.array([k_t - 1.0, rb_t - rbar_target])
            n_t = float(np.sqrt(F_t[free] @ F_t[free] / free.size
                                + r_t @ r_t))
            ray.append(round(n_t, 6))
            if n_t < cur:
                improved = True
                break
            step *= 0.5
        if not improved:
            print(f"bord it={it} STALL (|t0|={t0n:.2e} |dp|={dpn:.2e} "
                  f"detS={np.linalg.det(S):.2e} ray={ray})", flush=True)
            break
        u = u + step * (t0 + T @ dp)
        u[dd] = 0.0
        mu, c = mu + step * dp[0], c + step * dp[1]
    F, f, s, ds, kap, rbar, iz = full(u, mu, c)
    rows = np.array([kap - 1.0, rbar - rbar_target])
    nrm = float(np.sqrt(F[free] @ F[free] / free.size + rows @ rows))
    return {"u": u, "mu": mu, "c": c, "res": nrm, "iters": it + 1,
            "kap": kap, "rbar": rbar, "iz": iz, "mesh": mesh,
            "basis": basis, "f": f, "reg": reg}
if __name__ == "__main__":
    main()
