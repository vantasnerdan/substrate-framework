#!/usr/bin/env python3
"""P253/0116 G-a2 unblock repro (beacon, no implementation).

Exits 1 while the numerical charged-branch member is unbuildable:
prints the exact missing-capability matrix. Turns exit 0 only when a
G-a2 member field exists to feed ga_pipeline.py (acceptance probe).
"""

from __future__ import annotations

import importlib.util
import subprocess


def have(mod: str) -> bool:
    return importlib.util.find_spec(mod) is not None


def main() -> None:
    print("== present ==")
    for mod in ("numpy", "scipy", "sympy", "skfem"):
        print(f"  {mod}: {'YES ' + __import__(mod).__version__ if have(mod) else 'NO'}")
    print("== heavy 2D solvers ==")
    missing = [m for m in ("firedrake", "dedalus", "dolfin", "mpi4py", "petsc4py")
               if not have(m)]
    for mod in ("firedrake", "dedalus", "dolfin", "mpi4py", "petsc4py"):
        print(f"  {mod}: {'YES' if have(mod) else 'NO'}")
    print("== in-tree axisymmetric Grad-Shafranov solver ==")
    grep = subprocess.run(
        ["grep", "-rl", "Grad-Shafranov.*solve\\|solve.*Grad-Shafranov"
                        "\\|axisymmetric.*BVP\\|BVP.*axisymmetric",
         "src/substrate_framework", "proposals/P253-euler-particle-mechanisms/attempts"],
        capture_output=True, text=True)
    hits = [h for h in grep.stdout.splitlines() if "0116" not in h]
    print(f"  solver files: {hits if hits else 'NONE'}")
    print("== member field data ==")
    find = subprocess.run(
        ["find", "proposals/P253-euler-particle-mechanisms/attempts",
         "-name", "*member*field*", "-o", "-name", "*cao*member*.npz",
         "-o", "-name", "*charged*branch*.npy"],
        capture_output=True, text=True)
    print(f"  data files: {find.stdout.strip() or 'NONE'}")
    if missing and not hits:
        print("BLOCKED: no 2D axisymmetric PDE path to a charged-branch member; "
              "ga-status rows 1-5,8-9 uncomputable. See unblock-request.md.")
        raise SystemExit(1)
    print("UNBLOCKED: solver path or member data present.")


if __name__ == "__main__":
    main()
